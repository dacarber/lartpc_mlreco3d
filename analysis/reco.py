import numpy as np
import ROOT
from sklearn.decomposition import PCA
from scipy.optimize import curve_fit
from scipy.spatial.distance import cdist
import warnings

def load_range_reco(path, particle_type='muon', kinetic_energy=True):
    """
    Return a function maps the residual range of a track to the kinetic
    energy of the track. The mapping is based on the Bethe-Bloch formula
    and stored per particle type in TGraph objects. The TGraph::Eval
    function is used to perform the interpolation.

    Parameters
    ----------
    particle_type: A string with the particle name.
    kinetic_energy: If true (false), return the kinetic energy (momentum)
    
    Returns
    -------
    The kinetic energy or momentum according to Bethe-Bloch.
    """
    output_var = ('_RRtoT' if kinetic_energy else '_RRtodEdx')
    if particle_type in ['muon', 'pion', 'kaon', 'proton']:
        input_file = ROOT.TFile.Open(path, 'read')
        graph = input_file.Get(f'{particle_type}{output_var}')
        return np.vectorize(graph.Eval)
    else:
        print(f'Range-based reconstruction for particle "{particle_type}" not available.')

def fiducial_cut(points, meta, margin=5):
    """
    A boolean selection function that determines if the given 3D space
    points are entirely contained within the volume with the given
    fiducial padding in each direction.
    
    Parameters
    ----------
    points: The space-points of the particle object.
    meta: The output from the sparse_meta3d parser.
    margin: The fiducial padding (in voxels) from the edge of the image.
    
    Returns
    -------
    True if all points are within the margin from the detector edges.
    """
    r = [meta[6+i]*points[:,i] + meta[i] for i in [0,1,2]]
    x_cut = np.all(((r[0] > -358.49 + margin) & (r[0] < -61.94 - margin))
                   |((r[0] < 358.49 - margin) & (r[0] > 61.94 + margin)))
    y_cut = np.all((r[1] > -181.86 + margin) & (r[1] < 134.96 - margin))
    z_cut = np.all((r[2] > -894.95 + margin) & (r[2] < 894.95 - margin))
    return (x_cut and y_cut and z_cut)
    #return not (np.any(points >= 768 - margin) or np.any(points <= margin))

def format_particle_string(particle_counts):
    """
    Returns a string denoting the particle counts in the interaction.
    
    Parameters
    ----------
    particle_counts: A dictionary keyed by PID number (0-4) describing
    particle the type and the associate count in the interaction.
    
    Returns
    -------
    A string of the format "VphWeXmuYpiZp" where V, W, X, Y, and Z
    denote the number of photons, electrons, muons, pions, and protons.
    """
    if isinstance(particle_counts, list):
        return (f'{particle_counts[0]}ph'
                + f'{particle_counts[1]}e'
                + f'{particle_counts[2]}mu'
                + f'{particle_counts[3]}pi'
                + f'{particle_counts[4]}p')
    else: 
        return (f'{particle_counts["Photon"]}ph'
                + f'{particle_counts["Electron"]}e'
                + f'{particle_counts["Muon"]}mu'
                + f'{particle_counts["Pion"]}pi'
                + f'{particle_counts["Proton"]}p')
        
def length(points, cfg):
    """
    Segments the 3D space point into multiple track segments by
    crawling along the track and clustering points within a sphere
    of the configured radius. The total length is taken to be the
    sum of the segment lengths.
    
    Parameters
    ----------
    points: A list of 3D space points for the track.
    cfg: the configuration dictionary.
    
    Returns
    -------
    The calculated length of the track.
    """
    return 0

def calo_reconstruct(points, depositions, work_func=23.6e-6,
                     gain=86.83, R=0.66, tau=3000, drift_v=1.6):
    """
    Applies electron lifetime corrections and a basic energy
    reconstruction.

    Parameters
    ----------
    points: (N,3) list of 3D space points.
    depositions: (N,) list of corresponding depositions.
    work_func: Work function of argon [MeV].
    gain: Gain of detector [e-/ADC].
    R: Average recombination.
    tau: Electron lifetime [us].
    drift_v: Electron drift velocity [mm/us].
    
    Return
    ------
    The estimated energy of the cluster.
    """
    boundaries = [-210.215, 0, 210.215]
    drift = np.array([-358.49, -61.94, 61.94, 358.49])
    indices = np.digitize(points[:,0], boundaries)
    drift_t = 10*np.abs(drift[indices] - points[:,0]) / drift_v
    return work_func * gain * (1/R) * np.sum(depositions * np.exp(drift_t / tau))
    
def cathode_crosser_cut(x0, x1):
    """
    A boolean selection function that determines if the track crosses
    the cathode based upon the coordinates of the most extreme x-points.
    
    Parameters
    ----------
    x0: X-coordinate of the minimum point along the x-direction.
    x1: X-coordinate of the maximum point along the x-direction.
    
    Return
    ------
    True if x0 and x1 are in different TPCs (different cryostats?).
    """
    dividers = [-210.215, 0, 210.215]
    x0_tpc = np.digitize(x0, dividers)
    x1_tpc = np.digitize(x1, dividers)
    return (x0_tpc != x1_tpc)

def boundary_crossing_cut(r0, r1, which='y'):
    """
    A boolean selection function that determines if the track crosses
    a detector boundary in the given direction.
    
    Parameters
    ----------
    r0: The coordinate of the minimum point along the given direction.
    r1: The coordinate of the maximum point along the given direction.
    
    Returns
    -------
    True if r0 and r1 are separated by a detector boundary.
    """
    dividers = {'x': [-358.49, 358.49],
                'y': [-181.86, 134.96],
                'z': [-894.95, 894.95]}[which]
    c0 = np.digitize(r0, dividers)
    c1 = np.digitize(r1, dividers)
    return (c0 != 1 or c1 != 1)

def segmentize(points, depositions, cfg):
    """
    Divides the track into segments starting from one end of the track
    and crawling along to the other end. The ends are determined using
    a pca fit_transform and the clustering of each segment is done via
    a spherical method.
    
    Parameters
    ----------
    points: The 3D space points of the track object.
    depositions: The charge depositions corresponding to the points.
    cfg: A configuration dictionary containing settings for the method.
    
    Return
    ------
    A tuple containing the assigned segment index for each space point
    as the first entry, and the segment lengths as the second entry.
    """
    pca = PCA(n_components=2)
    not_segmented = np.repeat(True, len(points))
    segment_id = np.repeat(-1, len(points))
    segment_length = list()
    min_voxels = cfg['analysis']['min_voxels']
    radius = cfg['analysis']['radius']
    s = 0
    
    while np.sum(not_segmented) > min_voxels:
        projection = pca.fit_transform(points[not_segmented])
        if np.all(not_segmented):
            p0 = points[not_segmented][np.argmin(projection[:, 0])]
            p1 = points[not_segmented][np.argmax(projection[:, 0])]
            q0 = np.sum(depositions[cdist([p0], points)[0] < radius])
            q1 = np.sum(depositions[cdist([p1], points)[0] < radius])
            start = (p0 if q0 > q1 else p1 )
        else:
            start = points[not_segmented][np.argmin(cdist([end], points[not_segmented])[0])]
        mask = cdist([start], points)[0] < radius
        end = points[mask & not_segmented][np.argmax(cdist([start], points[mask & not_segmented])[0])]
        if np.sum(mask & not_segmented) > min_voxels:
            segment_id[mask & not_segmented] = s
            s += 1
            local_projection = pca.fit_transform(points[mask & not_segmented])
            segment_length.append(local_projection[:, 0].max() - local_projection[:, 0].min())
        not_segmented[mask] = np.repeat(False, np.sum(mask))
    if np.sum(not_segmented) > 1:
        segment_length.append(np.amax(cdist(points[not_segmented], points[not_segmented])))
        segment_id[not_segmented] = s
    return np.array(segment_id), 0.3 * np.array(segment_length)

def track_dqdx(points, depositions, cfg):
    """
    Calculates the residual range, dQ, segment length, and dQ/dx for
    each track segment as defined by the segmentize function.
    
    Parameters
    ----------
    points: The 3D space points for the track.
    depositions: The charge depositions corresponding to the points.
    cfg: A configuration dictionary containing settings for the method.
    
    Returns
    -------
    The residual range, segment dQ, and segment length for each track
    segment.
    """
    segment_ids, segment_lengths = segmentize(points, depositions, cfg)
    residual_range = np.cumsum(segment_lengths) - (segment_lengths / 2.0)
    ids = np.unique(segment_ids)
    ids = ids[ids != -1]
    dqs = [np.sum(depositions[segment_ids == i]) for i in ids]
    return residual_range, dqs, segment_lengths, segment_ids

def exponential_fit(residual_range, dqdx):
    """
    Fits an exponential to the dQ/dx vs. residual range points.
    
    Parameters
    ----------
    residual_range: The residual_range for each segment.
    dqdx: The dQ/dx for each segment.
    
    Returns
    -------
    A tuple (a,b) containing the fit parameters for a*e^bx.
    """
    warnings.filterwarnings('ignore')
    if len(residual_range) > 2:
        try:
            return curve_fit(lambda t,a,b: a*np.exp(b*t), residual_range, dqdx)[0]
        except RuntimeError:
            return [-9999,-9999]
    else:
        return [-9999,-9999]
    
def get_children(interaction):
    """
    Counts the daughter particles in the interaction by type.
    
    Parameters
    ----------
    interaction: The interaction to search for children.
    
    Returns
    -------
    A list containing the count of children of each type.
    """
    particles = {}
    for p in interaction.particles: particles[p.particle_asis.id()] = p.pid
    children = dict()
    for p in interaction.particles:
        children_counts = [0,0,0,0,0]
        for c in p.particle_asis.children_id():
            if c in particles and c != p.id:
                children_counts[particles[c]] += 1
        children[p.particle_asis.id()] = children_counts
    return children

def get_track_theta(points, depositions):
    """
    Calculates the approximate angle of the track relative to the
    z-direction. Uses only the first part of the track to make the
    determination.
    
    Parameters
    ----------
    points: The 3D space points of the track
    depositions: The charge depositions corresponding to the points.
    
    Return
    ------
    The angle of the track with respect to the z-direction
    """
    pca = PCA(n_components=2)    
    radius = 20
    try:
        projection = pca.fit_transform(points)
        p0 = points[np.argmin(projection[:, 0])]
        p1 = points[np.argmax(projection[:, 0])]
        q0 = np.sum(depositions[cdist([p0], points)[0] < radius])
        q1 = np.sum(depositions[cdist([p1], points)[0] < radius])
        start = (p0 if q1 > q0 else p1 )
        mask = cdist([start], points)[0] < radius
        if np.sum(mask) > 2:
            primary = pca.fit(points[mask]).components_[0]
            return (primary[2] / np.linalg.norm(primary))
        else:
            return -9999
    except ValueError:
        return -9999
    
def get_endpoints(points, depositions, radius=20):
    """
    Calculates the start/end points of the track using local charge
    density to guess at the Bragg peak.
    
    Parameters
    ----------
    points: The 3D space points of the track
    depositions: the charge depositions corresponding to the points.
    radius: Radius (in voxels) for local charge density calculation.
    
    Return
    ------
    A list with two numpy arrays of shape (3,) containing the start
    and end point (respectively).
    """
    pca = PCA(n_components=2)
    projection = pca.fit_transform(points)
    candidates = [points[np.argmin(projection[:,0])],
                  points[np.argmax(projection[:,0])]]
    local_density = list()
    for ci, c in enumerate(candidates):
        mask = (cdist([c], points)[0] < radius)
        if(np.sum(mask) > 10):
            local_projection = pca.fit_transform(points[mask])
            local_candidates = [points[mask][np.argmin(local_projection[:,0])],
                                points[mask][np.argmax(local_projection[:,0])]]
            candidates[ci] = local_candidates[np.argmin(cdist([c], local_candidates)[0])]
            mask = (cdist([candidates[ci]], points)[0] < radius)
        local_density.append(np.sum(depositions[mask]))
    if np.argmin(local_density) == 1:
        local_density.reverse()
        candidates.reverse()
    return candidates
    
def get_track_angles(points, depositions, vertex, radius=20):
    """
    Calculates the approximate angle of the track using a local
    PCA about each endpoint.
    
    Parameters
    ----------
    points: The 3D space points of the track
    depositions: The charge depositions corresponding to the points.
    vertex: The vertex of the parent interaction.
    radius: Radius used for local primary direction calculation.
    
    Return
    ------
    The primary components of the track (starting at both ends),
    a bool flagging the first endpoint as having the lowest local
    charge density, and the calculated endpoints.
    """
    pca = PCA(n_components=2)
    try:
        endpoints = get_endpoints(points, depositions, radius)
        ret = list()
        centroid = list()
        p0_localQ_lowest = True
        for p in endpoints:
            mask = cdist([p], points)[0] < radius
            if np.sum(mask) > 2:
                centroid.append(np.mean(points[mask], axis=0))
                primary = pca.fit(points[mask]).components_[0]
                ret.append(primary / np.linalg.norm(primary))
            else:
                ret.append(np.array([-9999.0, -9999.0, -9999.0]))
                centroid.append(np.array([0,0,0]))
        if vertex[0] > 0:
            if np.argmin(cdist([vertex], endpoints)[0]) == 1:
                endpoints.reverse()
                ret.reverse()
                centroid.reverse()
                p0_localQ_lowest = False
            for ri in range(len(ret)):
                v = centroid[ri]
                cosine = (np.dot((v-vertex), ret[ri])
                          / (np.linalg.norm(v-vertex)
                             * np.linalg.norm(ret[ri])))
                if cosine < -0.5: ret[ri] = -1*ret[ri]
        return ret, p0_localQ_lowest, endpoints
    except ValueError:
        return None, None, None
