import numpy as np
import ROOT
from sklearn.decomposition import PCA
from scipy.optimize import curve_fit
from scipy.spatial.distance import cdist
import warnings
import warnings, yaml
warnings.filterwarnings('ignore')
import analysis.reco as reco
import analysis.contain as contain

def shower_energy(particle,data,entry):
    tdrift = 0
    min_x, min_y, min_z = data['meta'][entry][0:3]
    max_x, max_y, max_z = data['meta'][entry][3:6]
    size_voxel_x, size_voxel_y, size_voxel_z = data['meta'][entry][6:9]
    out_of_detector = False
    for pt in range(len(particle.points)):
        absolute_pos = particle.points[pt] * size_voxel_x + min_x
        if absolute_pos[0] >= -359.63 and absolute_pos[0] <= -210.29:
            anode_dist = abs(absolute_pos[0] +359.63)
            tdrift =  10*anode_dist/ (1.6)
        if absolute_pos[0] >= -210.14 and absolute_pos[0] <= -60.8:
            anode_dist = abs(absolute_pos[0] +60.8)
            tdrift =  10*anode_dist/ (1.6)
        if absolute_pos[0] <= 210.14 and absolute_pos[0] >= 60.8:
            anode_dist = abs(absolute_pos[0] -60.8)
            tdrift =  10*anode_dist/ (1.6)
        if absolute_pos[0] <= 359.63 and absolute_pos[0] >= 210.14:
            anode_dist = abs(absolute_pos[0] -359.63)
            tdrift =  10*anode_dist/ (1.6)
        if absolute_pos[0] >= 359.63 and absolute_pos[0] <= -359.63:
            out_of_detector = True
    total_energy = np.sum(particle.depositions)*(23.6*10**(-6))*(85.25)*np.exp(tdrift/3000)*0.66**(-1)
    return total_energy
def true_energy(particle):
    total_energy = particle.asis.energy_deposit()
    return total_energy
def track_energy(particle):
    ana_cfg = yaml.load(open('/sdf/home/d/dcarber/DATADIR/ana_cfg.cfg', 'r').read(),Loader=yaml.Loader)
    bethe_bloch_energy = {2: reco.load_range_reco(ana_cfg['rrinput'],'muon',kinetic_energy = True), 3: reco.load_range_reco(ana_cfg['rrinput'],'pion',kinetic_energy = True), 4: reco.load_range_reco(ana_cfg['rrinput'],'proton',kinetic_energy = True)}
    rr, dqs, dxs, ids = reco.track_dqdx(particle.points, particle.depositions, ana_cfg)
    length = np.sum(dxs)
    recop = bethe_bloch_energy[particle.pid](length)
    return recop
def startpoint_dist(particle,interact):
    sp = particle.start_point
    inter_point = interact.vertex
    distance = cdist(sp.reshape(1,-1), inter_point.reshape(1,-1))
    return distance[0][0]
def directions(particle):
    direction = particle.start_dir
    #direction = cluster_direction(particle.points,particle.start_point,optimize = True)
    if direction[0] > 0:
        cluster_theta = (np.arctan((direction[1])/(direction[0])))
    elif direction[1] >= 0:
        cluster_theta = (np.arctan((direction[1])/(direction[0])) +np.pi)
    elif direction[1] < 0:
        cluster_theta = (np.arctan((direction[1])/(direction[0])) -np.pi)

    if direction[2] > 0:
        cluster_phi = (np.arctan((np.sqrt((direction[0])**(2)+(direction[1])**(2)))/(direction[2])))
    else:
        cluster_phi = (np.arctan((np.sqrt((direction[0])**(2)+(direction[1])**(2)))/(direction[2]))+np.pi)
    return cluster_theta, cluster_phi
def dqdx(particle, data, entry):

    pca = PCA(n_components=2)
    min_x, min_y, min_z = data['meta'][entry][0:3]
    max_x, max_y, max_z = data['meta'][entry][3:6]
    size_voxel_x, size_voxel_y, size_voxel_z = data['meta'][entry][6:9]
    tdrift = 0
    point_dist =[]
    min_segment_size = 3
    ppn_prediction = particle.start_point
    dist = cdist(particle.points, ppn_prediction.reshape(1,-1))#Distance of each point in particle from the start point'''
    mask = dist.squeeze() < 10
    selected_points = particle.points[mask]
    if selected_points.shape[0] < 2:# Makes sure there are enough points
        dqdx = None
        return dqdx
    absolute_pos = selected_points * size_voxel_x + min_x
    out_of_detector = False
    for pt in range(len(absolute_pos)):
        if absolute_pos[pt][0] >= -359.63 and absolute_pos[pt][0] <= -210.29:
            anode_dist = abs(absolute_pos[pt][0] +359.63)
            tdrift =  10*anode_dist/ (1.6)
        if absolute_pos[pt][0] >= -210.14 and absolute_pos[pt][0] <= -60.8:
            anode_dist = abs(absolute_pos[pt][0] +60.8)
            tdrift =  10*anode_dist/ (1.6)
        if absolute_pos[pt][0] <= 210.14 and absolute_pos[pt][0] >= 60.8:
            anode_dist = abs(absolute_pos[pt][0] -60.8)
            tdrift =  10*anode_dist/ (1.6)
        if absolute_pos[pt][0] <= 359.63 and absolute_pos[pt][0] >= 210.14:
            anode_dist = abs(absolute_pos[pt][0] -359.63)
            tdrift =  10*anode_dist/ (1.6)
        if absolute_pos[pt][0] >= 359.63 and absolute_pos[pt][0] <= -359.63:
            out_of_detector = True
    proj = pca.fit_transform(selected_points)
    dx = (proj[:,0].max() - proj[:,0].min())*.3
    if dx < min_segment_size:
        dqdx = None
        return dqdx
    dq = np.sum(particle.depositions[mask])*(23.6*10**(-6))*(85.25)*np.exp(tdrift/3000)*0.66**(-1)
    dqdx = dq/dx
    return dqdx
def true_directions(particle):
    true_mom_x = particle.asis.px()
    true_mom_y = particle.asis.py()
    true_mom_z = particle.asis.pz()
    true_vector = [true_mom_x,true_mom_y,true_mom_z]

    if true_mom_x > 0:
        true_theta = (np.arctan((true_mom_y)/(true_mom_x)))
    elif true_mom_y >= 0: 
        true_theta = (np.arctan((true_mom_y)/(true_mom_x)) + np.pi)
    elif true_mom_y < 0:
        true_theta = (np.arctan((true_mom_y)/(true_mom_x))- np.pi)

    if true_mom_z >0:
        true_phi = (np.arctan((np.sqrt((true_mom_x)**(2)+(true_mom_y)**(2)))/(true_mom_z)))
    else:
        true_phi = (np.arctan((np.sqrt((true_mom_x)**(2)+(true_mom_y)**(2)))/(true_mom_z))+np.pi)

    return true_theta, true_phi
