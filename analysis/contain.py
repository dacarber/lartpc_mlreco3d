'''
Containment function
Finds if voxel of a particle is contained within the TPC

Input:
    particle: particle object
    data: data from from lartpc_mlreco_3d
'''

def containment_mask_noghost(particle, data):
    '''
    Containment function
    Finds if voxel of a particle is contained within the TPC

    Input:
        particle: particle object
        data: data from from lartpc_mlreco_3d
    Return:
        contained_mask: Array of true and false that corresponds to voxel in the particle
    '''
    min_x, min_y, min_z = data['meta'][entry][0:3]
    max_x, max_y, max_z = data['meta'][entry][3:6]
    size_voxel_x, size_voxel_y, size_voxel_z = data['meta'][entry][6:9]
    contained_mask = []
    for v in particle.coords_noghost:
        absolute_pos = v * size_voxel_x + min_x
        if absolute_pos[0] <= -354.63 or (absolute_pos[0] >= -55.8 and absolute_pos[0] <= 55.8) or  absolute_pos[0] >= 354.63:
            contained_mask.append(False)
            continue
        elif absolute_pos[1] <= -176.86 or absolute_pos[1] >= 129.96:
            contained_mask.append(False)
            continue
        elif absolute_pos[2] <= -890.95 or absolute_pos[2] >= 889.95:
            contained_mask.append(False)
            continue
        else:
            contained_mask.append(True)
    return contained_mask

def contained_noghost(particle, data):
    '''
    Containment function
    Finds if a particle is contained within the TPC

    Input:
        particle: particle object
        data: data from from lartpc_mlreco_3d
    Return:
        out_of_detector: Boolean that says if the particle is contained or uncontained
    '''
    min_x, min_y, min_z = data['meta'][entry][0:3]
    max_x, max_y, max_z = data['meta'][entry][3:6]
    size_voxel_x, size_voxel_y, size_voxel_z = data['meta'][entry][6:9]
    out_of_detector = False
    for pt in particle.coords_noghost:
        absolute_pos = pt * size_voxel_x + min_x
        if absolute_pos[0] <= -359.63 or (absolute_pos[0] >= -60.8 and absolute_pos[0] <= 60.8) or  absolute_pos[0] >= 359.63:
            out_of_detector = True
            continue
        elif absolute_pos[1] <= -181.86 or absolute_pos[1] >= 134.96:
            out_of_detector = True
            continue
        elif absolute_pos[2] <= -895.95 or absolute_pos[2] >= 894.95:
            out_of_detector = True
            continue
    return out_of_detector
def containment_mask(particle, data):
    '''
    Containment function
    Finds if voxel of a particle is contained within the TPC

    Input:
        particle: particle object
        data: data from from lartpc_mlreco_3d
    Return:
        contained_mask: Array of true and false that corresponds to voxel in the particle
    '''
    min_x, min_y, min_z = data['meta'][entry][0:3]
    max_x, max_y, max_z = data['meta'][entry][3:6]
    size_voxel_x, size_voxel_y, size_voxel_z = data['meta'][entry][6:9]
    contained_mask = []
    for v in particle.points:
        absolute_pos = v * size_voxel_x + min_x
        if absolute_pos[0] <= -359.63 or (absolute_pos[0] >= -60.8 and absolute_pos[0] <= 60.8) or  absolute_pos[0] >= 359.63:
            contained_mask.append(False)
            continue
        elif absolute_pos[1] <= -181.86 or absolute_pos[1] >= 134.96:
            contained_mask.append(False)
            continue
        elif absolute_pos[2] <= -895.95 or absolute_pos[2] >= 894.95:
            contained_mask.append(False)
            continue
        else:
            contained_mask.append(True)
    return contained_mask

def contained(particle, data,margin = 5):
    '''
    Containment function
    Finds if a particle is contained within the TPC

    Input:
        particle: particle object
        data: data from from lartpc_mlreco_3d
    Return:
        out_of_detector: Boolean that says if the particle is contained or uncontained
    '''
    min_x, min_y, min_z = data['meta'][0][0:3]
    max_x, max_y, max_z = data['meta'][0][3:6]
    size_voxel_x, size_voxel_y, size_voxel_z = data['meta'][0][6:9]
    out_of_detector = False
    absolute_pos = [0,0,0]
    
    for p,pt in enumerate(particle.points):
        absolute_pos[0] = pt[0] #* size_voxel_x + min_x
        absolute_pos[1] = pt[1] #* size_voxel_y + min_y
        absolute_pos[2] = pt[2] #* size_voxel_z + min_z
        if particle.sources[p][0] == 1 and (particle.sources[p][1] == 2 or particle.sources[p][1] == 3): #WW
            if absolute_pos[0] >= (359.63-margin) or absolute_pos[0] <= (210.29+margin):
                out_of_detector = True
                #print(1,absolute_pos[0])
                break
            elif absolute_pos[1] <= (-181.86+margin) or absolute_pos[1] >= (134.96-margin):
                out_of_detector = True
                #print(2,absolute_pos[0])
                break
            elif absolute_pos[2] <= (-894.95+margin) or absolute_pos[2] >= (894.95-margin):
                out_of_detector = True
                #print(3,absolute_pos[0])
                break
        elif particle.sources[p][0] == 1 and (particle.sources[p][1] == 0 or particle.sources[p][1] == 1): #WE
            if absolute_pos[0] >= (210.29-margin) or absolute_pos[0] <= (60.8+margin):
                out_of_detector = True
                #print(4,absolute_pos[0])
                break
            elif absolute_pos[1] <= (-181.86+margin) or absolute_pos[1] >= (134.96-margin):
                out_of_detector = True
                #print(5,absolute_pos[0])
                break
            elif absolute_pos[2] <= (-894.95+margin) or absolute_pos[2] >= (894.95-margin):
                out_of_detector = True
                #print(6,absolute_pos[0])
                break
        elif particle.sources[p][0] == 0 and (particle.sources[p][1] == 2 or particle.sources[p][1] == 3): #EW
            if absolute_pos[0] <= (-210.29+margin) or absolute_pos[0] >= (-60.8-margin):
                out_of_detector = True
                #print(7,absolute_pos[0])
                break
            elif absolute_pos[1] <= (-181.86+margin) or absolute_pos[1] >= (134.96-margin):
                out_of_detector = True
                #print(8,absolute_pos[0])
                break
            elif absolute_pos[2] <= (-894.95+margin) or absolute_pos[2] >= (894.95-margin):
                out_of_detector = True
                #print(9,absolute_pos[0])
                break
        elif particle.sources[p][0] == 0 and (particle.sources[p][1] == 0 or particle.sources[p][1] == 1): #EE
            if absolute_pos[0] <= (-359.63+margin) or absolute_pos[0] >= (-210.29-margin):
                out_of_detector = True
                #print(10,absolute_pos[0])
                break
            elif absolute_pos[1] <= (-181.86+margin) or absolute_pos[1] >= (134.96-margin):
                out_of_detector = True
                #print(11,absolute_pos[0])
                break
            elif absolute_pos[2] <= (-894.95+margin) or absolute_pos[2] >= (894.95-margin):
                out_of_detector = True
                #print(12,absolute_pos[0])
                break
    return out_of_detector