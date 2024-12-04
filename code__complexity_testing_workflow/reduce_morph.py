import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def reducr_morphology(file_location):
    
    file_location
    out_file_location = file_location + '/output/'
    swc_file_location = file_location +'/components/morphologies/'
    chopp = 5
    import glob
    for file in glob.glob(swc_file_location + '/*.swc'):
        swc_file_location = file


    # Load data from HDF5 file
    with h5py.File(out_file_location + '/v_report_all.h5', 'r') as f:
        swc_ids_beg = np.array(f['/report/single_neuron/mapping/swc_ids_beg'])
        swc_ids_end = np.array(f['/report/single_neuron/mapping/swc_ids_end'])
    
    # Load morphology data from CSV file
    morphology = pd.read_csv(swc_file_location, delimiter=' ', header=None)
    morphology_xyz = morphology.iloc[:, 2:6].values
    morphology_idx = morphology.iloc[:, 0].values
    morphology_father = morphology.iloc[:, 6].values
    morphology_cat = morphology.iloc[:, 1].values



    # Find soma index
    soma_loca = np.where(morphology_cat == 1)[0]
    print(soma_loca)
    soma_loca = soma_loca[0]
    print(soma_loca)

    soma_idx = morphology_idx[soma_loca]

    # Replace -1 in morphology_father with soma_idx
    morphology_father[morphology_father == -1] = soma_idx

    # Find locations in morphology_idx that match morphology_father
    morphology_father = np.searchsorted(morphology_idx, morphology_father)

    morphology_father[soma_loca] = -1

    # Find locations in morphology_idx that match 
    swc_ids_beg = np.searchsorted(morphology_idx, swc_ids_beg)
    swc_ids_end = np.searchsorted(morphology_idx, swc_ids_end)


    # Extract h5_xyz
    h5_xyz = morphology_xyz[swc_ids_beg]
    
    # Plot h5_xyz
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(h5_xyz[:, 0], h5_xyz[:, 1], h5_xyz[:, 2], s=h5_xyz[:, 3] * 100.0, c='b', marker='.')
    # ax.set_aspect('equal', adjustable='box')
    plt.savefig(out_file_location +'/red_morph_1.png')

    # Calculate somadistance
    somadistance = morphology_xyz - morphology_xyz[soma_loca, :]
    somadistance = np.sqrt(np.sum(somadistance ** 2, axis=1))
    
    # Set somadistance to 0 where morphology_cat == 2
    somadistance[morphology_cat == 2] = 0
    
    # Find index of max distance
    max_dist_idx = np.argmax(somadistance)
    # Trace path
    path = np.full_like(np.arange(len(somadistance)), np.nan)
    targeting = max_dist_idx
    path_idx = 0

    while targeting != -1 :
        path[path_idx] = targeting
        targeting = int(morphology_father[targeting])
        path_idx += 1
    
    path = path[0:path_idx]


    # Calculate diff_dist


    morphology_father[soma_loca] = soma_loca


    diff_dist = np.zeros_like(morphology_xyz)
    diff_dist = morphology_xyz - morphology_xyz[morphology_father.astype(int), :]
    diff_dist = np.sqrt(np.sum(diff_dist ** 2, axis=1))
    
    # Get diff_dist for the path
    diff_dist_path = diff_dist[path.astype(int)]
    
    # Accumulate diff_dist_path
    for ii in range(len(diff_dist_path) - 1, 0, -1):
        diff_dist_path[ii - 1] += diff_dist_path[ii]

    # Determine chopped points
    chopped_point = []
    for ii in range(1, chopp):
        threshold = np.max(diff_dist_path) / chopp * ii
        print(threshold)
        asd = diff_dist_path < threshold
        if np.any(asd):
            
            diff_dist_path_t = np.array(diff_dist_path)
            conditions= np.array(asd)
            diff_dist_path_t[conditions] = 0

            idx = np.argmin(diff_dist_path_t)
            chopped_point.append(idx)
    chopped_point.append(0)
    
    # Convert chopped_point to path indices
    chopped_point = path[chopped_point]

    # Extract path_xyz
    path_xyz = morphology_xyz[path.astype(int)]

    # Determine indices
    ind_other = ~np.isin(morphology_idx, path)
    ind_axon = (morphology_cat == 2)
    ind_other = ind_other & (~ind_axon)

    # Plot results
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # ax.scatter(morphology_xyz[ind_axon, 0], morphology_xyz[ind_axon, 1], morphology_xyz[ind_axon, 2],
            #    c='lightgray', s=1, marker='.')

    ax.scatter(morphology_xyz[ind_other, 0], morphology_xyz[ind_other, 1], morphology_xyz[ind_other, 2],
               c='orangered', s=5, marker='.')

    ax.scatter(path_xyz[:, 0], path_xyz[:, 1], path_xyz[:, 2],
               c='red', s=10, marker='.')

    ax.scatter(morphology_xyz[chopped_point, 0], morphology_xyz[chopped_point, 1], morphology_xyz[chopped_point, 2],
               c='blue', s=100, marker='.')

    # ax.set_aspect('equal', adjustable='box')

    plt.savefig(out_file_location + '/red_morph.png')

    # Compute morphology_h5_idx
    morphology_h5_idx = np.full(morphology_xyz.shape[0], np.nan)
    for ii in range(len(swc_ids_beg)):
        morphology_h5_idx[swc_ids_beg[ii]:swc_ids_end[ii] + 1] = ii
    
    # Compute line_idx
    line_idx = morphology_h5_idx[chopped_point.astype(int)]

    return np.append([soma_loca], line_idx.astype(int), axis=0)
