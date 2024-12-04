import numpy as np
import pandas as pd


def dendrite_surface_area(swc_file_location):
    # swc_file_location = file_location +'/components/morphologies/'
    # import glob
    # for file in glob.glob(swc_file_location + '/*.swc'):
    #     swc_file_location = file


    morphology = pd.read_csv(swc_file_location, delimiter=' ', header=None)
    # morphology_xyz = morphology.iloc[:, 2:6].values
    # morphology_idx = morphology.iloc[:, 0].values -1
    # morphology_father = morphology.iloc[:, 6].values -1
    # morphology_cat = morphology.iloc[:, 1].values
    # ro,co=morphology_xyz.shape
    # diff_dist = np.zeros_like(morphology_xyz[:,0])
    # for ii in range(1, ro):
    #     diff_dist[ii] = np.sqrt(np.sum((morphology_xyz[ii,:] - morphology_xyz[morphology_father[ii],:]) ** 2))
    # diff_dist[morphology_cat == 2] = 0
    # return np.sum(  diff_dist * 3.14159 *  morphology_xyz[:,3] * 2.0 )

    # Extract columns from the morphology array
    morphology_xyz = morphology.iloc[:, 2:5].values
    morphology_r = morphology.iloc[:, 5].values

    morphology_idx = morphology.iloc[:, 0].values
    morphology_father = morphology.iloc[:, 6].values
    morphology_cat = morphology.iloc[:, 1].values



    # Find soma index
    soma_idx = np.where(morphology_cat == 1)[0]
    soma_idx = morphology_idx[soma_idx[0]]

    # Replace -1 in morphology_father with soma_idx
    morphology_father[morphology_father == -1] = soma_idx

    # Find locations in morphology_idx that match morphology_father
    morphology_father_loc = np.searchsorted(morphology_idx, morphology_father)

    # Compute distance to father
    father_xyz = morphology_xyz[morphology_father_loc]
    dist_to_father = np.sqrt(np.sum((morphology_xyz - father_xyz) ** 2, axis=1))

    # Determine which rows to include in the calculation
    flagg = morphology_cat != 2

    # Compute surface area
    surf_area = np.sum(dist_to_father[flagg] * 2.0 * np.pi * morphology_r[flagg])

    return surf_area