# Jinseok Ryu
# Electron Microscopy and Spectroscopy Lab.
# Seoul National University
# 20220705
# simple implementation of PCA or NMF for feature extraction from an EELS spectrum image
# some parameters in PCA and NMF were already pre-determined (you can change the parameters in this script)
# please visit the documentation page of Scikit-learn package for detailed description of PCA and NMF
# https://scikit-learn.org/stable/index.html


# ********************************************************************************
print("Execute Python script in GMS 3")

import numpy as np
import DigitalMicrograph as DM
from sklearn.decomposition import NMF, PCA

#import sys
#sys.argv.extend(['-a', ' '])
#import matplotlib.pyplot as plt

print("Libraries have been imported completely")
# ********************************************************************************

if (False == DM.IsScriptOnMainThread()):
    print('MatplotLib scripts require to be run on the main thread.')
    exit()

# ********************************************************************************
SI = DM.GetFrontImage()
print(SI)

origin0, scale0, unit0 = SI.GetDimensionCalibration(0, 0)
print(origin0, scale0, unit0)
origin1, scale1, unit1 = SI.GetDimensionCalibration(1, 0)
print(origin1, scale1, unit1)
origin2, scale2, unit2 = SI.GetDimensionCalibration(2, 0)
print(origin2, scale2, unit2)


SI_data = np.rollaxis(SI.GetNumArray(), 0, 3)
print(SI_data.shape)
# ********************************************************************************

# ********************************************************************************

crop_check = input("Do you want to crop spectra ? (Y or N)")

if crop_check == "Y":
    start_ind = int(input("initial index of the crop range: "))
    end_ind = int(input("final index of the crop range: "))
    cr_range = [start_ind, end_ind]
    SI_data_cropped = SI_data[:, :, cr_range[0]:cr_range[1]].copy()

elif crop_check == "N":
    SI_data_cropped = SI_data.copy()

else:
    print("Wrong input ! (only Y or N possible)")
    exit()

data_shape = SI_data_cropped.shape[:2]
depth = SI_data_cropped.shape[2]

dataset_input = SI_data_cropped.reshape(-1, depth).clip(min=0.0)
normalize_check = input("Do you want to normalize each spectrum ? (max normalization) (Y or N)")
print(dataset_input.shape)

if normalize_check == "Y":
    dataset_input = dataset_input / np.max(dataset_input, axis=1)[:, np.newaxis]
# ********************************************************************************

q_text = """Select one option.
1: PCA (principal component analysis)
2: NMF (non-negative factorization method)"""

decomp_check = int(input(q_text))

num_comp = int(input("How many loading vectors do you want to extract ?"))

if decomp_check == 1:
    # ********************************************************************************
    pca_num_comp = num_comp
    skl_pca = PCA(n_components=pca_num_comp, whiten=False, svd_solver="auto")
    pca_coeffs = skl_pca.fit_transform(dataset_input)
    pca_comps = skl_pca.components_


    pca_reconstructed = np.dot(pca_coeffs[:, :pca_num_comp], pca_comps[:pca_num_comp]) + skl_pca.mean_

    print(pca_coeffs.shape)
    print(pca_comps.shape)
    print(pca_reconstructed.shape)
    # ********************************************************************************

    # ********************************************************************************
    pca_explained = DM.CreateImage(skl_pca.explained_variance_ratio_.copy())
    pca_explained.SetName("Explained variance ratio")

    pca_comps_tmp = np.rollaxis(pca_comps.reshape(-1, 1, depth), 2, 0)
    pca_comps_dm = DM.CreateImage(pca_comps_tmp.copy())
    pca_comps_dm.SetName("PCA loading vectors")
    pca_comps_dm.SetDimensionCalibration(1, 1, 1, "loading vector", 0)
    if crop_check == "Y":
        pca_comps_dm.SetDimensionCalibration(2, origin2+start_ind*scale2, scale2, unit2, 0)
    else:
        pca_comps_dm.SetDimensionCalibration(2, origin2, scale2, unit2, 0)

    pca_coeffs_tmp = np.reshape(pca_coeffs, (data_shape[0], data_shape[1], pca_num_comp, 1))
    pca_coeffs_dm = DM.CreateImage(pca_coeffs_tmp.copy())
    pca_coeffs_dm.SetName("PCA coefficient maps")
    pca_coeffs_dm.SetDimensionCalibration(1, 1, 1, "loading vector", 0)
    pca_coeffs_dm.SetDimensionCalibration(2, origin0, scale0, unit0, 0)
    pca_coeffs_dm.SetDimensionCalibration(3, origin1, scale1, unit1, 0)

    pca_rec_tmp = np.rollaxis(np.reshape(pca_reconstructed, (data_shape[0], data_shape[1], -1)), 2, 0)
    pca_rec_dm = DM.CreateImage(pca_rec_tmp.copy())
    pca_rec_dm.SetName("PCA reconstructed SI")
    pca_rec_dm.SetDimensionCalibration(0, origin0, scale0, unit0, 0)
    pca_rec_dm.SetDimensionCalibration(1, origin1, scale1, unit1, 0)
    
    if crop_check == "Y":
        pca_rec_dm.SetDimensionCalibration(2, origin2+start_ind*scale2, scale2, unit2, 0)
    else:
        pca_rec_dm.SetDimensionCalibration(2, origin2, scale2, unit2, 0)
    
    # ********************************************************************************

    # ********************************************************************************
    pca_explained.ShowImage()
    pca_comps_dm.ShowImage()
    pca_coeffs_dm.ShowImage()
    pca_rec_dm.ShowImage()
    # ********************************************************************************

elif decomp_check == 2:
    # ********************************************************************************
    nmf_num_comp = num_comp
    skl_nmf = NMF(n_components=nmf_num_comp, init="nndsvda", solver="mu", max_iter=1000, verbose=True, beta_loss="frobenius", l1_ratio=0.0, alpha=0.0)

    nmf_coeffs = skl_nmf.fit_transform(dataset_input)
    print(nmf_coeffs[:, [1,2]].shape)
    nmf_comps = skl_nmf.components_
    
    # ********************************************************************************
    nmf_comps_tmp = np.rollaxis(nmf_comps.reshape(-1, 1, depth), 2, 0)
    nmf_comps_dm = DM.CreateImage(nmf_comps_tmp.copy())
    nmf_comps_dm.SetName("NMF loading vectors")
    nmf_comps_dm.SetDimensionCalibration(1, 1, 1, "loading vector", 0)
    if crop_check == "Y":
        nmf_comps_dm.SetDimensionCalibration(2, origin2+start_ind*scale2, scale2, unit2, 0)
    else:
        nmf_comps_dm.SetDimensionCalibration(2, origin2, scale2, unit2, 0)

    nmf_coeffs_tmp = np.reshape(nmf_coeffs, (data_shape[0], data_shape[1], nmf_num_comp, 1))
    nmf_coeffs_dm = DM.CreateImage(nmf_coeffs_tmp.copy())
    nmf_coeffs_dm.SetName("NMF coefficient maps")
    nmf_coeffs_dm.SetDimensionCalibration(1, 1, 1, "loading vector", 0)
    nmf_coeffs_dm.SetDimensionCalibration(2, origin0, scale0, unit0, 0)
    nmf_coeffs_dm.SetDimensionCalibration(3, origin1, scale1, unit1, 0)
    
    nmf_comps_dm.ShowImage()
    nmf_coeffs_dm.ShowImage()
    
    # ********************************************************************************

    nmf_reconstructed = np.dot(nmf_coeffs, nmf_comps)
    print(nmf_coeffs.shape)
    print(nmf_comps.shape)
    print(nmf_reconstructed.shape)

    nmf_rec_tmp = np.rollaxis(np.reshape(nmf_reconstructed, (data_shape[0], data_shape[1], -1)), 2, 0)
    nmf_rec_dm = DM.CreateImage(nmf_rec_tmp.copy())
    nmf_rec_dm.SetName("NMF reconstructed SI")
    nmf_rec_dm.SetDimensionCalibration(0, origin0, scale0, unit0, 0)
    nmf_rec_dm.SetDimensionCalibration(1, origin1, scale1, unit1, 0)
    if crop_check == "Y":
        nmf_rec_dm.SetDimensionCalibration(2, origin2+start_ind*scale2, scale2, unit2, 0)
    else:
        nmf_rec_dm.SetDimensionCalibration(2, origin2, scale2, unit2, 0)
    
    nmf_rec_dm.ShowImage()
    # ********************************************************************************

else:
    print("Wrong input ! (only 1 or 2 possible)")
    exit()