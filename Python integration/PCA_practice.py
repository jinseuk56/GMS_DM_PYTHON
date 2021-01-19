# Jinseok Ryu
# Electron Microscopy and Spectroscopy Lab.
# Seoul National University
# 20201209
# basic operation of PCA for EELS-SI


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


# ********************************************************************************
def zero_one_rescale(spectrum):
    """
    get rid of negative values
    rescale a spectrum [0, 1]
    """
    spectrum = spectrum.clip(min=0.0)
    min_val = np.min(spectrum)
    
    rescaled = spectrum - min_val
    
    if np.max(rescaled) != 0:
        rescaled = rescaled / np.max(rescaled)
    
    return rescaled
# ********************************************************************************


# ********************************************************************************
SI = DM.GetFrontImage()
print(SI)

SI_data = np.rollaxis(SI.GetNumArray(), 0, 3)
print(SI_data.shape)
# ********************************************************************************


# ********************************************************************************
cr_range = [250, 550]
SI_data_cropped = SI_data[:, :, cr_range[0]:cr_range[1]]
print(SI_data_cropped.shape)

data_shape = SI_data_cropped.shape[:2]
depth = SI_data_cropped.shape[2]

dataset_input = SI_data_cropped.reshape(-1, depth)
for i in range(len(dataset_input)):
	dataset_input[i] = zero_one_rescale(dataset_input[i])
print(dataset_input.shape)
# ********************************************************************************

# ********************************************************************************
pca_num_comp = 5
pca_rec = np.arange(pca_num_comp)
#pca_rec = np.array([0, 1])
skl_pca = PCA(n_components=pca_num_comp, whiten=False, svd_solver="auto")
pca_coeffs = skl_pca.fit_transform(dataset_input)
pca_comps = skl_pca.components_
pca_reconstructed = np.dot(pca_coeffs[:, pca_rec], pca_comps[pca_rec]) + skl_pca.mean_

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

pca_coeffs_tmp = np.rollaxis(np.reshape(pca_coeffs, (data_shape[0], data_shape[1], pca_num_comp)), 2, 0)
pca_coeffs_dm = DM.CreateImage(pca_coeffs_tmp.copy())
pca_coeffs_dm.SetName("PCA coefficient maps")


pca_rec_tmp = np.rollaxis(np.reshape(pca_reconstructed, (data_shape[0], data_shape[1], -1)), 2, 0)
pca_rec_dm = DM.CreateImage(pca_rec_tmp.copy())
pca_rec_dm.SetName("PCA reconstructed SI")
# ********************************************************************************

# ********************************************************************************
pca_explained.ShowImage()
pca_comps_dm.ShowImage()
pca_coeffs_dm.ShowImage()
pca_rec_dm.ShowImage()
# ********************************************************************************