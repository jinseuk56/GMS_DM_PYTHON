# Jinseok Ryu
# Electron Microscopy and Spectroscopy Lab.
# Seoul National University
# 20201209
# basic operation of NMF for EELS-SI


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
nmf_num_comp = 5
nmf_rec_num = nmf_num_comp
skl_nmf = NMF(n_components=nmf_num_comp, init="nndsvda", solver="mu", max_iter=1000, 
             verbose=True, beta_loss="frobenius", l1_ratio=0.1, alpha=1.0)

nmf_coeffs = skl_nmf.fit_transform(dataset_input)
nmf_comps = skl_nmf.components_
nmf_reconstructed = skl_nmf.inverse_transform(nmf_coeffs[:, :nmf_rec_num])
print(nmf_coeffs.shape)
print(nmf_comps.shape)
print(nmf_reconstructed.shape)
# ********************************************************************************

# ********************************************************************************
nmf_comps_tmp = np.rollaxis(nmf_comps.reshape(-1, 1, depth), 2, 0)
nmf_comps_dm = DM.CreateImage(nmf_comps_tmp.copy())
nmf_comps_dm.SetName("NMF loading vectors")

nmf_coeffs_tmp = np.rollaxis(np.reshape(nmf_coeffs, (data_shape[0], data_shape[1], nmf_num_comp)), 2, 0)
nmf_coeffs_dm = DM.CreateImage(nmf_coeffs_tmp.copy())
nmf_coeffs_dm.SetName("NMF coefficient maps")

nmf_rec_tmp = np.rollaxis(np.reshape(nmf_reconstructed, (data_shape[0], data_shape[1], -1)), 2, 0)
nmf_rec_dm = DM.CreateImage(nmf_rec_tmp.copy())
nmf_rec_dm.SetName("NMF reconstructed SI")
# ********************************************************************************

# ********************************************************************************
nmf_comps_dm.ShowImage()
nmf_coeffs_dm.ShowImage()
nmf_rec_dm.ShowImage()
# ********************************************************************************