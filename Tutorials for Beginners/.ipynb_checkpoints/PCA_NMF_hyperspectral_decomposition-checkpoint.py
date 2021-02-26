# Jinseok Ryu
# Electron Microscopy and Spectroscopy Lab.
# Seoul National University
# 20210226
# basic applications of PCA or NMF to hyperspectral data (e.g. EELS-SI)


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

if ( False == DM.IsScriptOnMainThread() ):
	print('MatplotLib scripts require to be run on the main thread.')
	exit()

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

dataset_input = SI_data_cropped.reshape(-1, depth)
for i in range(len(dataset_input)):
	dataset_input[i] = zero_one_rescale(dataset_input[i])
print(dataset_input.shape)
# ********************************************************************************

q_text = """Select one option.
1: PCA (principal component analysis)
2: NMF (non-negative factorization method)"""

decomp_check = int(input(q_text))

num_comp = int(input("How many loading vectors do you want to extract ?"))
num_rec = int(input("How many loading vectors do you want to use when reconstructing the data ?"))

if decomp_check == 1:
	# ********************************************************************************
	pca_num_comp = num_comp
	pca_rec = np.arange(num_rec)
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

	pca_coeffs_tmp = np.reshape(pca_coeffs, (data_shape[0], data_shape[1], pca_num_comp, 1))
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

elif decomp_check == 2:
	# ********************************************************************************
	nmf_num_comp = num_comp
	nmf_rec = np.arange(num_rec)
	#nmf_rec = np.array([0,1])
	skl_nmf = NMF(n_components=nmf_num_comp, init="nndsvda", solver="mu", max_iter=1000, 
				 verbose=True, beta_loss="frobenius", l1_ratio=0.0, alpha=0.0)

	nmf_coeffs = skl_nmf.fit_transform(dataset_input)
	print(nmf_coeffs[:, [1,2]].shape)
	nmf_comps = skl_nmf.components_
	nmf_reconstructed = np.dot(nmf_coeffs[:, nmf_rec], nmf_comps[nmf_rec])
	print(nmf_coeffs.shape)
	print(nmf_comps.shape)
	print(nmf_reconstructed.shape)
	# ********************************************************************************

	# ********************************************************************************
	nmf_comps_tmp = np.rollaxis(nmf_comps.reshape(-1, 1, depth), 2, 0)
	nmf_comps_dm = DM.CreateImage(nmf_comps_tmp.copy())
	nmf_comps_dm.SetName("NMF loading vectors")

	nmf_coeffs_tmp = np.reshape(nmf_coeffs, (data_shape[0], data_shape[1], nmf_num_comp, 1))
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

else:
	print("Wrong input ! (only 1 or 2 possible)")
	exit()