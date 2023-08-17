# Jinseok Ryu
# Electron Microscopy and Spectroscopy Lab.
# Seoul National University
# 20230710
# Simple implementation of PCA or NMF for feature extraction from a spectrum image including 4D-STEM data
# Some parameters for PCA and NMF were already pre-determined (you can change the parameters in this script)
# If the unit of the scan shape is '??m', it must be converted into another length unit, e.g., 'nm'
# It is because there is an error in reading greek alphabets from the calibration information.
# Please visit the documentation page of Scikit-learn package for detailed description of PCA and NMF
# https://scikit-learn.org/stable/index.html


# ********************************************************************************
print("Execute Python script in GMS 3")

import numpy as np
from scipy import optimize
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

def threed_roll_axis(img):
    stack = np.rollaxis(img, 0, 3)
    return stack
    
def fourd_roll_axis(stack):
    stack = np.rollaxis(np.rollaxis(stack, 2, 0), 3, 1)
    return stack

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
    
def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*np.exp(-(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments"""
    total = data.sum()
    X, Y = np.indices(data.shape) # row, col
    x = (X*data).sum()/total # row
    y = (Y*data).sum()/total # col
    col = data[:, int(y)]
    width_x = np.sqrt(np.abs((np.arange(col.size)-y)**2*col).sum()/col.sum()) # row
    row = data[int(x), :]
    width_y = np.sqrt(np.abs((np.arange(row.size)-x)**2*row).sum()/row.sum()) # col
    height = data.max()
    return height, x, y, width_x, width_y

def fitgaussian(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    errorfunction = lambda p: np.ravel(gaussian(*p)(*np.indices(data.shape)) -
                                 data)
    p, success = optimize.leastsq(errorfunction, params)
    return p


def gaussian_center(image, cbox_edge=0):
    y, x = np.indices(image.shape)
    if not cbox_edge:
        center = np.array([(y.max()-y.min())/2.0, (x.max()-x.min())/2.0])
        
    else:
        cbox_outy = int(image.shape[0]/2 - cbox_edge/2)
        cbox_outx = int(image.shape[1]/2 - cbox_edge/2)
        center_box = image[cbox_outy:-cbox_outy, cbox_outx:-cbox_outx]
        fit_params = fitgaussian(center_box)
        (_, center_y, center_x, _, _) = fit_params
        center = [center_y+cbox_outy, center_x+cbox_outx]
        
    return center
# ********************************************************************************

SI = DM.GetFrontImage()
print(SI)

SI_data = SI.GetNumArray()
SI_data = SI_data.clip(min=0.0)
num_dim = len(SI_data.data.shape)

if num_dim == 3:
    origin0, scale0, unit0 = SI.GetDimensionCalibration(0, 0)
    print(origin0, scale0, unit0)
    origin1, scale1, unit1 = SI.GetDimensionCalibration(1, 0)
    print(origin1, scale1, unit1)
    origin2, scale2, unit2 = SI.GetDimensionCalibration(2, 0)
    print(origin2, scale2, unit2)
    
    SI_data = threed_roll_axis(SI_data)
    print(SI_data.shape)
    
    e_range = np.arange(origin2, scale2*SI_data.shape[2]+origin2, scale2)
    print(e_range[0], e_range[-1])
    print(e_range.shape)
    
    crop_check = input("Do you want to crop spectra ? (Y or N)")

    if crop_check == "Y":
        start_eV = eval(input("initial energy loss (eV) of the crop range: "))
        end_eV = eval(input("final energy loss (eV) of the crop range: "))
        start_ind = find_nearest(e_range, start_eV)
        end_ind = find_nearest(e_range, end_eV)
        print(start_ind, end_ind)
        cr_range = [start_ind, end_ind]
        SI_data_cropped = SI_data[:, :, cr_range[0]:cr_range[1]].copy()

    elif crop_check == "N":
        SI_data_cropped = SI_data

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


elif num_dim == 4:
    origin0, scale0, unit0 = SI.GetDimensionCalibration(0, 0)
    print(origin0, scale0, unit0)
    origin1, scale1, unit1 = SI.GetDimensionCalibration(1, 0)
    print(origin1, scale1, unit1)
    origin2, scale2, unit2 = SI.GetDimensionCalibration(2, 0)
    print(origin2, scale2, unit2)
    origin3, scale3, unit3 = SI.GetDimensionCalibration(3, 0)
    print(origin3, scale3, unit3)

    SI_data = fourd_roll_axis(SI_data)
    
    pacbed = np.mean(SI_data, axis=(0,1))

    #find center position
    q_text_1 = """Select one option for finding the center position.
    1: Gaussian fitting - PACBED
    2: Center of mass - PACBED"""

    q_check_1 = int(input(q_text_1))

    if q_check_1 == 1:
        cbox_edge = int(input("size of the fitting box (data index): "))
        ct = gaussian_center(pacbed, cbox_edge=cbox_edge)
        print("center position")
        print(ct)

    elif q_check_1 == 2:
        cbox_edge = int(input("size of the fitting box (data index): "))
        cbox_outy = int(pacbed.shape[0]/2 - cbox_edge/2)
        cbox_outx = int(pacbed.shape[1]/2 - cbox_edge/2)
        center_box = pacbed[cbox_outy:-cbox_outy, cbox_outx:-cbox_outx]
        Y, X = np.indices(center_box.shape)
        com_y = np.sum(center_box * Y) / np.sum(center_box)
        com_x = np.sum(center_box * X) / np.sum(center_box)
        ct = [com_y+cbox_outy, com_x+cbox_outx]
        print("center position")
        print(ct)
        
    else:
        print("*"*50)
        print("wrong input !")
        print("*"*50)
        exit()
    
    fb = int(input("size of the box flattened (length of the center square box): "))

    dataset_input = []
    for i in range(SI_data.shape[0]):
        for j in range(SI_data.shape[1]):
            flat_temp = SI_data[i, j, int(ct[0]-fb/2):int(ct[0]+fb/2), int(ct[1]-fb/2):int(ct[1]+fb/2)].flatten()
            dataset_input.append(flat_temp)


    dataset_input = np.asarray(dataset_input).reshape(-1, fb*fb)
    print(dataset_input.shape)
    
    data_shape = (SI_data.shape[0], SI_data.shape[1], fb, fb)
    depth = fb*fb
    
    normalize_check = input("Do you want to normalize each spectrum ? (max normalization) (Y or N)")
    print(dataset_input.shape)

    if normalize_check == "Y":
        dataset_input = dataset_input / np.max(dataset_input, axis=1)[:, np.newaxis]


# ********************************************************************************

# ********************************************************************************

q_text = """Select one option.
1: PCA (principal component analysis)
2: NMF (non-negative matrix factorization)"""

decomp_check = int(input(q_text))

num_comp = int(input("How many loading vectors do you want to extract ?"))

if decomp_check == 1:
    pca_num_comp = num_comp
    skl_pca = PCA(n_components=pca_num_comp, whiten=False)
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
    
    if num_dim == 3:
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
            
    elif num_dim == 4:
        pca_comps_tmp = pca_comps.reshape(-1, 1, fb, fb)
        pca_comps_tmp = fourd_roll_axis(pca_comps_tmp)
        pca_comps_dm = DM.CreateImage(pca_comps_tmp.copy())
        pca_comps_dm.SetName("PCA loading vectors")
        pca_comps_dm.SetDimensionCalibration(1, 1, 1, "loading vector", 0)
        pca_comps_dm.SetDimensionCalibration(2, origin2, scale2, unit2, 0)
        pca_comps_dm.SetDimensionCalibration(3, origin3, scale3, unit3, 0)

        pca_coeffs_tmp = np.reshape(pca_coeffs, (data_shape[0], data_shape[1], pca_num_comp, 1))
        pca_coeffs_dm = DM.CreateImage(pca_coeffs_tmp.copy())
        pca_coeffs_dm.SetName("PCA coefficient maps")
        pca_coeffs_dm.SetDimensionCalibration(1, 1, 1, "loading vector", 0)
        pca_coeffs_dm.SetDimensionCalibration(2, origin0, scale0, unit0, 0)
        pca_coeffs_dm.SetDimensionCalibration(3, origin1, scale1, unit1, 0)
        
        pca_rec_tmp = np.reshape(pca_reconstructed, (data_shape[0], data_shape[1], fb, fb))
        pca_rec_tmp = fourd_roll_axis(pca_rec_tmp)
        pca_rec_dm = DM.CreateImage(pca_rec_tmp.copy())
        pca_rec_dm.SetName("PCA reconstructed SI")
        pca_rec_dm.SetDimensionCalibration(0, origin0, scale0, unit0, 0)
        pca_rec_dm.SetDimensionCalibration(1, origin1, scale1, unit1, 0)
        pca_rec_dm.SetDimensionCalibration(2, origin2, scale2, unit2, 0)
        pca_rec_dm.SetDimensionCalibration(3, origin3, scale3, unit3, 0)
    
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
    skl_nmf = NMF(n_components=nmf_num_comp, init="nndsvda", solver="mu", max_iter=1000, verbose=True, beta_loss="frobenius")
    nmf_coeffs = skl_nmf.fit_transform(dataset_input)
    print(nmf_coeffs[:, [1,2]].shape)
    nmf_comps = skl_nmf.components_
    
    # ********************************************************************************
    
    if num_dim == 3:
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
        nmf_coeffs_dm.SetDimensionCalibration(0, 1, 1, "loading vector", 0)
        nmf_coeffs_dm.SetDimensionCalibration(1, origin0, scale0, unit0, 0)
        nmf_coeffs_dm.SetDimensionCalibration(2, origin1, scale1, unit1, 0)

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
            
    elif num_dim == 4:
        nmf_comps_tmp = nmf_comps.reshape(-1, 1, fb, fb)
        nmf_comps_tmp = fourd_roll_axis(nmf_comps_tmp)
        nmf_comps_dm = DM.CreateImage(nmf_comps_tmp.copy())
        nmf_comps_dm.SetName("NMF loading vectors")
        nmf_comps_dm.SetDimensionCalibration(1, 1, 1, "loading vector", 0)
        nmf_comps_dm.SetDimensionCalibration(2, origin2, scale2, unit2, 0)
        nmf_comps_dm.SetDimensionCalibration(3, origin3, scale3, unit3, 0)

        nmf_coeffs_tmp = np.reshape(nmf_coeffs, (data_shape[0], data_shape[1], nmf_num_comp, 1))
        nmf_coeffs_dm = DM.CreateImage(nmf_coeffs_tmp.copy())
        nmf_coeffs_dm.SetName("NMF coefficient maps")
        nmf_coeffs_dm.SetDimensionCalibration(1, 1, 1, "loading vector", 0)
        nmf_coeffs_dm.SetDimensionCalibration(2, origin0, scale0, unit0, 0)
        nmf_coeffs_dm.SetDimensionCalibration(3, origin1, scale1, unit1, 0)

        # ********************************************************************************

        nmf_reconstructed = np.dot(nmf_coeffs, nmf_comps).reshape(data_shape[0], data_shape[1], fb, fb)
        print(nmf_coeffs.shape)
        print(nmf_comps.shape)
        print(nmf_reconstructed.shape)

        nmf_rec_tmp = fourd_roll_axis(nmf_reconstructed)
        nmf_rec_dm = DM.CreateImage(nmf_rec_tmp.copy())
        nmf_rec_dm.SetName("NMF reconstructed SI")
        nmf_rec_dm.SetDimensionCalibration(0, origin0, scale0, unit0, 0)
        nmf_rec_dm.SetDimensionCalibration(1, origin1, scale1, unit1, 0)
        nmf_rec_dm.SetDimensionCalibration(2, origin2, scale2, unit2, 0)
        nmf_rec_dm.SetDimensionCalibration(3, origin3, scale3, unit3, 0)
        
    nmf_comps_dm.ShowImage()
    nmf_coeffs_dm.ShowImage()
    nmf_rec_dm.ShowImage()
    # ********************************************************************************

else:
    print("Wrong input ! (only 1 or 2 possible)")
    exit()