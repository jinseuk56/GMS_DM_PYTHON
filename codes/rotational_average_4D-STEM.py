# Jinseok Ryu
# Electron Microscopy and Spectroscopy Lab.
# Seoul National University
# last update : 20210226
# rotational average of each DP in 4D-STEM data


# ********************************************************************************
print("Execute Python script in GMS 3")

import sys
sys.argv.extend(['-a', ' '])
import DigitalMicrograph as DM
import numpy as np
from scipy import optimize

print("Libraries have been imported completely")
# ********************************************************************************

if ( False == DM.IsScriptOnMainThread() ):
	print('MatplotLib scripts require to be run on the main thread.')
	exit()
	
# refer to https://scipy-cookbook.readthedocs.io/items/FittingData.html

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


# refer to "github.com/mkolopanis/python/blob/master/radialProfile.py"

def radial_average_with_center(image, center=None):
    y, x = np.indices(image.shape)
    if not center:
        center = np.array([(y.max()-y.min())/2.0, (x.max()-x.min())/2.0])
        
    r = np.hypot(y - center[0], x - center[1])
    #plt.imshow(r, cmap="Accent")
    #plt.show()

    # Get sorted radii
    ind = np.argsort(r.flat)
    r_sorted = r.flat[ind]
    i_sorted = image.flat[ind]

    # Get the integer part of the radii (bin size = 1)
    r_int = np.around(r_sorted)

    # Find all pixels that fall within each radial bin.
    deltar = r_int[1:] - r_int[:-1]  # Assumes all radii represented
    rind = np.where(deltar)[0]       # location of changed radius
    nr = rind[1:] - rind[:-1]        # number of radius bin
    #print(nr)
    
    csim = np.cumsum(i_sorted, dtype=float)
    tbin = csim[rind[1:]] - csim[rind[:-1]]

    radial_prof = tbin / nr

    return radial_prof 

def fourd_roll_axis(stack):
    stack = np.rollaxis(np.rollaxis(stack, 2, 0), 3, 1)
    return stack

fd = DM.GetFrontImage()
print(fd)

origin0, scale0, unit0 = fd.GetDimensionCalibration(0, 0)
print(origin0, scale0, unit0)
origin1, scale1, unit1 = fd.GetDimensionCalibration(1, 0)
print(origin1, scale1, unit1)
origin2, scale2, unit2 = fd.GetDimensionCalibration(2, 0)
print(origin2, scale2, unit2)
origin3, scale3, unit3 = fd.GetDimensionCalibration(3, 0)
print(origin3, scale3, unit3)

print("loading 4D-STEM data")
stack_4d_cropped = fourd_roll_axis(fd.GetNumArray())
stack_4d_cropped = np.nan_to_num(stack_4d_cropped)
print(stack_4d_cropped.shape)
print(np.max(stack_4d_cropped))
print(np.min(stack_4d_cropped))
print(np.mean(stack_4d_cropped))

print("maximum-normalizing")
stack_4d_cropped = stack_4d_cropped / np.max(stack_4d_cropped)
print(np.max(stack_4d_cropped))
print(np.min(stack_4d_cropped))
print(np.mean(stack_4d_cropped))

pacbed = np.mean(stack_4d_cropped, axis=(0,1))

#find center position
q_text = """Select one option for finding the center position.
1: Gaussian fitting - PACBED
2: Center of mass - PACBED"""

q_check = int(input(q_text))

if q_check == 1:
	cb = int(input("size of the fitting box (data index): "))
	ct = gaussian_center(pacbed, cbox_edge=cb)
	print("center position")
	print(ct)

elif q_check == 2:
	Y, X = np.indices(pacbed.shape)
	com_y = np.sum(pacbed * Y) / np.sum(pacbed)
	com_x = np.sum(pacbed * X) / np.sum(pacbed)
	ct = [com_y, com_x]
		
	print("center position")
	print(ct)
		
else:
	print("*"*50)
	print("wrong input !")
	print("*"*50)
	exit()


# radial average of DPs (not variance, intensity direcltly, RDF?)
print("calculating the rotational average of each DP")
radial_avg_stack = []
len_profile = []
for i in range(stack_4d_cropped.shape[0]):
    for j in range(stack_4d_cropped.shape[1]):
        radial_temp = radial_average_with_center(stack_4d_cropped[i, j], center=ct)
        len_profile.append(len(radial_temp))
        radial_avg_stack.append(radial_temp)

print(np.unique(len_profile))
for i in range(len(radial_avg_stack)):
    radial_avg_stack[i] = radial_avg_stack[i][:np.min(len_profile)]

radial_avg_stack = np.asarray(radial_avg_stack).reshape(stack_4d_cropped.shape[0], stack_4d_cropped.shape[1], -1)
print(radial_avg_stack.shape)

radial_avg_stack = np.rollaxis(radial_avg_stack, 2, 0)
print(radial_avg_stack.shape)
ravg_dm = DM.CreateImage(radial_avg_stack.copy())
ravg_dm.SetName("rotational average spectrum image")
ravg_dm.SetDimensionCalibration(0, origin0, scale0, unit0, 0)
ravg_dm.SetDimensionCalibration(1, origin1, scale1, unit1, 0)
ravg_dm.SetDimensionCalibration(2, 0, scale2, unit2, 0)
ravg_dm.ShowImage()