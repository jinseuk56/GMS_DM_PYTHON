# Jinseok Ryu
# Electron Microscopy and Spectroscopy Lab.
# Seoul National University
# last update : 20210520
# save the front image (multi-dimensional) as a tiff stack file


# ********************************************************************************
print("Execute Python script in GMS 3")

import DigitalMicrograph as DM
import numpy as np
import tifffile
import sys
sys.argv.extend(['-a', ' '])

import tkinter.filedialog as tkf

print("Libraries have been imported completely")
# ********************************************************************************

def threed_roll_axis(img):
    stack = np.rollaxis(img, 2, 0)
    return stack
    
def fourd_roll_axis(stack):
    stack = np.rollaxis(np.rollaxis(stack, 2, 0), 3, 1)
    return stack

fd = DM.GetFrontImage()
print(fd)

fd = fd.GetNumArray()
print(fd.shape)

if len(fd.shape) == 4:
	fd = fourd_roll_axis(fd)
	
elif len(fd.shape) == 3:
	fd = threed_roll_axis(fd)

tifffile.imsave(tkf.asksaveasfilename(), fd)