# Jinseok Ryu
# Electron Microscopy and Spectroscopy Lab.
# Seoul National University
# last update : 20210528
# load .blo (DigiSTAR)

# ********************************************************************************
print("Execute Python script in GMS 3")

import sys
sys.argv.extend(['-a', ' '])
import tkinter.filedialog as tkf
import DigitalMicrograph as DM
import numpy as np
import hyperspy.api as hys

print("Libraries have been imported completely")
# ********************************************************************************

if ( False == DM.IsScriptOnMainThread() ):
    print('MatplotLib scripts require to be run on the main thread.')
    exit()

def fourd_roll_axis(stack):
    stack = np.rollaxis(np.rollaxis(stack, 2, 0), 3, 1)
    return stack

img_adr = tkf.askopenfilename()
print(img_adr)

img = hys.load(img_adr)
print(img)

origin0 = img.axes_manager[0].offset
scale0 = img.axes_manager[0].scale
unit0 = img.axes_manager[0].units

origin1 = img.axes_manager[1].offset
scale1 = img.axes_manager[1].scale
unit1 = img.axes_manager[1].units

origin2 = img.axes_manager[2].offset
scale2 = img.axes_manager[2].scale
unit2 = img.axes_manager[2].units

origin3 = img.axes_manager[3].offset
scale3 = img.axes_manager[3].scale
unit3 = img.axes_manager[3].units

data = fourd_roll_axis(img.data)
print(data.shape)

data_dm = DM.CreateImage(data.copy())
data_dm.SetName("imported")

data_dm.SetDimensionCalibration(0, origin0, scale0, unit0, 0)
data_dm.SetDimensionCalibration(1, origin1, scale1, unit1, 0)
data_dm.SetDimensionCalibration(2, origin2, scale2, unit2, 0)
data_dm.SetDimensionCalibration(3, origin3, scale3, unit3, 0)

data_dm.ShowImage()