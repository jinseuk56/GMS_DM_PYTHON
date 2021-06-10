# Jinseok Ryu
# Electron Microscopy and Spectroscopy Lab.
# Seoul National University
# 20210520
# load velox STEM-EDS-SI (spectrum image)

import sys
sys.argv.extend(['-a', ' '])
import DigitalMicrograph as DM

import numpy as np
import tkinter.filedialog as tkf
import hyperspy.api as hys


def threed_roll_axis(img):
    stack = np.rollaxis(img, 2, 0)
    return stack
   
eds_adr = tkf.askopenfilename()
print(eds_adr)

eds = hys.load(eds_adr)[-1]
print(eds)

origin0 = eds.axes_manager[0].offset
scale0 = eds.axes_manager[0].scale
unit0 = eds.axes_manager[0].units

origin1 = eds.axes_manager[1].offset
scale1 = eds.axes_manager[1].scale
unit1 = eds.axes_manager[1].units

step = eds.axes_manager[2].scale
left = eds.axes_manager[2].offset
e_unit = eds.axes_manager[2].units

dm_out = eds.data.copy()
dm_out = dm_out.astype(np.uint8)
dm_out = threed_roll_axis(dm_out)
data_dm = DM.CreateImage(dm_out.copy())

data_dm.SetDimensionCalibration(0, origin0, scale0, unit0, 0)
data_dm.SetDimensionCalibration(1, origin1, scale1, unit1, 0)
data_dm.SetDimensionCalibration(2, left, step, e_unit, 0)

data_dm.ShowImage()
