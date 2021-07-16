# Jinseok Ryu
# Electron Microscopy and Spectroscopy Lab.
# Seoul National University
# 20210715
# remove the zero-loss peaks (ZLP) from a low-loss EELS spectrum image by deconvolving a reference ZLP
# the calibration values must be the same between the EELS-SI and reference ZLP
# it needs the FWHM of the reference ZLP

# ********************************************************************************
print("Execute Python script in GMS 3")

import numpy as np
import DigitalMicrograph as DM
import hyperspy.api as hys
import tifffile
import tkinter.filedialog as tkf

import sys
sys.argv.extend(['-a', ' '])

print("Libraries have been imported completely")
# ********************************************************************************

if (False == DM.IsScriptOnMainThread()):
    print('MatplotLib scripts require to be run on the main thread.')
    exit()

# ********************************************************************************

def threed_roll_axis(img):
    stack = np.rollaxis(img, 2, 0)
    return stack

si_adr = tkf.askopenfilename()
print(si_adr)

zlp_adr = tkf.askopenfilename()
print(zlp_adr)


si = hys.load(si_adr, signal_type="EELS")
ref_zlp = hys.load(zlp_adr, signal_type="EELS")


origin0 = si.axes_manager[0].offset
scale0 = si.axes_manager[0].scale
unit0 = si.axes_manager[0].units

origin1 = si.axes_manager[1].offset
scale1 = si.axes_manager[1].scale
unit1 = si.axes_manager[1].units

scale = ref_zlp.axes_manager[0].scale
offset = ref_zlp.axes_manager[0].offset + scale
n_channel = ref_zlp.axes_manager[0].size
print(scale, offset, n_channel)
e_range = np.arange(offset, n_channel*scale+offset, scale)
print(e_range[0], e_range[-1])

v_range = 1 / e_range

FWHM = eval(input("Input the FWHM of the reference ZLP(eV) : "))

FT_r = np.fft.fft(ref_zlp.data)
FT_g = np.max(ref_zlp.data)*np.exp(-np.pi**2*(FWHM/1.665)**2*v_range**2)
SSD = []

for spec in si.data.reshape(-1, n_channel):
    tmp_j = np.fft.fft(spec)
    tmp_k = np.multiply(tmp_j, FT_g) / FT_r
    SSD.append(np.abs(np.fft.ifft(tmp_k)))

SSD = np.asarray(SSD).reshape(si.data.shape[0], si.data.shape[1], si.data.shape[2])
SSD = threed_roll_axis(SSD)

data_dm = DM.CreateImage(SSD.copy())
data_dm.SetName("ZLP removed SI")
data_dm.SetDimensionCalibration(0, origin0, scale0, unit0, 0)
data_dm.SetDimensionCalibration(1, origin1, scale1, unit1, 0)
data_dm.SetDimensionCalibration(2, offset, scale, "eV", 0)

data_dm.ShowImage()