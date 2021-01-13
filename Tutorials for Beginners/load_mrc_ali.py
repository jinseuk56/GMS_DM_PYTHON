# Jinseok Ryu
# Electron Microscopy and Spectroscopy Lab.
# Seoul National University
# last update : 20210110
# load .mrc or .ali (tilt serise for electron tomography)

# ********************************************************************************
print("Execute Python script in GMS 3")

import sys
sys.argv.extend(['-a', ' '])
import tkinter.filedialog as tkf
import DigitalMicrograph as DM
import numpy as np
import hyperspy.api as hys
import mrcfile

print("Libraries have been imported completely")
# ********************************************************************************

file_adr = tkf.askopenfilename()
print(file_adr)

if file_adr[-3:] == "ali":
    img = hys.load(file_adr)
    print(img)
    data = img.data.copy()
    tilt_angles = img.original_metadata["fei header"]["a_tilt"][:len(data)]
    print("tilt angles (degree)")
    print(tilt_angles)
    tilt_rad = tilt_angles * np.pi / 180
    
elif file_adr[-3:] == "mrc":
    mrc_file = mrcfile.open(file_adr)
    print(mrc_file)
    data = mrc_file.data.copy()
    tilt_angles = []
    for arr in mrc_file.extended_header:
        tilt_angles.append(arr[10])
    tilt_angles = np.asarray(tilt_angles)
    print("tilt angles (degree)")
    print(tilt_angles)
    tilt_rad = tilt_angles * np.pi / 180
   
data_dm = DM.CreateImage(data.copy())
data_dm.SetName("tilt series")
print(data_dm)

data_dm.ShowImage()