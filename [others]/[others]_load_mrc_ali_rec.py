# Jinseok Ryu
# Electron Microscopy and Spectroscopy Lab.
# Seoul National University
# last update : 20210322
# load .mrc or .ali (tilt serise for electron tomography)
# -> it needs "mrcfile" and "Hyperspy" (Python packages)
# load .rec (reconstruction result from Inspect 3D)
# -> it requires the shape of the reconstructed box
# -> it is assumed that the data type of .rec files is 16-bit integer

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

if ( False == DM.IsScriptOnMainThread() ):
    print('MatplotLib scripts require to be run on the main thread.')
    exit()

file_adr = tkf.askopenfilename()
print(file_adr)

if file_adr[-3:] == "ali":
    img = hys.load(file_adr)
    print(img)
    
    origin0 = img.axes_manager[0].offset
    scale0 = img.axes_manager[0].scale
    unit0 = img.axes_manager[0].units
    print(origin0, scale0, unit0)

    origin1 = img.axes_manager[1].offset
    scale1 = img.axes_manager[1].scale
    unit1 = img.axes_manager[1].units
    print(origin1, scale1, unit1)

    origin2 = img.axes_manager[2].offset
    scale2 = img.axes_manager[2].scale
    unit2 = img.axes_manager[2].units
    print(origin2, scale2, unit2)
    print(type(origin2), type(scale2), type(unit2))

    data = img.data.copy()
    tilt_angles = img.original_metadata["fei header"]["a_tilt"][:len(data)]
    print("tilt angles (degree)")
    print(tilt_angles)
    tilt_rad = tilt_angles * np.pi / 180
    
    print(type(tilt_angles[0]), type(tilt_angles[1]-tilt_angles[0]), type("degree"))
    
    data_dm = DM.CreateImage(data.copy())
    data_dm.SetName("tilt series")
    
    data_dm.SetDimensionCalibration(0, origin1, scale1, unit1, 0)
    data_dm.SetDimensionCalibration(1, origin2, scale2, unit2, 0)
    data_dm.SetDimensionCalibration(2, float(tilt_angles[0]), float(tilt_angles[1]-tilt_angles[0]), "degree", 0)
    
    data_dm.ShowImage()
    
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
    data_dm.SetDimensionCalibration(2, float(tilt_angles[0]), float(tilt_angles[1]-tilt_angles[0]), "degree", 0)
    print(data_dm)
    data_dm.ShowImage()
    
elif file_adr[-3:] == "rec":
    s1 = eval(input("1st dimensions (positive integer): "))
    s2 = eval(input("2st dimensions (positive integer): "))
    s3 = eval(input("3st dimensions (positive integer): "))
    original_shape = (s1, s2, s3)
    stack_tmp = np.fromfile(file_adr, dtype=np.int16)
    orthoslices = stack_tmp[-np.prod(original_shape):].reshape(original_shape)
    data_dm = DM.CreateImage(orthoslices.copy())
    data_dm.SetName("orthoslices")
    print(data_dm)
    data_dm.ShowImage()