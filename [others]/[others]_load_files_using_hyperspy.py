# Jinseok Ryu
# Electron Microscopy and Spectroscopy Lab.
# Seoul National University
# 20230710
# Load files using the io module of the Hyperspy package
# it supports the file formats that the Hyperspy package does (but, it may not work properly)
# You can import one or many files
# The imported data will be returned as DM file(s) on the activated workspace
# The memory problem (lack of available memory) may arise when you import a STEM-EDS spectrum image (or a large multi-dimensional data)

import sys
sys.argv.extend(['-a', ' '])
import DigitalMicrograph as DM

import os
import numpy as np
import tkinter.filedialog as tkf
import hyperspy
import hyperspy.api as hys

if ( False == DM.IsScriptOnMainThread() ):
    print('MatplotLib scripts require to be run on the main thread.')
    exit()

def threed_roll_axis(img):
    stack = np.rollaxis(img, 2, 0)
    return stack
   
def fourd_roll_axis(stack):
    stack = np.rollaxis(np.rollaxis(stack, 2, 0), 3, 1)
    return stack
    
def transform_to_DM(img, file_adr, datatype=False, roll_axis=True):
    n_dim = len(img.data.shape)
    calibration_info = []
    for i in range(n_dim):
        origin = img.axes_manager[i].offset
        scale = img.axes_manager[i].scale
        unit = img.axes_manager[i].units
        if not isinstance(unit, str):
            unit = "undefined"
        calibration_info.append([origin, scale, unit.replace(" ", "")])
    
    dm_out = img.data.copy()
    
    if datatype:
        dm_out = dm_out.astype(datatype)
        
    if roll_axis:
        if n_dim == 3:
            dm_out = threed_roll_axis(dm_out)
            
        if n_dim == 4:
            dm_out = fourd_roll_axis(dm_out)
            
    dm_out = DM.CreateImage(dm_out.copy())
    
    for i in range(n_dim):
        dm_out.SetDimensionCalibration(i, calibration_info[i][0], calibration_info[i][1], calibration_info[i][2], 0)
    
    try:
        dm_out.SetName(img.metadata.General.original_filename[:-4]+"_"+img.metadata.General.title)
    except:
        dm_out.SetName(os.path.basename(file_adr).split("\\").split(".")[0])
    dm_out.ShowImage()
    
file_adr = tkf.askopenfilenames()
print(file_adr)

for i in range(len(file_adr)):
    data_loaded = hys.load(file_adr[i])
    print(data_loaded)
    if isinstance(data_loaded, list):
        for j in range(len(data_loaded)):
            if data_loaded[j].metadata.General.title == "EDS" and len(data_loaded[j].data.shape) == 3:
                transform_to_DM(data_loaded[j], file_adr[i], datatype=np.uint8, roll_axis=True)
            else:
                transform_to_DM(data_loaded[j], file_adr[i], roll_axis=False)
            
    else:
        if data_loaded.metadata.General.title == "EDS" and len(data_loaded.data.shape) == 3:
            transform_to_DM(data_loaded, file_adr[i], datatype=np.uint8, roll_axis=True)
        else:
            transform_to_DM(data_loaded, file_adr[i], roll_axis=False)
