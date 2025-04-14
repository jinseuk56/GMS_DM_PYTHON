# Jinseok Ryu
# Electron Microscopy and Spectroscopy Lab.
# Seoul National University
# 20230710
# load one or many tiff files

import os
import sys
sys.argv.extend(['-a', ' '])
import DigitalMicrograph as DM
import numpy as np
import tifffile
import tkinter.filedialog as tkf

def threed_roll_axis(img):
    stack = np.rollaxis(img, 2, 0)
    return stack
    
def fourd_roll_axis(stack):
    stack = np.rollaxis(np.rollaxis(stack, 2, 0), 3, 1)
    return stack

img_adr = tkf.askopenfilenames()
print("number of the selected files: ", len(img_adr))
print(*img_adr, sep="\n")

q_check = input("Do you want to reshape the data? (Y/N): ")
if q_check == "Y":
    q_dim = input("Write the data shape in parenthesis, e.g., (256, 256, 128, 128)")
    data_shape = eval(q_dim)
    print("final data shape :", data_shape)


for i in range(len(img_adr)):
    adr = img_adr[i]
    dataname = os.path.basename(adr).split('.')[0]
    if adr[-3:] == "tif" or adr[-4:]=="tiff":
        data = tifffile.imread(adr)
        if q_check == "Y":
            data = np.reshape(data, data_shape)
        print(data.shape)

    else:
        print("wrong input !")
        exit()

    if len(data.shape) == 4:
        data = fourd_roll_axis(data)

    elif len(data.shape) == 3:
        data = threed_roll_axis(data)


    data_dm = DM.CreateImage(data.copy())
    data_dm.SetName(dataname)
    data_dm.ShowImage()
