# Jinseok Ryu
# Electron Microscopy and Spectroscopy Lab.
# Seoul National University
# last update : 20210323
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

data = fourd_roll_axis(img.data)
print(data.shape)

data_dm = DM.CreateImage(data.copy())
data_dm.SetName("imported")
data_dm.ShowImage()