# Ingyu Yoo
# Electron Microscopy and Spectroscopy Lab.
# Seoul National University
# last update : 20210308
# DPC imaging for 4D-STEM Data
# This code is based on the following publication & code
# The core functions of GetDPC, "https://github.com/hachteja/GetDPC", were copied
# J.A. Hachtel, J.C. Idrobo, and M. Chi, Adv. Struct. Chem. Imag. 4 (2018)
# The jupyter notebook for demonstration of GetDPC were re-arranged so that it could be also applied to 4D-STEM data in GMS 3


#####################################################################
print("Execute Python script in GMS 3")

import DigitalMicrograph as DM
import typing
from scipy import optimize
import numpy as np
import sys
sys.argv.extend(['-a', ' '])
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.patches import Rectangle
from pylab import cm
from matplotlib.colors import hsv_to_rgb

print("Libraries have been imported completely")
#####################################################################

##refer to https://github.com/hachteja/GetDPC/blob/master/getdpc/GetDPC.py

def CalibrateRonchigram(dat4d: np.ndarray, conv: float = 32, t: float = 0.3) -> typing.Tuple[np.ndarray, np.ndarray, np.ndarray, float, float, np.ndarray, np.ndarray]:
    """Find true center of Ronchigram, and pixels/mrad calibration

    :param dat4d: 4D Dataset, 2-spatial, 2-diffraction dimensions
    :param conv: Convergence Angle of Electron Probe in mrad
    :param t: Threshhold for BF Disk (fraction of 1)
    :return: center, calibrations
    """
    R = np.average(dat4d, axis=(0, 1))
    Rn = (R - np.amin(R)) / np.ptp(R)
    BFdisk = np.ones(R.shape) * (Rn > t)
    absct = t * np.ptp(R)
    rxx, ryy = np.meshgrid(np.arange(0, Rn.shape[1]), np.arange(0, Rn.shape[0]))
    rcx, rcy = np.sum(BFdisk * rxx / np.sum(BFdisk)), np.sum(BFdisk * ryy / np.sum(BFdisk))
    edge = (np.sum(np.abs(np.gradient(BFdisk)), axis=0)) > t
    pixcal = np.average(np.sqrt((rxx - rcx) ** 2 + (ryy - rcy) ** 2)[edge]) / conv
    return R, rcx, rcy, pixcal, BFdisk, absct, edge
    
def GetDetectorImage(dat4d: np.ndarray, RCX: float, RCY: float, RCal: float, RI: float = 0, RO: float = 32) -> np.ndarray:
    """Reconstruct a detector image from the 4D Dataset

    :param dat4d: 4D Dataset, 2-spatial, 2-diffraction dimensions
    :param RCX: X Center of the Ronchigram (pixels)
    :param RCY: Y Center of the Ronchigram (pixels)
    :param RCal: Calibration of the Ronchigram (pixels/mrad)
    :param RI: Inner Radius for CoM Measurement (mrad)
    :param RO: Outer Radius for CoM Measurement (mrad)
    :return detector image as ndarray
    """
    X, Y = np.meshgrid((np.arange(0, dat4d.shape[3]) - RCX) / RCal, (np.arange(0, dat4d.shape[2]) - RCY) / RCal)
    return np.average(dat4d * ((X ** 2 + Y ** 2 >= RI ** 2) & (X ** 2 + Y ** 2 < RO ** 2)), axis=(2, 3))

def GetiCoM(dat4d: np.ndarray, RCX: float, RCY: float, RCal: float, RI: float = 0, RO: float = 32) -> typing.Tuple[np.ndarray, np.ndarray]:
    """Get Ronchigram Center of Mass Shifts from 4D Dataset

    :param dat4d: 4D Dataset, 2-spatial, 2-diffraction dimensions
    :param RCX: X Center of the Ronchigram (pixels)
    :param RCY: Y Center of the Ronchigram (pixels)
    :param RCal: Calibration of the Ronchigram (pixels/mrad)
    :param RI: Inner Radius for CoM Measurement (mrad)
    :param RO: Outer Radius for CoM Measurement (mrad)
    :return iCoM as ndarray
    """
    X, Y = np.meshgrid((np.arange(0, dat4d.shape[3]) - RCX) / RCal, (np.arange(0, dat4d.shape[2]) - RCY) / RCal)
    maskeddat4d = dat4d * ((X ** 2 + Y ** 2 >= RI ** 2) & (X ** 2 + Y ** 2 < RO ** 2))
    return np.average(maskeddat4d * X, axis=(2, 3)), np.average(maskeddat4d * Y, axis=(2, 3))

def GetPLRotation(dpcx: np.ndarray, dpcy: np.ndarray, *,  order: int = 3, outputall: bool = False) -> float:
    """Find Rotation from PL Lenses by minimizing curl/maximizing divergence of DPC data

    :param dpcx: X-Component of DPC Data (2D numpy array)
    :param dpcy: Y-Component of DPC Data (2D numpy array)
    :param order: Number of times to iterated calculation (int)
    :param outputall: Output Curl and Divergence curves for all guesses in separate array (bool)
    :return: The true PL Rotation value (Note: Can potentially be off by 180 degrees, determine by checking signs of charge/field/potential)
    """
    def DPC_ACD(dpcx,dpcy,tlow,thigh):        
        A,C,D=[],[],[]
        for t in np.linspace(tlow,thigh,10,endpoint=False):            
            rdpcx,rdpcy=dpcx*np.cos(t)-dpcy*np.sin(t),dpcx*np.sin(t)+dpcy*np.cos(t)        
            gXY,gXX=np.gradient(rdpcx);gYY,gYX=np.gradient(rdpcy)        
            C.append(np.std(gXY-gYX));D.append(np.std(gXX+gYY));A.append(t)
        R=np.average([A[np.argmin(C)],A[np.argmax(D)]])
        return R,A,C,D
    RotCalcs=[]
    RotCalcs.append(DPC_ACD(dpcx,dpcy,0,np.pi))
    for i in range(1,order): 
        RotCalcs.append(DPC_ACD(dpcx,dpcy,RotCalcs[i-1][0]-np.pi/(10**i),RotCalcs[i-1][0]+np.pi/(10**i)))
    if outputall: return RotCalcs
    else: return RotCalcs[-1][0]

def GetElectricFields(dpcx: np.ndarray, dpcy: np.ndarray, *, rotation: float = 0, LegPix: int = 301, LegRad: float = 0.85) -> typing.Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert dpcx and dpcy maps to to a color map where the color corresponds to the angle

    :param dpcx: X-Component of DPC Data (2D numpy array)
    :param dpcy: Y-Component of DPC Data (2D numpy array)
    :param rotation: Optional rotation radians
    :param LegPix: Number of Pixels in Color Wheel Legend
    :param LegRad: Radius of Color Wheel in Legend (0-1)
    :return: The electric fields as a 2D numpy array
    """
    EX = -dpcx
    EY = -dpcy
    rEX = EX * np.cos(rotation) - EY * np.sin(rotation)
    rEY = EX * np.sin(rotation) + EY * np.cos(rotation)

    EMag = np.sqrt(rEX ** 2 + rEY ** 2)

    XY = np.zeros(rEX.shape + (3,), dtype=float)
    M = np.amax(EMag)
    EMagScale = EMag / M
    for i in range(rEX.shape[0]):
        for j in range(rEX.shape[1]):
            XY[i, j] = np.angle(np.complex(rEX[i, j], rEY[i, j])) / (2 * np.pi) % 1, 1, EMagScale[i, j]
    EDir=hsv_to_rgb(XY)
    x, y = np.meshgrid(np.linspace(-1, 1, LegPix, endpoint=True), np.linspace(-1, 1, LegPix, endpoint=True))
    X, Y = x * (x ** 2 + y ** 2 < LegRad ** 2), y * (x ** 2 + y ** 2 < LegRad ** 2)
    XYLeg = np.zeros(X.shape + (3,), dtype=float)
    RI = np.sqrt(X ** 2 + Y ** 2) / np.amax(np.sqrt(X ** 2 + Y ** 2))
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            XYLeg[i, j] = np.angle(np.complex(X[i, j], Y[i, j])) / (2 * np.pi) % 1, 1, RI[i, j]
    EDirLeg=hsv_to_rgb(XYLeg)
    return EMag, EDir, EDirLeg

def GetChargeDensity(dpcx: np.ndarray, dpcy: np.ndarray, *, rotation: float = 0) -> np.ndarray:
    """Calculate Charge Density from the Divergence of the Ronchigram Shifts

    :param dpcx: X-Component of DPC Data (2D numpy array)
    :param dpcy: Y-Component of DPC Data (2D numpy array)
    :param rotation: Optional rotation radians
    :return: The charge density as a 2D numpy array
    """

    rdpcx = dpcx * np.cos(rotation) + dpcy * np.sin(rotation)
    rdpcy = -dpcx * np.sin(rotation) + dpcy * np.cos(rotation)
    gxx, gyy = np.gradient(rdpcx)[1], np.gradient(rdpcy)[0]

    return - gxx - gyy

def GetPotential(dpcx: np.ndarray, dpcy: np.ndarray, *, rotation: float = 0, hpass: float = 0, lpass: float = 0) -> np.ndarray:
    """Convert X and Y Shifts (E-Field Vector) Into Atomic Potential By Inverse Gradient

    Note: This method is vulnerable to edge induced artifacts that a small degree of high-pass filtering
    can clear up without significantly affecting the atomic-level contrast

    :param dpcx: X-Component of DPC Data (2D numpy array)
    :param dpcy: Y-Component of DPC Data (2D numpy array)
    :param rotation: Optional rotation radians
    :param hpass: Optional constant to provide variable high-pass filtering
    :param lpass: Optional constant to provide variable low-pass filtering
    :return: The potential as a 2D numpy array
    """

    rdpcx = dpcx * np.cos(rotation) + dpcy * np.sin(rotation)
    rdpcy = -dpcx * np.sin(rotation) + dpcy * np.cos(rotation)
    fCX = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(rdpcx)))
    fCY = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(rdpcy)))
    KX = fCX.shape[1]
    KY = fCY.shape[0]
    kxran = np.linspace(-1, 1, KX, endpoint=True)
    kyran = np.linspace(-1, 1, KY, endpoint=True)
    kx, ky = np.meshgrid(kxran, kyran)
    fCKX = fCX * kx
    fCKY = fCY * ky
    fnum = (fCKX + fCKY)
    fdenom = np.pi * 2 * (0 + 1j) * (hpass + (kx ** 2 + ky ** 2) + lpass * (kx ** 2 + ky ** 2) ** 2)
    fK = np.divide(fnum, fdenom)

    return np.real(np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(fK))))

### 0 . Get 4d data & virtual detector
image_4d = DM.GetFrontImage()
print("load 4D-STEM data")
dat4d = np.rollaxis(np.rollaxis(image_4d.GetNumArray(), 2, 0), 3, 1)
print(dat4d.shape)
print(np.max(dat4d))
print(np.min(dat4d))
print(np.mean(dat4d))

### 1. Calibrate Ronchigram to find calibration and Ronchigram center
dat4d_shape = dat4d.shape

dummy = 10.0

R, RonchiCenterX, RonchiCenterY, Ronchipixcal, BFdiskIm, absct, BFEdgeIm = CalibrateRonchigram(dat4d,conv=dummy,t=0.3)
print("center of the Ronchigram: %f (x), %f (y)"%(RonchiCenterX, RonchiCenterY))
Ronchipixcal = eval(input("Write a calibration value for the detector angle (mrad per pixel)"))
imcal = eval(input("Write a calibration value for STEM image (Angstrom per pixel)"))

showtype = int(input(""" Choose a result type
1) DM file
2) Matplotlib figure """))

print("R shape", R.shape)
f,a=plt.subplots(1,4,dpi=200)
a[0].imshow(dat4d[0,0],cmap=cm.nipy_spectral)
a[0].add_patch(Circle((RonchiCenterX,RonchiCenterY),radius=2,ec='w',fc='None',lw=1.))
a[0].set_title('Single Ronchigram',fontsize=6)
a[0].add_patch(Rectangle((4,dat4d.shape[2]-7),50/Ronchipixcal,3,fc='w',ec='None'))
a[0].text(4+25/Ronchipixcal,dat4d.shape[2]-7,'50 mrad',fontweight='bold',color='w',fontsize=6,ha='center',va='bottom')

a[1].imshow(np.average(dat4d,axis=(0,1)),cmap=cm.nipy_spectral)
a[1].add_patch(Circle((RonchiCenterX,RonchiCenterY),radius=2,ec='w',fc='None',lw=1.))
a[1].set_title('Average Ronchigram',fontsize=6)
a[1].add_patch(Rectangle((4,dat4d.shape[2]-7),50/Ronchipixcal,3,fc='w',ec='None'))

a[2].imshow(BFdiskIm,cmap=cm.nipy_spectral)
a[2].add_patch(Circle((RonchiCenterX,RonchiCenterY),radius=2,ec='w',fc='None',lw=1.))
a[2].set_title('Average BF Disk',fontsize=6)
a[2].add_patch(Rectangle((4,dat4d.shape[2]-7),50/Ronchipixcal,3,fc='w',ec='None'))

a[3].imshow(BFEdgeIm,cmap=cm.nipy_spectral)
a[3].add_patch(Circle((RonchiCenterX,RonchiCenterY),radius=2,ec='w',fc='None',lw=1.))
a[3].set_title('Average BF Disk Edge',fontsize=6)
a[3].add_patch(Rectangle((4,dat4d.shape[2]-7),50/Ronchipixcal,3,fc='w',ec='None'))
plt.setp(a, xticks=[],yticks=[])


"""
### 2. Use 4D Dataset to reconstruct images from arbitrary detector
BF_outer = 10.0 #mrad
BF = GetDetectorImage(dat4d,RonchiCenterX,RonchiCenterY,Ronchipixcal,0.,BF_outer)
print("BF shape",BF.shape)
inner_ind = int(input("\n
To get an as-acquired stem image profile, inner and outer angles for a virtual detector are needed.\n
First, write the index of the inner angle of the virtual detector (positive integer): \n
"))
outer_ind = int(input("Next, write the index of the outer angle of the virtual detector (positive integer):"))
ADF_STEM = GetDetectorImage(dat4d,RonchiCenterX,RonchiCenterY,Ronchipixcal,inner_ind*Ronchipixcal,outer_ind*Ronchipixcal)
if showtype == 2 :
    f2,a2=plt.subplots(1,2)
    a2[0].imshow(ADF_STEM,cmap=cm.inferno)
    a2[0].add_patch(Rectangle((4,dat4d.shape[0]-7),0.5/imcal,3,fc='w',ec='None'))
    a2[0].text(4+0.25/imcal,dat4d.shape[0]-7,'5 nm',fontweight='bold',color='w',fontsize=5,ha='center',va='bottom')
    a2[0].set_title('As-Acquired STEM ({}-{} mrad)'.format(inner_ind*Ronchipixcal, outer_ind*Ronchipixcal),fontsize=6)
    a2[1].imshow(BF,cmap=cm.inferno)
    a2[1].add_patch(Rectangle((4,dat4d.shape[0]-7),0.5/imcal,3,fc='w',ec='None'))
    a2[1].set_title('Reconstructed BF (0-%d mrad)'%BF_outer,fontsize=8)
    plt.setp(a2, xticks=[],yticks=[])
elif showtype == 1 :
    im2_1 = DM.CreateImage(ADF_STEM.copy())
    im2_1.SetName("As-Acquired STEM ({}-{} mrad)".format(inner_ind*Ronchipixcal, outer_ind*Ronchipixcal))
    im2_2 = DM.CreateImage(BF.copy())
    im2_2.SetName('Reconstructed BF (0-%d mrad)'%BF_outer)
    im2_1.ShowImage()
    im2_2.ShowImage()
"""

### 3. Calculate Center of Mass Shifts (No Rotation)

CoMX,CoMY=GetiCoM(dat4d,RonchiCenterX,RonchiCenterY,Ronchipixcal)
Checkbox2 = input(" Do you want to apply PLRotation? Y/N")

if Checkbox2 == "Y":
    RotationCalcs=GetPLRotation(CoMX,CoMY,order=4,outputall=True)
    PLRotation=RotationCalcs[-1][0]
    rCoMX, rCoMY = CoMX * np.cos(PLRotation) - CoMY * np.sin(PLRotation), CoMX * np.sin(PLRotation) + CoMY * np.cos(PLRotation)
    Checkpoint = input(" Do you want to compare Center of Mass with PLRotation and without PLRotation? Y/N")

    if Checkpoint == "Y":
        if showtype == 2:
            f3,a3=plt.subplots(2,2,dpi=200)
            a3[0,0].imshow(CoMX)
            a3[0,0].add_patch(Rectangle((4,dat4d.shape[0]-7),0.5/imcal,3,fc='w',ec='None'))
            a3[0,0].text(4+0.25/imcal,dat4d.shape[0]-7,'5 A',fontweight='bold',color='w',fontsize=5,ha='center',va='bottom')
            a3[0,0].set_title('CoM Shifts-X (No PL Rotation)',fontsize=6)
            a3[0,1].imshow(CoMY)
            a3[0,1].add_patch(Rectangle((4,dat4d.shape[0]-7),0.5/imcal,3,fc='w',ec='None'))
            a3[0,1].text(4+0.25/imcal,dat4d.shape[0]-7,'5 A',fontweight='bold',color='w',fontsize=5,ha='center',va='bottom')
            a3[0,1].set_title('CoM Shifts-Y (No PL Rotation)',fontsize=6)
            a3[1,0].imshow(rCoMX)
            a3[1,0].add_patch(Rectangle((4,dat4d.shape[0]-7),0.5/imcal,3,fc='w',ec='None'))
            a3[1,0].text(4+0.25/imcal,dat4d.shape[0]-7,'5 A',fontweight='bold',color='w',fontsize=5,ha='center',va='bottom')
            a3[1,0].set_title('CoM Shifts-X (w/ PL Rotation)',fontsize=6)
            a3[1,1].imshow(rCoMY)
            a3[1,1].add_patch(Rectangle((4,dat4d.shape[0]-7),0.5/imcal,3,fc='w',ec='None'))
            a3[1,1].text(4+0.25/imcal,dat4d.shape[0]-7,'5 A',fontweight='bold',color='w',fontsize=5,ha='center',va='bottom')
            a3[1,1].set_title('CoM Shifts-Y (w/ PL Rotation)',fontsize=6)
            plt.setp(a3, xticks=[],yticks=[])
        elif showtype == 1:
            im3_1 = DM.CreateImage(CoMX.copy())
            im3_1.SetName("CoM Shifts-X (No PL Rotation)")
            im3_1.ShowImage()
            im3_2 = DM.CreateImage(CoMY.copy())
            im3_2.SetName("CoM Shifts-Y (No PL Rotation)")
            im3_2.ShowImage()
            im3_3 = DM.CreateImage(rCoMX.copy())
            im3_3.SetName("CoM Shifts-X (w/ PL Rotation)")
            im3_3.ShowImage()
            im3_4 = DM.CreateImage(rCoMY.copy())
            im3_4.SetName("CoM Shifts-Y (w/ PL Rotation)")
            im3_4.ShowImage()	

elif Checkbox2 == "N":
    rCoMX, rCoMY = CoMX, CoMY
    
    
### 4. Calculate Electric Field Magnitude and Direction from Rotated CoM Shifts
### NOTE: Without accounting for PL rotation, E-Field does not point radially away from nuclei.
### This is not physical so the PL rotation calculation is required
### 5. Calculate PL Rotation & Calculate Center of Mass Shifts (With Calculated Rotation)

checkbox3 = list(input("""
Write one or many numbers (separting them with commas) described below :
1) Electric field / 2) Charge density / 3) Potential
"""))
while "," in checkbox3:
    checkbox3.remove(",")
                    
### 6. Calculate Electric Field Magnitude and Direction from Rotated CoM Shifts

if "1" in checkbox3:
    EIm,EDirIm,EDirLeg=GetElectricFields(rCoMX,rCoMY)
    
    if showtype ==2:
        f6,a6=plt.subplots(1,3,dpi=200)
        a6[0].imshow(EIm)
        a6[0].add_patch(Rectangle((4,dat4d.shape[0]-7),0.5/imcal,3,fc='w',ec='None'))
        a6[0].text(4+0.25/imcal,dat4d.shape[0]-7,'5 A',fontweight='bold',color='w',fontsize=5,ha='center',va='bottom')
        a6[0].set_title('Electric Field Magnitude',fontsize=6)
        a6[1].imshow(EDirIm)
        a6[1].add_patch(Rectangle((4,dat4d.shape[0]-7),0.5/imcal,3,fc='w',ec='None'))
        a6[1].text(4+0.25/imcal,dat4d.shape[0]-7,'5 A',fontweight='bold',color='w',fontsize=5,ha='center',va='bottom')
        a6[1].set_title('Electric Field Directions',fontsize=6)
        a6[2].imshow(EDirLeg)
        a6[2].set_title('Field Directions Legend',fontsize=6)
        plt.setp(a6, xticks=[],yticks=[])
    elif showtype == 1:
        f6,a6=plt.subplots(1,2,dpi=200)
        a6[0].imshow(EDirIm)
        a6[0].add_patch(Rectangle((4,dat4d.shape[0]-7),0.5/imcal,3,fc='w',ec='None'))
        a6[0].text(4+0.25/imcal,dat4d.shape[0]-7,'5 A',fontweight='bold',color='w',fontsize=5,ha='center',va='bottom')
        a6[0].set_title('Electric Field Directions',fontsize=6)
        a6[1].imshow(EDirLeg)
        a6[1].set_title('Field Directions Legend',fontsize=6)
        plt.setp(a6, xticks=[],yticks=[])
        im5_1 = DM.CreateImage(EIm.copy())
        im5_1.SetName('Electric Field Magnitude')
        im5_1.ShowImage()
        
        

### 7. Calculate Charge Density from Divergence of CoM Shifts

if "2" in checkbox3:
    RhoIm = GetChargeDensity(rCoMX, rCoMY)
    if showtype == 2:
        f7,a7=plt.subplots()
        a7.imshow(RhoIm,cmap=cm.seismic,vmin=-np.amax(np.abs(RhoIm)),vmax=np.amax(np.abs(RhoIm)))
        a7.add_patch(Rectangle((4,dat4d.shape[0]-7),0.5/imcal,3,fc='w',ec='None'))
        a7.text(4+0.25/imcal,dat4d.shape[0]-7,'5 A',fontweight='bold',color='w',fontsize=5,ha='center',va='bottom')
        a7.set_title('Charge Density',fontsize=8)
        plt.setp(a7, xticks=[],yticks=[])
        
    elif showtype == 1:
        im7 = DM.CreateImage(RhoIm.copy())
        im7.SetName("Charge Density")
        im7.ShowImage()
    

### 8. Calculate Potential from Inverse Gradient of CoM Shifts
### Add High Pass Filtering to remove Edge Effects

if "3" in checkbox3:
    hpassnum1, hpassnum2= eval(input("""
3) Calculate Potential
Write two positive real numbers between 0 and 1 with a comma for high-pass filtering"""))
    VImhp1 =GetPotential(rCoMX,rCoMY, hpass=hpassnum1)    
    VImhp2 =GetPotential(rCoMX,rCoMY,hpass=hpassnum2)
    if showtype == 2:
        f8,a8=plt.subplots(1,2,dpi=200)
        a8[0].imshow(VImhp1,cmap=cm.hot)
        a8[0].add_patch(Rectangle((4,dat4d.shape[0]-7),0.5/imcal,3,fc='w',ec='None'))
        a8[0].text(4+0.25/imcal,dat4d.shape[0]-7,'5 A',fontweight='bold',color='w',fontsize=5,ha='center',va='bottom')
        a8[0].set_title('Atomic Potential (w/ High Pass {})'.format(hpassnum1),fontsize=6)
        a8[1].imshow(VImhp2,cmap=cm.hot)
        a8[1].add_patch(Rectangle((4,dat4d.shape[0]-7),0.5/imcal,3,fc='w',ec='None'))
        a8[1].text(4+0.25/imcal,dat4d.shape[0]-7,'5 A',fontweight='bold',color='w',fontsize=5,ha='center',va='bottom')
        a8[1].set_title('Atomic Potential (w/ High Pass {})'.format(hpassnum2),fontsize=6)
        plt.setp(a8, xticks=[],yticks=[])
        
    elif showtype == 1:
        im8_1 = DM.CreateImage(VImhp1.copy())
        im8_1.SetName('Atomic Potential (w/ High Pass {})'.format(hpassnum1))
        im8_2 = DM.CreateImage(VImhp2.copy())
        im8_2.SetName('Atomic Potential (w/ High Pass {})'.format(hpassnum2))
        im8_1.ShowImage()
        im8_2.ShowImage()

plt.show()