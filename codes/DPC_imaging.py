# Ingyu Yoo, Jinseok Ryu
# Electron Microscopy and Spectroscopy Lab.
# Seoul National University
# last update : 20210630
# differential phase contrast imaging for 4D-STEM data

import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors

class dpc_python():
    
    def __init__(self, f_stack, ang_per_pixel, mrad_per_pixel):
        self.original_stack = f_stack
        self.original_shape = f_stack.shape
        self.original_pacbed = np.mean(self.original_stack, axis=(0, 1))
        self.ang_per_pixel = ang_per_pixel
        self.mrad_per_pixel = mrad_per_pixel
        
        print("the shape of the data =", self.original_shape)
    
    def find_center(self, com=True, gaussian=False):
        
        if com and gaussian:
            print("Warning! Choose only one option to find the center")
        
        if not com and not gaussian:
            print("Warning! Choose at least one option to find the center")
        
        Y, X = np.indices(self.original_pacbed.shape)
        com_y = np.sum(self.original_pacbed * Y) / np.sum(self.original_pacbed)
        com_x = np.sum(self.original_pacbed * X) / np.sum(self.original_pacbed)
        self.com_ct = [com_y, com_x]
        
        (_, center_y, center_x, _, _) = fitgaussian(self.original_pacbed)
        self.gauss_ct = [center_y, center_x]
        
        if com:
            self.ct=self.com_ct
        
        else:
            self.ct=self.gauss_ct
        
    def disk_extract(self, buffer_size=0):
        grad = np.gradient(self.original_pacbed)
        grad_map = grad[0]**2 + grad[1] **2
        grad_map = grad_map / np.max(grad_map)
        
        max_ind = np.unravel_index(np.argmax(grad_map, axis=None), grad_map.shape)
        self.least_R = ((max_ind[0]-self.ct[0])**2 + (max_ind[1]-self.ct[1])**2)**(1/2)
        
        print("radius of the BF disk = %.2f mrad"%(self.mrad_per_pixel*self.least_R))
        
        self.ct_ind  = np.around(self.ct).astype(int)
        self.cropped_size = np.around(self.least_R + buffer_size).astype(int)
        
        print("radius of the RoI = %.2f mrad"%(self.mrad_per_pixel*self.cropped_size))
        
        if self.cropped_size > self.ct_ind[0] or self.cropped_size > self.ct_ind[1]:
            self.cropped_size = np.min(ct_ind)
        
        self.c_ct = [self.cropped_size, self.cropped_size]
        self.c_stack = self.original_stack[:, :, self.ct_ind[0]-self.cropped_size:self.ct_ind[0]+self.cropped_size+1, 
                               self.ct_ind[1]-self.cropped_size:self.ct_ind[1]+self.cropped_size+1].copy()
        self.c_shape = self.c_stack.shape
        self.c_pacbed = np.mean(self.c_stack, axis=(0, 1))
        
    def virtual_stem(self, BF, ADF):
        self.BF_detector = radial_indices(self.original_pacbed.shape, BF, self.mrad_per_pixel, center=self.ct)
        self.BF_stem = np.sum(np.multiply(self.original_stack, self.BF_detector), axis=(2, 3))
        
        self.ADF_detector = radial_indices(self.original_pacbed.shape, ADF, self.mrad_per_pixel, center=self.ct)
        self.ADF_stem = np.sum(np.multiply(self.original_stack, self.ADF_detector), axis=(2, 3))
        
    def DPC(self, correct_rotation=True, n_theta=100, hpass=0.0, lpass=0.0):
        """
        Hachtel, J.A., J.C. Idrobo, and M. Chi, Adv Struct Chem Imaging, 2018. 4(1): p. 10. (https://github.com/hachteja/GetDPC)
        Lazic, I., E.G.T. Bosch, and S. Lazar, Ultramicroscopy, 2016. 160: p. 265-280.
        Savitzky, B.H., et al., arXiv preprint arXiv:2003.09523, 2020. (https://github.com/py4dstem/py4DSTEM)
        """
        
        Y, X = np.indices(test.c_pacbed.shape)
        self.ysh = np.sum(test.c_stack * Y, axis=(2, 3)) / np.sum(test.c_stack, axis=(2, 3)) - test.c_ct[0]
        self.xsh = np.sum(test.c_stack * X, axis=(2, 3)) / np.sum(test.c_stack, axis=(2, 3)) - test.c_ct[1]
        
        self.ysh -= np.mean(self.ysh)
        self.xsh -= np.mean(self.xsh)
        
        if correct_rotation:
            theta = np.linspace(-np.pi/2, np.pi/2, n_theta, endpoint=True)
            self.div = []
            self.curl = []
            for t in theta:
                r_ysh = self.xsh * np.sin(t) + self.ysh * np.cos(t)
                r_xsh = self.xsh * np.cos(t) - self.ysh * np.sin(t)

                gyy, gyx = np.gradient(r_ysh)
                gxy, gxx = np.gradient(r_xsh)
                shift_divergence = gyy + gxx
                shift_curl = gyx - gxy

                self.div.append(np.mean(shift_divergence**2))
                self.curl.append(np.mean(shift_curl**2))
                
            self.c_theta = theta[np.argmin(self.curl)]
            tmp_ysh = self.xsh * np.sin(self.c_theta) + self.ysh * np.cos(self.c_theta)
            tmp_xsh = self.xsh * np.cos(self.c_theta) - self.ysh * np.sin(self.c_theta)
            
            self.ysh = tmp_ysh
            self.xsh = tmp_xsh
            
        self.E_mag = np.sqrt(self.ysh**2 + self.xsh**2)
        self.E_field_y = -self.ysh / np.max(self.E_mag)
        self.E_field_x = -self.xsh / np.max(self.E_mag)
        
        self.charge_density = np.gradient(self.E_field_y)[0] + np.gradient(self.E_field_x)[1]
        
        self.potential = get_icom(self.ysh, self.xsh, hpass, lpass)

#####################################################################################
# https://scipy-cookbook.readthedocs.io/items/FittingData.html
def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*np.exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    X, Y = np.indices(data.shape) # row, col
    x = (X*data).sum()/total # row
    y = (Y*data).sum()/total # col
    col = data[:, int(y)]
    width_x = np.sqrt(np.abs((np.arange(col.size)-y)**2*col).sum()/col.sum()) # row
    row = data[int(x), :]
    width_y = np.sqrt(np.abs((np.arange(row.size)-x)**2*row).sum()/row.sum()) # col
    height = data.max()
    return height, x, y, width_x, width_y

def fitgaussian(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    errorfunction = lambda p: np.ravel(gaussian(*p)(*np.indices(data.shape)) -
                                 data)
    p, success = optimize.leastsq(errorfunction, params)
    return p
#####################################################################################

def radial_indices(shape, radial_range, scale, center=None):
    y, x = np.indices(shape)
    if not center:
        center = np.array([(y.max()-y.min())/2.0, (x.max()-x.min())/2.0])
    
    r = np.hypot(y - center[0], x - center[1]) * scale
    ri = np.ones(r.shape)
    
    if len(np.unique(radial_range)) > 1:
        ri[np.where(r < radial_range[0])] = 0
        ri[np.where(r > radial_range[1])] = 0
        
    else:
        r = np.round(r)
        ri[np.where(r != round(radial_range[0]))] = 0
    
    return ri

def get_icom(ysh, xsh, hpass=0, lpass=0):
    
    FT_ysh = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(ysh)))
    FT_xsh = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(xsh)))
    
    ky = np.fft.fftshift(np.fft.fftfreq(FT_ysh.shape[0])).reshape(-1, 1)
    kx = np.fft.fftshift(np.fft.fftfreq(FT_xsh.shape[1])).reshape(1, -1)

    k2 = ky**2 + kx**2
    zero_ind = np.where(k2 == 0.0)
    k2[zero_ind] = 1.0

    FT_phase = (FT_ysh*ky + FT_xsh*kx) / (2*np.pi*1j*(hpass+k2) + lpass*k2)
    FT_phase[zero_ind] = 0.0

    Iicom = np.real(np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(FT_phase))))
    
    return Iicom

# Run with GMS 3
print("Execute Python script in GMS 3")
import DigitalMicrograph as DM
import sys
sys.argv.extend(['-a', ' '])

image_4d = DM.GetFrontImage()
print("load 4D-STEM data")
origin0, scale0, unit0 = image_4d.GetDimensionCalibration(0, 0)
print(origin0, scale0, unit0)
origin1, scale1, unit1 = image_4d.GetDimensionCalibration(1, 0)
print(origin1, scale1, unit1)
origin2, scale2, unit2 = image_4d.GetDimensionCalibration(2, 0)
print(origin2, scale2, unit2)
origin3, scale3, unit3 = image_4d.GetDimensionCalibration(3, 0)
print(origin3, scale3, unit3)

image_4d = np.nan_to_num(image_4d)
f_stack = np.rollaxis(np.rollaxis(image_4d.GetNumArray(), 2, 0), 3, 1)
print(f_stack.shape)

# Should know what ang_per_pixel & mrad_per_pixel is.
test = dpc_python(f_stack, scale0, scale2)

# Create number for DPC analysis
while True:
    center_num = int(input("""Select one option for finding the center position.
    1) Center of Mass-PACBED
    2) Gaussian fitting-PACBED"""))
    crop_num = int(input("""Do you want to crop diffraction data?
    1) Yes 2) No"""))
    rot_num = int(input("""Do you want to find PL Rotation
    1) Yes 2) No"""))
    plot_num = int(input("""Select one option for result data type.
    1) DigitalMicrograph file
    2) Matplotlib subplot file"""))
    if center_num in [1,2] and crop_num in [1,2] and rot_num in [1,2]:
        break
    else:
        print("Warning! Choose only between 1 and 2")

# Create number for image plot
x = True
while x == True:
    img_num = list(input("""Choose between four options. (Multiple options are allowed)
1) Get PACBED (with center indication)
2) Virtual imaging (BF & ADF image)
3) Electric field
4) Charge density & Potential"""))
    while "," in img_num:
        img_num.remove(",")
    x = False
    for i in range(len(img_num)):
        if int(img_num[i]) not in [1,2,3,4]:
            print("Warning! Choose only option 1 to 4")
            x = True

if "4" in img_num:
    high_num = float(input("""High pass filter : (Choose the number between 0 to 1)
    If you don't need high pass filter, choose 0"""))
    low_num = float(input("""Low pass filter : (Choose the number between 0 to 1)
    If you don't need low pass filter, choose 0"""))
else:
    high_num, low_num = 0.0, 0.0

if center_num ==1:
    test.find_center(com=True, gaussian=False)
elif center_num ==2:
    test.find_center(com=False, gaussian=True)

buffer_size = 15
test.disk_extract(buffer_size)

if crop_num ==1:
    print("Crop diffraction 2D with buffer_size:{} pixel".format(buffer_size))
    print("Cropped PACBED shape: {}".format(test.c_stack.shape))
    print("Least radius: {}".format(test.least_R))
    
if crop_num ==2:
    test.c_ct = test.ct
    test.c_stack = test.original_stack
    test.c_pacbed = test.original_pacbed
    test.c_shape =  test.original_stack.shape

BF_det = [0, 25]
ADF_det = [50, 80]
test.virtual_stem(BF_det, ADF_det)

if rot_num == 1:
    correct_rotation = True
if rot_num ==2:
    correct_rotation = False
test.DPC(correct_rotation, n_theta=100, hpass= high_num, lpass= low_num)
if rot_num ==1:
    print("Optimized angle = {} degree".format(test.c_theta*180/np.pi))

if plot_num ==1:
    if "1" in img_num:
        im1 = DM.CreateImage(test.c_pacbed.copy())
        im1.SetName("PACBED")
        im1.ShowImage()
        # Center Indication test.c_ct[1], test.c_ct[0]

    if "2" in img_num:
        im2_1 = DM.CreateImage(test.BF_stem.copy())
        im2_1.SetName("BF-STEM image")
        im2_1.ShowImage()
        im2_2 = DM.CreateImage(test.ADF_stem.copy())
        im2_2.SetName("ADF-STEM image")
        im2_2.ShowImage()

    if "3" in img_num:
        im3_1 = DM.CreateImage(test.ADF_stem.copy())
        im3_1.SetName("ADF-STEM image")
        im3_1.ShowImage()
        im3_2 = DM.CreateImage(test.E_field_y.copy())
        im3_2.SetName("Electric field y_component")
        im3_2.ShowImage()
        im3_3 = DM.CreateImage(test.E_field_x.copy())
        im3_3.SetName("Electric field x_component")
        im3_3.ShowImage()
        im3_4 = DM.CreateImage(test.E_mag.copy())
        im3_4.SetName("Electric field magnitude")
        im3_4.ShowImage()

    if "4" in img_num:
        
        dir_num = int(input("""Electric field direction plot is only available in matplotlib. Proceed?
        Select one option.
        1) Yes (Get matplotlib subplot)
        2) No (Skip electric field direction plot)"""))
        if dir_num == 1:
            RY, RX = np.indices(test.c_shape[:2])
            fig, ax = plt.subplots(1, 1, figsize=(5, 5))
            ax.imshow(test.ADF_stem, cmap="gray")
            ax.quiver(RX.flatten(), RY.flatten(), test.E_field_x.flatten(), test.E_field_y.flatten(), color=cm.jet(mcolors.Normalize()(test.E_mag.flatten())))
            ax.set_title("Electric field direction")
            ax.axis("off")
            fig.tight_layout()
            plt.show()
        
        im4_1 = DM.CreateImage(test.charge_density.copy())
        im4_1.SetName("Charge density")
        im4_1.ShowImage()
        im4_2 = DM.CreateImage(test.potential.copy())
        im4_2.SetName("Potential")
        im4_2.ShowImage()

if plot_num == 2:
    if "1" in img_num:
        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        ax.imshow(test.c_pacbed, cmap="jet")
        ax.set_title("PACBED")
        ax.scatter(test.c_ct[1], test.c_ct[0], s=15, c="k")
        ax.axis("off")
        fig.tight_layout()

    if "2" in img_num:
        fig, ax = plt.subplots(2, 2, figsize=(10, 10))   
        ax[0][0].imshow(test.original_pacbed, cmap="jet")
        ax[0][0].imshow(test.BF_detector, cmap="gray", alpha=0.5)
        ax[0][0].set_title("BF detector")
        ax[0][0].axis("off")
        ax[0][1].imshow(test.BF_stem, cmap="afmhot")
        ax[0][1].set_title("BF-STEM image")
        ax[0][1].axis("off")
        ax[1][0].imshow(test.original_pacbed, cmap="jet")
        ax[1][0].imshow(test.ADF_detector, cmap="gray", alpha=0.5)
        ax[1][0].set_title("ADF detector")
        ax[1][0].axis("off")
        ax[1][1].imshow(test.ADF_stem, cmap="afmhot")
        ax[1][1].set_title("ADF-STEM image")
        ax[1][1].axis("off")
        fig.tight_layout()

    if "3" in img_num:
        fig, ax = plt.subplots(2, 2, figsize=(10, 10))
        ax[0][0].imshow(test.ADF_stem, cmap="afmhot")
        ax[0][0].set_title("ADF-STEM image")
        ax[0][0].axis("off")
        ax[0][1].imshow(test.E_field_y, cmap="gray")
        ax[0][1].set_title("Electric field y_component")
        ax[0][1].axis("off")
        ax[1][0].imshow(test.E_field_x, cmap="gray")
        ax[1][0].set_title("Electric field x_component")
        ax[1][0].axis("off")
        ax[1][1].imshow(test.E_mag, cmap="inferno")
        ax[1][1].set_title("Electric field magnitude")
        ax[1][1].axis("off")
        fig.tight_layout()

    if "4" in img_num:
        RY, RX = np.indices(test.c_shape[:2])
        fig, ax = plt.subplots(2, 2, figsize=(10, 10))
        ax[0][0].imshow(test.ADF_stem, cmap="gray")
        ax[0][0].quiver(RX.flatten(), RY.flatten(), test.E_field_x.flatten(), test.E_field_y.flatten(), color=cm.jet(mcolors.Normalize()(test.E_mag.flatten())))
        ax[0][0].set_title("Electric field direction")
        ax[0][0].axis("off")
        ax[0][1].axis("off")
        ax[1][0].imshow(test.ADF_stem, cmap="gray")
        ax[1][0].imshow(test.charge_density, cmap="RdBu_r")
        ax[1][0].set_title("Charge density")
        ax[1][0].axis("off")
        ax[1][1].imshow(test.ADF_stem, cmap="gray")
        ax[1][1].imshow(test.potential, cmap="RdBu_r")
        ax[1][1].set_title("Potential")
        ax[1][1].axis("off")
        fig.tight_layout()
    plt.show()
