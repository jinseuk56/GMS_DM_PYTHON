# Ingyu Yoo, Jinseok Ryu
# Electron Microscopy and Spectroscopy Lab.
# Seoul National University
# last update : 20210906
# differential phase contrast imaging for 4D-STEM data
# You can adjust the parameters (pass filter) when acquiring the iDPC image.

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
    
    def find_center(self):
        Y, X = np.indices(self.original_pacbed.shape)
        com_y = np.sum(self.original_pacbed * Y) / np.sum(self.original_pacbed)
        com_x = np.sum(self.original_pacbed * X) / np.sum(self.original_pacbed)
        
        self.ct = [com_y, com_x]
        print("the center of the pacbed (y, x) =", self.ct)
        
    def disk_extract(self):
        grad = np.gradient(self.original_pacbed)
        grad_map = grad[0]**2 + grad[1] **2
        grad_map = grad_map / np.max(grad_map)
        
        max_ind = np.unravel_index(np.argmax(grad_map, axis=None), grad_map.shape)
        self.least_R = ((max_ind[0]-self.ct[0])**2 + (max_ind[1]-self.ct[1])**2)**(1/2)
        
        print("radius of the BF disk = %.2f mrad"%(self.mrad_per_pixel*self.least_R))
        
        self.ct_ind  = np.around(self.ct).astype(int)
        self.cropped_size = np.min(self.ct_ind)
        print("radius of the RoI = %.2f mrad"%(self.mrad_per_pixel*self.cropped_size))
        
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
        
    def DPC(self, correct_rotation=True, n_theta=100, hpass=0.1, lpass=0.1):
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

    FT_phase = (FT_ysh*ky + FT_xsh*kx) / (2*np.pi*1j*(hpass+k2+lpass*k2))
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

test = dpc_python(f_stack, scale0, scale2)
test.find_center()
test.disk_extract()

BF_det = [0, int(test.least_R)]
ADF_det = [int(test.least_R), int(np.min(test.ct))]
test.virtual_stem(BF_det, ADF_det)

test.DPC()
print("Optimized angle = {} degree".format(test.c_theta*180/np.pi))

im1 = DM.CreateImage(test.c_pacbed.copy())
im1.SetName("PACBED")
im1.ShowImage()

im2_1 = DM.CreateImage(test.BF_stem.copy())
im2_1.SetName("BF-STEM image")
im2_1.ShowImage()
im2_2 = DM.CreateImage(test.ADF_stem.copy())
im2_2.SetName("ADF-STEM image")
im2_2.ShowImage()

im3_1 = DM.CreateImage(test.E_field_y.copy())
im3_1.SetName("Electric field y_component")
im3_1.ShowImage()
im3_2 = DM.CreateImage(test.E_field_x.copy())
im3_2.SetName("Electric field x_component")
im3_2.ShowImage()
im3_3 = DM.CreateImage(test.E_mag.copy())
im3_3.SetName("Electric field magnitude")
im3_3.ShowImage()
      
im4_1 = DM.CreateImage(test.charge_density.copy())
im4_1.SetName("Charge density")
im4_1.ShowImage()
im4_2 = DM.CreateImage(test.potential.copy())
im4_2.SetName("Potential")
im4_2.ShowImage()


py_figure = input("Do you also want the Python figures of the results? (Y/N)")
if py_figure == "Y":
    fig1, ax1 = plt.subplots(1, 1, figsize=(5, 5))
    ax1.imshow(test.c_pacbed, cmap="jet")
    ax1.set_title("PACBED")
    ax1.scatter(test.c_ct[1], test.c_ct[0], s=15, c="k")
    ax1.axis("off")
    fig1.tight_layout()

    fig2, ax2 = plt.subplots(2, 2, figsize=(10, 10))   
    ax2[0][0].imshow(test.original_pacbed, cmap="jet")
    ax2[0][0].imshow(test.BF_detector, cmap="gray", alpha=0.5)
    ax2[0][0].set_title("BF detector")
    ax2[0][0].axis("off")
    ax2[0][1].imshow(test.BF_stem, cmap="gray")
    ax2[0][1].set_title("BF-STEM image")
    ax2[0][1].axis("off")
    ax2[1][0].imshow(test.original_pacbed, cmap="jet")
    ax2[1][0].imshow(test.ADF_detector, cmap="gray", alpha=0.5)
    ax2[1][0].set_title("ADF detector")
    ax2[1][0].axis("off")
    ax2[1][1].imshow(test.ADF_stem, cmap="gray")
    ax2[1][1].set_title("ADF-STEM image")
    ax2[1][1].axis("off")
    fig2.tight_layout()

    fig3, ax3 = plt.subplots(2, 2, figsize=(10, 10))
    ax3[0][0].imshow(test.BF_stem, cmap="afmhot")
    ax3[0][0].set_title("BF-STEM image")
    ax3[0][0].axis("off")
    ax3[0][1].imshow(test.E_field_y, cmap="gray")
    ax3[0][1].set_title("Electric field y_component")
    ax3[0][1].axis("off")
    ax3[1][0].imshow(test.E_field_x, cmap="gray")
    ax3[1][0].set_title("Electric field x_component")
    ax3[1][0].axis("off")
    ax3[1][1].imshow(test.E_mag, cmap="inferno")
    ax3[1][1].set_title("Electric field magnitude")
    ax3[1][1].axis("off")
    fig3.tight_layout()

    RY, RX = np.indices(test.c_shape[:2])
    fig4, ax4 = plt.subplots(2, 2, figsize=(10, 10))
    ax4[0][0].imshow(test.ADF_stem, cmap="gray")
    ax4[0][0].quiver(RX.flatten(), RY.flatten(), test.E_field_x.flatten(), test.E_field_y.flatten(), color=cm.jet(mcolors.Normalize()(test.E_mag.flatten())))
    ax4[0][0].set_title("Electric field direction")
    ax4[0][0].axis("off")
    ax4[0][1].axis("off")
    ax4[1][0].imshow(test.ADF_stem, cmap="gray")
    ax4[1][0].imshow(test.charge_density, cmap="RdBu_r")
    ax4[1][0].set_title("Charge density")
    ax4[1][0].axis("off")
    ax4[1][1].imshow(test.ADF_stem, cmap="gray")
    ax4[1][1].imshow(test.potential, cmap="RdBu_r")
    ax4[1][1].set_title("Potential")
    ax4[1][1].axis("off")
    fig4.tight_layout()
    plt.show()