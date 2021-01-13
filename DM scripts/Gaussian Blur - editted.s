// editting Gaussian Blur DM Script(by D. R. G. Mitchell, October. 2016)

image GaussianConvolution(image sourceimg, number standarddev)
	{
	
		// Trap for zero (or negative) standard deviations
		
		if(standarddev<=0) return sourceimg
	
	
		// get the size of the source image. If it is not a power of 2 in dimension

		number xsize, ysize
		getsize(sourceimg, xsize, ysize)
	
		// Create the gaussian kernel

		image kernelimg:=realimage("",4,xsize, ysize)
		number xmidpoint=xsize/2
		number ymidpoint=ysize/2
		kernelimg=1/(2*pi()*standarddev**2)*exp(-1*(((icol-xmidpoint)**2+(irow-ymidpoint)**2)/(2*standarddev**2)))
		//showimage(kernelimg)

		// Carry out the convolution in Fourier space

		compleximage fftkernelimg:=realFFT(kernelimg)
		//showimage(fftkernelimg)
		compleximage FFTSource:=realfft(sourceimg)
		//showimage(fftsource)
		compleximage FFTProduct:=FFTSource*fftkernelimg.modulus().sqrt()
		//showimage(FFTProduct)
		realimage invFFT:=realIFFT(FFTProduct)
		return invFFT

	}

// Main program to calculate and display a gaussian blurred image

// Check an image is displayed

number nodocs=countdocumentwindowsoftype(5)
if(nodocs<1)
	{
		showalert("Please ensure an image is displayed.",2)
		exit(0)
	}

// Source the image

number standarddev
image front:=getfrontimage()
string imgname=getname(front)


// Prompt for the Standard Deviation of the blur to be used and call the above function

if(!getnumber("Select the standard deviation of the Gaussian Kernel:",3,standarddev)) exit(0)
image gaussblur:=GaussianConvolution(front, standarddev)


// Display the image and copy across the calibrations and tag groups
showimage(gaussblur)
showimage(front-gaussblur)
setname(gaussblur, imgname+" - Gaussian Blur - sigma ("+standarddev+")")

imagecopycalibrationfrom(gaussblur, front)

taggroup fronttags=front.imagegettaggroup()
taggroup gausstags=gaussblur.imagegettaggroup()
taggroupcopytagsfrom(gausstags, fronttags)