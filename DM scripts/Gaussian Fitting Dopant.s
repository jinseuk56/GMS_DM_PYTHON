// 2D Gaussian Fitting Part


// 2D Gaussian fitting functions which use a Monte Carlo approach. Use this method to find the centre
// of a blob image - such as a diffraction spot or atomic column.

 
// D. R. G. Mitchell, adminnospam@dmscripting.com (remove the nospam)
// version:20150720, v1.0, July 2015.


// Function to draw a 2D Gaussian

image Create2DGaussian(number xsize, number ysize, number centrex, number centrey, number intensity, number sigmax, number sigmay)
	{
		image gaussian=realimage("", 4,xsize, ysize)
		if(centrex>xsize || centrex<0 || centrey>ysize || centrey<0) return gaussian
		gaussian=intensity*exp(-(icol-centrex)**2/(2*sigmax**2))*exp(-(irow-centrey)**2/(2*sigmay**2))
		return gaussian
	}


// Function to compute chi squared by fitting an original image to a proposed fit image.
// The centralpcnt is the area (as a percentage of the original size) over which the 
// comparison is made (25-75%). This region is centred on the estimated centre of the passed in image
// In this case, fit is the passed in Gaussian image.

number computechisqrd(image original, image fit, number centrex, number centrey, number centralpcnt)
	{
		number xsize, ysize
		getsize(original, xsize, ysize)
		
		
		/*// Compare only subarea of the image - centred on the target centre.
		
		if(centralpcnt<25) centralpcnt=25
		if(centralpcnt>75) centralpcnt=75
		*/
		
		
		// Keep the size of the central region sensible >=10 pixels
		
		number xdim=xsize
		//if(xdim<10) xdim=10
		
		number ydim=ysize
		//if(ydim<10) ydim=10
		
		
		// Define the search area centred on the target centre
		

		/*
		if(top<0) top=0
		if(left<0) left=0
		if(bottom>ysize-1) bottom=ysize-1
		if(right>xsize-1) right=xsize-1
		*/
		

		// measure the difference between the original and the fit (gaussian) images within their subareal regions
		// and compute chisqrd
		
		image difference=(original-fit)**2/(original+0.000001)
		number chisqrd=sum(difference)
		return chisqrd
	}


// Uses a random walk to find the optimum sigma value which fits a 2D Gaussian to a blob image (front)
// Centre x and y are the estimated centres of the blob and intensity its intensity - normally normalised to 1. Initialsigma
// is the starting value of sigma (usually estimated as 1/4 of the  x dimension of the image. Iterations are the number of interations
// for the walk. Subareapcnt is the percentage of the sub-area part of the image within which to hunt for the match. Sigma is the final
// sigma value. Minchisqrd is the current chi sqrd value of prior fitting - eg a centre optimisation. This enables this function
// to only apply changes which improve the fit on previous fits. The function returns the chisqrd of the optimum fit.

// This function has dependencies on Create2DGaussian() and ComputeChiSqrd() functions

number FindSigmaMonteCarlo(image front, number centrex, number centrey, number subareapcnt, number intensity, number initialsigma, number iterations, number &minchisqrd, number &sigma)
	{
		// Source the image size

		number xsize, ysize, chisqrd
		getsize(front, xsize, ysize)
		number i
		

		// random walk - vary the initial sigma with a random gaussian function and compute the fit for the resulting sigma
		
		for(i=0; i<iterations; i++)
			{
				number sigmatest=initialsigma+gaussianrandom()
				if(sigmatest<1) sigmatest=1
				image testgaussian=Create2DGaussian(xsize, ysize, centrex, centrey, intensity, sigmatest, sigmatest)			
				chisqrd=computechisqrd(front, testgaussian, centrex, centrey, subareapcnt)

				if(chisqrd<minchisqrd) 
					{
						minchisqrd=chisqrd
						sigma=sigmatest
						initialsigma=sigmatest
					}	

			}
		return minchisqrd
	}


// Monte Carlo function to vary the centre of a Gaussian function and to choose those variants which improve the fit to the
// experimental image (front). Sigma is the sigma of the model Gaussian function. Iterations are the number of iterations for
// the random walk. Intensity in the Gaussian intensity - usually 1, centrex and y are the initial centre coordinates. Subareapcnt 
// percentage of the image (as a subarea centred on the guessed centre) which will be compared for the fit. Minchisqrd is the 
// minimum chi squared value obtained from previous fitting (eg Sigma optimisation). This ensures that this function will 
// only apply changes which improve on preceding fits. Fitcentrex and y are the fitted coordinates of an improved fit.

// This function has dependencies on Create2DGaussian() and ComputeChiSqrd() functions


number FindCentreMonteCarlo(image front, number sigma, number intensity, number iterations, number centrex, number centrey, number subareapcnt, number &minchisqrd, number &fitcentrex, number &fitcentrey)
	{
		number xsize, ysize, chisqrd, newx, newy
		getsize(front, xsize, ysize)
		image testgaussian=front*0
	
	
		// If no centre information is provided  estimate the centre as the geometric centre
		
		number i, randcentrex, randcentrey
		fitcentrex=centrex
		fitcentrey=centrey
	
	
		// Do a random walk looking for an improved fit
		
		for(i=0; i<iterations; i++)
			{
				// vary the centres with a random walk
				
				randcentrex=gaussianrandom()
				randcentrey=gaussianrandom()
				newx=randcentrex+fitcentrex
				newy=randcentrey+fitcentrey
				
				
				// Ignore any values which take the centre outside the bounds of the image
				
				if(newx<0) randcentrex=0
				if(newx>xsize-1) randcentrex=0
				if(newy<0) randcentrey=0
				if(newy>ysize-1) randcentrey=0
			
			
				// test the fit with the new centre
				
				testgaussian=Create2DGaussian(xsize, ysize, newx,newy, intensity, sigma, sigma)
				chisqrd=computechisqrd(front, testgaussian, newx, newy, subareapcnt)

				if(chisqrd<minchisqrd) 
					{
						minchisqrd=chisqrd
						fitcentrex=newx
						fitcentrey=newy
					}	
			}
		
		return minchisqrd
	}

//Gaussian Fitting Function

image GaussianFitting(image ClusterRawData, image tmp_target, number sx, number sy, number ys, number scale)
{
number NoFittingAtomCount, FittingAtomCount
number top, left, bottom, right
number atomcnt1, atomcnt2
NofittingAtomcount = 0
FittingAtomCount = 0
for(atomcnt1=0;atomcnt1<ys;atomcnt1++)
{
	if(getpixel(ClusterRawData, 2, atomcnt1)!=0)
	{
		if(getpixel(ClusterRawData,1,atomcnt1)/scale<3 || getpixel(ClusterRawData,0,atomcnt1)/scale<3 || getpixel(ClusterRawData,1,atomcnt1)/scale>sy-4 || getpixel(ClusterRawData,0,atomcnt1)/scale>sx-4)
		{
		//result("no fitting "+getpixel(ClusterRawData,0,atomcnt1)/scale+" "+getpixel(ClusterRawData,1,atomcnt1)/scale+"\n")
		NoFittingAtomCount += 1
		}
		else
		{
		//result("fitting "+getpixel(ClusterRawData,0,atomcnt1)/scale+" "+getpixel(ClusterRawData,1,atomcnt1)/scale+"\n")
		FittingAtomCount += 1
		}
		//result(top+", "+left+", "+bottom+", "+right+"\n")
	}
}
//result(NoFittingAtomCount+", "+FittingAtomCount+"\n")

image NoFittingAtom := realimage("", 4, 2, NoFittingAtomCount)
image FittingAtom := realimage("", 4, 2, FittingAtomCount)
number count1 = 0
number count2 = 0

for(atomcnt2=0;atomcnt2<ys;atomcnt2++)
{

	if(getpixel(ClusterRawData, 2, atomcnt2)!=0)
	{
		
		if(getpixel(ClusterRawData,1,atomcnt2)/scale<3 || getpixel(ClusterRawData,0,atomcnt2)/scale<3 || getpixel(ClusterRawData,1,atomcnt2)/scale>sy-4 || getpixel(ClusterRawData,0,atomcnt2)/scale>sx-4)
		{
		setpixel(NoFittingAtom, 0, count1, getpixel(ClusterRawData,0,atomcnt2)/scale)
		setpixel(NoFittingAtom, 1, count1, getpixel(ClusterRawData,1,atomcnt2)/scale)
		count1 += 1
		}

//showimage(NoFittingAtom)
//NoFittingAtom.setdisplaytype(5)

		else
		{
		top = getpixel(ClusterRawData, 1, atomcnt2)/scale-3
		bottom = getpixel(ClusterRawData, 1, atomcnt2)/scale+4
		left = getpixel(ClusterRawData, 0, atomcnt2)/scale-3
		right = getpixel(ClusterRawData, 0, atomcnt2)/scale+4
		//result(top+", "+left+", "+bottom+", "+right+"\n")

		
		// Gaussian Fitting of the constructed ROIs
		
		image front = tmp_target[top, left, bottom, right]
		number xsize, ysize
		getsize(front, xsize, ysize)


		//Normalise the image so intensities run from 0 to 1

		number minval, maxval
		minmax(front, minval, maxval)
		number range=maxval-minval
		front=(front-minval)/range


		// Set up some initial values

		number centrex = round(xsize/2) 
		number centrey = round(ysize/2)
		
		// or to produce a better result

		number subareapcnt=100 // different from the original script

		number fitcentrex, fitcentrey
		number i, chisqrd
		number intensity=1 // set to 1 for normalised images
		number initialsigma
		number iterations=50 // The number of random steps taken to refine the value


		// Note the processing time (for a ca 50 x 50 blob image) is a function of iterations
		// multiplied by refinementloops (below). Values of 50 iterations x 5 refinements takes about 0.4s.

		number sigma=xsize/4 // Sigma of the Gaussian - estimate an initial value before refinement
		initialsigma=sigma
		number minchisqrd=1e99



		// Loop to refine sigma, then Centre then Sigma etc

		number refinementloops=10
		number firstchi, secondchi
		number counter
			for(i=0; i<refinementloops; i++)
			{
			// Refine sigma

			minchisqrd=FindSigmaMonteCarlo(front, centrex, centrey, subareapcnt, intensity, initialsigma, iterations, minchisqrd, sigma)
			initialsigma=sigma		
			//result("\n\n\nSigma refinement  : "+(i+1)+" Sigma = "+format(sigma, "%3.3f")+" Chi Squared = "+format(minchisqrd, "%3.3f"))
			firstchi=minchisqrd
		
		
			// Refine the centre

			image gaussian=Create2DGaussian(xsize, ysize, centrex, centrey, intensity, sigma,sigma)
			minchisqrd=FindCentreMonteCarlo(front, sigma, intensity, iterations, centrex, centrey, subareapcnt, minchisqrd, fitcentrex, fitcentrey)
			centrex=fitcentrex
			centrey=fitcentrey
			//result("\nCentre refinement : "+(i+1)+" Centre x = "+format(fitcentrex, "%3.3f")+"  Centre y = "+format(fitcentrey, "%3.3f")+"  Chi Squared = "+format(minchisqrd, "%3.3f"))

		
			// If the chisqrd value for both sigma and centre optimisation are the same, increment a counter.
			// If that counter reaches 3, then assume the iteration is optimised and curtail further iteration.
				
			secondchi=minchisqrd
				if(firstchi==secondchi) counter=counter+1
				else counter=0
		
				if(counter==3) i=refinementloops
			}
		setpixel(FittingAtom, 0, count2, round(centrex)+left)
		setpixel(FittingAtom, 1, count2, round(centrey)+top)
		count2 += 1
		}
	}	
}

//showimage(NoFittingAtom)
//NoFittingAtom.setdisplaytype(5)
//showimage(FittingAtom)
//FittingAtom.setdisplaytype(5)

image RealTargetPosition := realimage("Real Target Imgage", 4, 2, NoFittingAtomCount+FittingAtomCount)
number atomcnt3
for(atomcnt3=0;atomcnt3<NoFittingAtomCount+FittingAtomCount;atomcnt3++)
{
	if(atomcnt3<NoFittingAtomCount)
	{
	setpixel(RealTargetPosition, 0, atomcnt3, getpixel(NoFittingAtom, 0, atomcnt3))
	setpixel(RealTargetPosition, 1, atomcnt3, getpixel(NoFittingAtom, 1, atomcnt3))
	}
	else
	{
	setpixel(RealTargetPosition, 0, atomcnt3, getpixel(FittingAtom, 0, atomcnt3-NoFittingAtomCount))
	setpixel(RealTargetPosition, 1, atomcnt3, getpixel(FittingAtom, 1, atomcnt3-NoFittingAtomCount))
	}
}
showimage(RealTargetPosition)
RealTargetPosition.setdisplaytype(5)

image RealTargetImg := realimage("RealTargetImg", 4, sx, sy)
number atomcnt4
for(atomcnt4=0;atomcnt4<NoFittingAtomCount+FittingAtomCount;atomcnt4++)
{
setpixel(RealTargetImg, getpixel(RealTargetPosition, 0, atomcnt4), getpixel(RealTargetPosition, 1, atomcnt4), 100)
}
return RealTargetImg
}

	
image ClusterRawData := B
number xs, ys
image TargetImg := C
number sx, sy
number scale = 0.0233468
getsize(ClusterRawData, xs, ys)
getsize(TargetImg, sx, sy)
image tmp_target = TargetImg
image RealTargetImg

RealTargetImg = GaussianFitting(ClusterRawData, tmp_target, sx, sy, ys, scale)
showimage(RealTargetImg)



