// 2019. 11. 4
// JINSEOK RYU, Electron Microscopy and Spectroscopy Lab. Seoul National University
// modifying "Assign ROI as EELS Background" decribed in detail below

// Script to fit a background to an EELS spectrum and extract the edge.
 // This emulates the function of the EELS/Assign ROI As/Background function
 // The script uses two event listeners to listen for an ROI change and an ROI removal
 // If the ROI is changed the background fit is adjusted accordingly. If the ROI
 // is removed, then the corresponding Edge and Background fits are also removed.
 // Only the background ROI is listened for - its ID is sourced and stored in the
 // listener as SpecialROIID. Any other ROIs added are ignored.
 
 // D. R. G. Mitchell, adminnospam@dmscripting.com (remove the nospam to make this address work)
 // January 2015, version:150110, www.dmscripting.com
 
 // Thanks to Steffen Meyer for assistance with the roi_removed function and also Vincent Hou for a 
 // really useful demonstration code on DMSUG, showing how to use Event Maps.
 
 
 // Global variables - the tokens which identify the listeners

number token1, token2, token3, token4

// This class object contains the functions to set and monitor the ROI. Two listeners use this class
// one listening for ROI changes and the other ROI removal. 

class ImageDisplayEventListener : object
{
	// Some variables
		
	ROI 	theroi
	image SI, front, normalized_line
	number 	left, right, SpecialROIID
	image spectrum		
			
	// Linear regression function expecting the data in two images - xdata for the xdata
	// and ydata for the y data. The coefficients (y=a+b*x) a and b are passed in by reference.d

	void linearregression(object self, image xdata, image ydata, number &a, number &b)
		{
			// web source: efunda - the least squares line

			// There is no error trapping for incorrectly formatted data (for speed). Carry out
			// any checks before calling this function
		
			number xsize, ysize
			getsize(xdata, xsize, ysize)

			image xsqrd=xdata**2
			image xtimesy=xdata*ydata

			number sumofy, sumofx, sumofxtimesy, sumofxsqrd

			sumofy=sum(ydata)
			sumofx=sum(xdata)
			sumofxtimesy=sum(xtimesy)
			sumofxsqrd=sum(xsqrd)

			number denominator=xsize*sumofxsqrd-sumofx**2

			a=((sumofy*sumofxsqrd)-(sumofx*sumofxtimesy))/denominator
			b=((xsize*sumofxtimesy)-(sumofx*sumofy))/denominator
		}


	// Following from Egerton's EELS book, the fitting routine involves converting both the 
	// intensity (y) data and the energy (x) data into log. This log-log plot is very close to a straight
	// line and linear regression may then be used to find the a and b coefficients to compute AE^-r
	// where A = the a coefficient and r=-b.

	// The entire spectrum is passed in as spectrum, along with the bounds of the ROI defining the
	// region for background fitting. The background fit parameters A and r are passed in by reference
	// This function calls the above linear regression function, which must be present and defined ahead
	// this function

	image fiteelsbackground(object self, image spectrum, number roistart, number roiend, number &a, number &b)
		{
		
			// Requires the LinearRegression() function above
		
			// getinformation on the passed in spectrum

			number origin, scale, fitrange
			fitrange=roiend-roistart
			string units
			spectrum.imagegetdimensioncalibration(0,origin, scale, units, 0)


			// Set up the log data of the intensity on the region defined by the ROI

			image fitregionintensity=spectrum[0,roistart, 1, roiend]
			ImageSetDimensionCalibration(fitregionintensity,0,origin-roistart,scale, "eV",1)
			fitregionintensity=log(fitregionintensity)


			// Set up the log data of the energy scale on the region defined by the ROI

			image fitregioninx=fitregionintensity*0
			fitregioninx=(icol-origin+roistart)*scale


			// Call the above function to fit a straight line to the log-log plot
			// since the coefficients are passed in by reference into the FitEEELSBackground function
			// the parameters (as well as the background fit) can be obtained from the calling function

			self.linearregression(fitregioninx, fitregionintensity, a, b)


			// On the basis of the above fit, compute the backgroud fit

			image background=spectrum*0
			background=a+b*(icol-origin)*scale // A.E**-r
			background=exp(background)
			return background
		}


	// ROIChanged responds to shifts in the ROI on the image/spectrum
	// The event listener ROIChangeListener calls this function when a change in the ROI is detected
	// As this is called by an event listener, various flags are returned. These are not used here, but they must be defined
	// in the function. In principal it should be possible to check for both removal and change with a single listener
	// and work out what the change is from the returned flags. However, in this case it was simpler to create two listeners
	// on listening for change the other listening for removal. Each listener calls its own function.

	void ROIChanged( Object self, Number event_flags, ImageDisplay imgdisp, Number roi_change_flags, Number roi_disp_change_flags, ROI theroi )
		{		
			// Check to make sure only the original background ROI is responded to. Changes in other ROIs will be ignored
				
			
			imgdisp = spectrum.imagegetimagedisplay(0)
			theroi = imgdisp.imagedisplaygetroi(0)
			number left, right
			theroi.roigetrange(left, right)
			
			if(right > 1)
			{
				
			
			// Calculate the fit and compute the background and edge spectra
				
			number A, r
			image background=self.fiteelsbackground(spectrum, left, right, A, r)
			imgdisp.lineplotimagedisplaysetlegendshown(1)
			
			background[0,0,1,left]=0
			image edge=spectrum-background
			edge[0,0,1,left]=0
			
			normalized_line = edge / max(edge)

				
			// remove any existing edge or background slices
				
			number noslices=imgdisp.imagedisplaycountslices()
			number i
			for(i=noslices-1; i>-1; i--)
				{
					object sliceid=imgdisp.imagedisplaygetsliceidbyindex(i)
					string slicename=imgdisp.imagedisplaygetslicelabelbyid(sliceid)
					if(slicename=="Bkgd" || slicename=="Edge") imgdisp.imagedisplaydeleteslicewithid(sliceid)
				}
				
			imgdisp.imagedisplayaddimage(background, "Bkgd")
			imgdisp.imagedisplayaddimage(edge, "Edge")			
			}
		}


	// Responds when the ROI is removed - see comments above on ROIChanged function
		
	void ROIRemoved( Object self, Number event_flags, ImageDisplay imgdisp, Number roi_change_flags, Number roi_disp_change_flags, ROI theroi )
		{		
			// Remove the listeners and remove the background and edge fits
			
			imgdisp.ImageDisplayRemoveEventListener(token1)
			imgdisp.ImageDisplayRemoveEventListener(token2)
			imgdisp.ImageDisplayRemoveEventListener(token3)
			imgdisp.ImageDisplayRemoveEventListener(token4)
			
			number noslices=imgdisp.imagedisplaycountslices()
			number i
			for(i=noslices-1; i>-1; i--)
				{
					object sliceid=imgdisp.imagedisplaygetsliceidbyindex(i)
					string slicename=imgdisp.imagedisplaygetslicelabelbyid(sliceid)
					if(slicename=="Bkgd" || slicename=="Edge") imgdisp.imagedisplaydeleteslicewithid(sliceid)
				}
		}



// This function sources the ROI on the passed in image
	
	object init(object self, image spectrum_image, image front, image normalized)
		{
			//result("\n initiation activated")
			
			
			SI := spectrum_image
			normalized_line := normalized			
			spectrum:=front
			imagedisplay imgdisp=spectrum.ImageGetImageDisplay(0)
			theroi = imgdisp.ImageDisplayGetROI(0)
			
			// Store the ID of the background ROI to avoid confusion with and changes/deletions
			// of other ROIs on the spectrum
			
			SpecialROIID=theroi.roigetid()


			// Set the ROI to the appropriate colour and form
			
			theroi.roisetvolatile(0)
			theroi.roisetcolor(1,0,0)
			theroi.roisetlabel("bkgd")
			
			
			// Call the ROIChanged function to ensure that the fit is applied at first use
			// (when there is no change). Note the ROIChanged function is called by the 
			// Event Listener and has various flags associated with it. These are not used
			// and so dummy values are used to keep the function syntax acceptable
			
			number dummy=0
			self.roichanged(dummy, imgdisp, dummy, dummy, theroi)
			return self
		}
		
	// Constructor
	
	ImageDisplayEventListener(object self)
		{
		result("\n ############################################")
		result("\n event listener object activated successfully")
		result("\n ############################################")
		}
		
	// Destructor
	
	~ImageDisplayEventListener(object self)
		{
		result("\n #############################################")
		result("\n event listener object terminated successfully")
		result("\n #############################################")
		}
}



// function to check the displayed image is appropriate and create and apply
// the event listeners

void main()
	{
		image SI, front, normalized_line
		number sx,sy,sz,xscale,yscale,zscale,xorigin,yorigin,zorigin
		string prompt

		if(!(getoneimagewithprompt(prompt, "Select a 3D image", SI))) exit(0)
		if(!(getoneimagewithprompt(prompt, "Select a spectrum with a ROI", front))) exit(0)

		SI.get3dsize(sx,sy,sz)

		xscale = SI.imagegetdimensionscale(0)
		yscale = SI.imagegetdimensionscale(1)
		zscale = SI.imagegetdimensionscale(2)

		xorigin = SI.imagegetdimensionorigin(0)
		yorigin = SI.imagegetdimensionorigin(1)
		zorigin = SI.imagegetdimensionorigin(2)

		normalized_line := realimage("dummy", 4, sz, 1)
		normalized_line.setname("normalized line")
			
		showimage(normalized_line)
		
		normalized_line.ImageSetDimensionOrigin( 0, zorigin ) 
		normalized_line.ImageSetDimensionScale( 0, zscale ) 
		normalized_line.ImageSetDimensionUnitString(0, "eV" )
		
		number nodocs=countdocumentwindowsoftype(5)
		if(nodocs<1)
			{
				showalert("Ensure a spectrum with a Region of Interest on it is shown.",2)
				exit(0)
			}


		// source the front-most image and check that it is a 1D profile
	
		imagedisplay SIdisp = SI.ImagegetImagedisplay(0)
		number xsize, ysize
		getsize(front, xsize, ysize)
		if(ysize>1)
			{
				showalert("This only works on EELS Spectra.",2)
				exit(0)
			}


		// Check to make sure there is an ROI on the selected image

		imagedisplay imgdisp=front.imagegetimagedisplay(0)
		number norois=imgdisp.imagedisplaycountrois()
		if(norois<1)
			{
				 showalert("Please ensure a ROI is present on the front-most image.",2)
				 exit(0)
			 }


		// Allocate the listener objects to the front image

		// Listener for ROI removal

		string messagemap1="roi_removed:ROIRemoved"
		object ROIRemovalListener=alloc(ImageDisplayEventListener).init(SI, front, normalized_line)


		// Listener for ROI change

		string messagemap2="roi_changed:ROIChanged"
		object ROIChangeListener=alloc(ImageDisplayEventListener).init(SI, front, normalized_line)
		

		// Add the event listeners

		token1 = imgdisp.ImageDisplayAddEventListener( RoiRemovalListener, messagemap1)
		token2 = imgdisp.ImageDisplayAddEventListener( RoiChangeListener, messagemap2)
		token3 = SIdisp.ImageDisplayAddEventListener( RoiRemovalListener, messagemap1)
		token4 = SIdisp.ImageDisplayAddEventListener( RoiChangeListener, messagemap2)
	}


// Main script

main()
