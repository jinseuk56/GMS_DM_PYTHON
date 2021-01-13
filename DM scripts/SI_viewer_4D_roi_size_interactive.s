// 2019. 11. 5
// JINSEOK RYU, Electron Microscopy and Spectroscopy Lab., Seoul National University
// SI viewer for 4D data (outer 2D + inner 2D)

// 2019. 11. 11
// it is also interactive to the roi size, but it gets slower when the size of the roi is larger

number el1, el2

class imagedisplayeventlistener : object
{
	ROI theroi
	image fourd, twod
	number sx, sy, isx, isy
	number rt, rl, rb, rr
	number specialROIID

	void ROIchanged(object self, number ef, imagedisplay imgdisp, number rcf, number rdcf, ROI theroi)
	{
	
		number thisroiid=theroi.roigetid()
		if(thisroiid!=SpecialROIID) return

		theroi.roigetrectangle(rt, rl, rb, rr)
		//result("\n ROI changed")
		//result("\n"+rt+" "+rl+" "+rb+" "+rr)
	
		number i, j
		twod = 0
		for(i=rl;i<rr;i++){
			for(j=rt;j<rb;j++){
				twod += fourd.sliceN(4, 2, i, j, 0, 0, 2, isx, 1, 3, isy, 1)
			}
		}
	}
	
	void ROIremoved(object self, number ef, imagedisplay imgdisp, number rcf, number rdcf, ROI theroi)
	{
		number thisroiid=theroi.roigetid()
		if(thisroiid!=SpecialROIID) return
		
		imgdisp.imagedisplayremoveeventlistener(el1)
		imgdisp.imagedisplayremoveeventlistener(el2)	
	}
	
	object init(object self, image II, image EI)
	{
	
		fourd := II
		twod := EI
		
		sx = fourd.ImageGetDimensionSize(0)
		sy = fourd.ImageGetDimensionSize(1)
		isx = fourd.ImageGetDimensionSize(2)
		isy = fourd.ImageGetDimensionSize(3)
	
		imagedisplay imgdisp = II.imagegetimagedisplay(0)
		theroi = imgdisp.imagedisplaygetroi(0)
		
		specialROIID = theroi.roigetid()
		
		theroi.roisetvolatile(0)
		theroi.roisetcolor(1,0,0)
		theroi.roisetlabel("selected area")
		
		number dummy = 0
		self.ROIchanged(dummy, imgdisp, dummy, dummy, theroi)
		return self
	}
	
	
	ImageDisplayEventListener(object self)
		{
		result("\n ############################################")
		result("\n event listener object activated successfully")
		result("\n ############################################")
		}
	
	~ImageDisplayEventListener(object self)
		{
		result("\n #############################################")
		result("\n event listener object terminated successfully")
		result("\n #############################################")
		}

}


void main()
{

	image fourd, twod, temp
	number sx, sy, isx, isy
	
	fourd.getfrontimage()
	if(fourd.imagegetnumdimensions() != 4){
		showalert("the front-most image must a 4D image", 2)
		exit(0)
	}
	
	imagedisplay imgdisp = fourd.imagegetimagedisplay(0)
	number noroi = imgdisp.imagedisplaycountrois()
	if(noroi < 1){
		showalert("please ensure a ROI is present on the front-most image", 2)
		exit(0)
	}

	sx = fourd.ImageGetDimensionSize(0)
	sy = fourd.ImageGetDimensionSize(1)
	isx = fourd.ImageGetDimensionSize(2)
	isy = fourd.ImageGetDimensionSize(3)
	
	twod := realimage("2D in 4D", 4, isx, isy)
	twod = 0
	
	showimage(twod)
	
	
	string msmap1 = "roi_changed:ROIchanged"
	object ROIchange = alloc(imagedisplayeventlistener).init(fourd, twod)
	
	string msmap2 = "roi_removed:ROIremoved"
	object ROIremove = alloc(imagedisplayeventlistener).init(fourd, twod)
	
	el1 = imgdisp.ImageDisplayAddEventListener(ROIchange, msmap1)
	el2 = imgdisp.ImageDisplayAddEventListener(ROIremove, msmap2)
}

main()