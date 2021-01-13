/*  Jusang Lee Edited to import new Velox .emd files (jusang0154@snu.ac.kr), see line 423,429 index*3
*	
*	This script import FEI's new (Ver. >1.2) Velox .emd files.
*	The metadata are dumped in image as hierachical tags
*	and calibrations are inserted in the unit of nm.
*
*	Ver. 0.1	Last update: 2016-10-06
*	Dr. Mingjian Wu 	FAU Erlangen-Nuremberg
*	mingjian.wu@fau.de
*
*   This script requires HDF5 library implemented by Tore Niermann
*		https://github.com/niermann/gms_plugin_hdf5
*	the following library should be present in DM plug-in folder:
*		Win64:	hdf5_GMS2X_amd64.dll
*		Win32	hdf5_GMS2X_x86.dll
*	
*	Limitations:	
*   Spectrum RAW DATA are not yet supported in the current version.
*   It is possible to read the images stored in SI, however.
*
*	This script DOES NOT work with OLD Velox (Ver. 1.1.x) files!
*	Dr. Christoph Gammer (Leobon) had implemented a script for that.
*	Please ask him for the corresponding script.
*/



//============== Parsing Velox metadata to image tags ==============
//	Thanks to Bernhard Schaffer, adopted from
//		http://stackoverflow.com/questions/39712122/
//			how-to-wrap-nested-string-metadata-into-taggroup
//==================================================================
Class MetaStr2TagGroup
{
    // Find next string bracketed by " in input string, starting search
    // at given index. Returns string and end-position of search
    string FindNextKeyName( object self, string input, number & searchPos )
    {
        number totalLength = len( input )
        number start = 0, end = 0
        while( searchPos < totalLength )
        {
            searchpos++
            if ( "\"" == input.mid(searchpos-1,1) )
            {
                if ( !start ) 
                    start = searchpos-1
                else
                    {
                    end = searchpos-1
                    return input.mid(start+1,end-start-1)
                    }
            }
        }
        return ""
    }

    // Returns the next of either "{" ,  "}" or """ after a collon ":"
    string GetNextIndicator( object self, string input, number & searchPos )
    {
        number totalLength = len( input )
        while( searchPos < totalLength )
        {
            searchpos++        
            if ( "{" == input.mid(searchpos-1,1) )
                return "{"
            if ( "}" == input.mid(searchpos-1,1) )
                return "}"
            if ( "\"" == input.mid(searchpos-1,1) )
                return "\""
        }
        return ""
    }

    // In a tag-path string, find location of last colon 
    number findLastColon( object self, string input )
    {
        number totalLength = len( input )
        number lastPos = -1
        number searchPos = 0
        while( searchPos < totalLength )
        {
            searchpos++
            if ( ":" == input.mid(searchpos-1,1) )
                lastPos = searchpos-1
        }
        return lastPos
    }

    // Parse text string and write it to image tag
    // Also returns a TagGroup object
    TagGroup ParseText2ImageTag( object self, string input, image &img )
    {
    
    //result("input string:\t"+input+"\n")
    
    //The following part replace all the "]." and "[" by "-". DM doesn't accept them as tag strings
    number IndexDot =0
    string newinput=""
    //result("Input Length\t"+ len(input)+"\n")
    while( IndexDot !=-1)
	{

		indexDot = input.find( "]." )
//		result(indexDot+"\t")
		newinput+=input.left(indexdot)+"-"
		input=input.right(len(input)-indexdot-2)
		
		//input=input.left(indexdot)+"-"+input.right(len(input)-indexdot-2)
//		result(input+"\n")

	}  
	
	input=Newinput
    Newinput=""
    indexdot=0
    while( IndexDot !=-1)
	{
		indexDot = input.find( "[" )
//		result(indexDot+"\t")
		newinput+=input.left(indexdot)+"-"
		input=input.right(len(input)-indexdot-1)
		
		//input=input.left(indexdot)+"-"+input.right(len(input)-indexdot-1)
//		result(input+"\n")

	}  

	input=Newinput
	
    //result(input+"\n")
    
		TagGroup rootTG = NewTagGroup()
		string currentPath = ""
        number totalLength = len( input )
        number searchPos = 0 
        number searchPos2
        string keyName, indicator
        while( searchPos < totalLength )
        {
            // search for new key or closing bracket, whatever first 
            searchPos2 = searchPos
            indicator = self.GetNextIndicator( input, searchPos2 )
            keyName = self.FindNextKeyName( input, searchPos )
            if ( ( "}" == indicator ) && (searchpos2<searchPos ) )
            {
                // decrease hierarchy
                number cutPos = self.findLastColon( currentPath )
                currentPath = left( currentPath, cutPos )
                //result("\n DEC ")
                searchPos = searchPos2
                //Result(searchPos+"\tCurrentPath:\t"+CurrentPath+"\n")
            }
            else	// Either add value or new  sub-tagGroup
            {
				//Result(searchPos+"\tkeyName:\t"+keyname+"\n")
                if ( "" == keyname ) break; 
					// No more keys found
                indicator = self.GetNextIndicator( input, searchPos )
                //Result(searchPos+"\tNextIndicator:\t"+Indicator+"\n")
                if ( "" == indicator ) break;   
					// No more indicator found -- should not happen!

                if ( "{" == indicator ){
					// increase hierachy		
					
                    currentPath += ":" + keyname
                   // Result("\tCurrentPath:\t"+CurrentPath+"\n")				
                    rootTg.TagGroupSetTagAsTagGroup( currentPath, NewTagGroup() )
                } else if ( "\"" == indicator ){
                    // Add value
                    searchPos--
                    string valStr = self.FindNextKeyName( input, searchPos )
                    rootTg.TagGroupSetTagAsString( currentPath + ":" + keyname, valStr )
                    img.SetStringNote(currentPath + ":" + keyname, valStr)
                }
            }
        }
        return rootTg
    }
}


//========================== defining functions ========================

String img2str (image img){
// Velox stores the metadata as ASCII valued image dataset,
// this function convert the image to a single string
	Number sx, sy
	object stream = NewStreamFromBuffer( 0 )
	img.GetSize(sx,sy)
	ImageWriteImageDataToStream( img, stream, 0 )
	stream.StreamSetPos(0,0)
	Return StreamReadAsText( stream, 0, sx*sy )
}

String VeloxEMDGetDataTypeName(Number dtype)
{
    if (dtype == 1)
        return "Int16"
    else if (dtype == 2)
        return "Float32"
    else if (dtype == 3)
        return "Complex64"
    else if (dtype == 6)
        return "UInt8"
    else if (dtype == 7)
        return "Int32"
    else if (dtype == 9)
        return "UInt8"
    else if (dtype == 10)
        return "UInt16"
    else if (dtype == 11)
        return "UInt32"
    else if (dtype == 12)
        return "Float64"
    else if (dtype == 12)
        return "Complex128"
    else if (dtype == 39)
        return "Int64"
    else if (dtype == 30)
        return "UInt64"
    else
        return "<Unknown>"
}

void VeloxEMDContentWalker(TagGroup content, TagGroup &names, TagGroup &infos)
{
    String type, name

    if (!TagGroupGetTagAsString(content, "Type", type) || \
		!TagGroupGetTagAsString(content, "Name", name))
        return;
	//result ("Name: "+name+" \n Type:" +type+"\n\n")
    if ((type == "DataSet")) {
        // If data type is not supported, no "DataType" tag exists, and the query fails
        Number dtype
        if (!TagGroupGetTagAsNumber(content, "DataType", dtype))
            return

        Number rank
        if (!TagGroupGetTagAsNumber(content, "Rank", rank) || rank < 0 || rank > 4)
            return          // DM only supports up to 4D

        // Create size string
        String size_text = ""
        if (rank > 0) {
            TagGroup size
            if (!TagGroupGetTagAsTagGroup(content, "Size", size))
                return

            Number n
            for (n = 0; n < rank; n++) {
                Number dim
                TagGroupGetIndexedTagAsLong(size, n, dim)
                if (n > 0)
                    size_text = size_text + "x" + dim
                else
                    size_text = "" + dim
            }
        }

        // Create info string
        String info_text = name.PathExtractBaseName(0) + " - " + \
			VeloxEMDGetDataTypeName(dtype) + "[" + size_text + "]"
        TagGroupInsertTagAsString(names, -1, name)
        if (name.PathExtractBaseName(0) == "Data")
			TagGroupInsertTagAsString(infos, -1, info_text)
    } else if (type == "Group") {
        TagGroup group_content
        if (!TagGroupGetTagAsTagGroup(content, "Contents", group_content))
            return

        Number num = TagGroupCountTags(group_content)
        Number n
        for (n = 0; n < num; n++) {
            TagGroup child
            if (TagGroupGetIndexedTagAsTagGroup(group_content, n, child))
                VeloxEMDContentWalker(child, names, infos)
        }
    }
}

Number H5IMPORT_select_dataset(String filename, TagGroup infos)
{
    TagGroup dialog = DLGCreateDialog(filename)
    dialog.DLGTableLayout(1, 2, 0)

    dialog.DLGAddElement(DLGCreateLabel( "Select dataset(s) to import:\n"+\
			"Note: only 3D data are images [SliceNo * X_sise * Y_size]. \n"+\
			"Wrong dimension will case error! Check the list below."))\
        .DLGExpand("X") \
        .DLGFill("X") \
        .DLGAnchor("West")
    TagGroup entries
    TagGroup choice = DLGCreateChoice(entries)
    dialog.DLGAddElement(choice) \
        .DLGExpand("X") \
        .DLGFill("X") \
        .DLGAnchor("West")

    entries.DLGAddChoiceItemEntry("<import all>")

    Number num = TagGroupCountTags(infos)
    Number n
    for (n = 0; n < num; n++) {
        String info
        TagGroupGetIndexedTagAsString(infos, n, info)
        entries.DLGAddChoiceItemEntry(info)
    }

    Object frame = alloc(uiframe).init(dialog)
    if (!frame.Pose())
        exit(0)

    return choice.DLGGetValue() - 1
}

number VeloxSetScaleInNM(image &imgout, number &xstep, number &ystep, string &CalString)
{
	if(CalString=="m")
	{
		xstep*=(10**9)  
		ystep*=(10**9)  
	}
	else if(CalString=="µm")
	{
		xstep*=(10**3)  
		ystep*=(10**3)  
	}
	else if(CalString=="nm")
	{
	}
	else if(CalString=="pm")
	{
		xstep/=(10**3)  
		ystep/=(10**3)  
	}
	else
	{
		return(0)
	}

	CalString="nm"
	imgout.setscale(xstep,ystep)
	imgout.setunitstring(CalString)
	//result(CalString+"\n")
	//result(xstep+";"+ystep+"\n")
	return(1)
}

Image VeloxEMDGetOneImage(String filename, String name)
{
    String path = name.PathExtractDirectory(0)
    String base = filename.PathExtractBaseName(0)
    //result("DataName: "+ base +"\nPathName:"+path+"\n")

	// Read the image data from $path\Data
	Image array := h5_read_dataset(filename, path+"Data")
	Number dim1, dim2, dim3
	array.Get3DSize(dim1, dim2, dim3)
	Image Img
	if(dim1 != 1){
		Img := array.slice3(0,0,0, 1,dim2,1, 2,dim3,1, 0,dim1,1).ImageClone()
		//Image img = imageclone(imgRE)
		}
	if(dim1 == 1){
		//Img := array.slice3(0,0,0, 1,dim2,1, 2,dim3,1, 0,dim1,1).slice2(0,0,0, 0,dim2,1, 1,dim3,1).ImageClone()
		Img := array.slice2(0,0,0, 1,dim2,1, 2,dim3,1).ImageClone()
		//Image img = imageclone(imgRE)
		}
	
	// Read the metadata from $path\Metadata
	Image MetaImg := filename.h5_read_dataset(path+"Metadata")
	// For image stack, only the metadata from the first slice is used
	Number metad1, metad2
	MetaImg.GetSize(metad1, metad2)
	Image MetaImg2 := MetaImg.slice2(0,0,0,0,1,1,1,metad2,1)
	String MetadataStr = MetaImg2.img2str()
	// Write metadata string to image tag and return a TagGroup object
	//Tan

	TagGroup MD = alloc(MetaStr2TagGroup).ParseText2ImageTag(MetadataStr, Img)
	//MD.TagGroupopenBrowserwindow("Metadata",0)

	// Get pixel size and set correct calibration
	String UnitXStr, UnitYStr	
	Number PixelSizeX, PixelSizeY
	if(!img.GetStringNote("BinaryResult:PixelUnitX", UnitXStr)) UnitXStr=""
	if(!img.GetStringNote("BinaryResult:PixelUnitX", UnitYStr)) UnitYStr=""
	if(!img.GetNumberNote("BinaryResult:PixelSize:width", \
		PixelSizeX)) PixelSizeX=0
	if(!img.GetNumberNote("BinaryResult:PixelSize:height",\
		PixelSizeY)) PixelSizeY=0
	if(UnitXStr!="" && PixelSizeX!=0 && PixelSizeY!=0)
		VeloxSetScaleInNM(Img, PixelSizeX, PixelSizeY, UnitXStr)
	// Set name using the name of detector
	String ImgName
	if(!img.GetStringNote("BinaryResult:Detector", ImgName)) ImgName=""
	img.SetName(base+"_"+ImgName)

	Return Img
}

Void Main(String filename)
{
	Image VeloxImg
    // Open file
    if(!h5_is_file(filename)) {
        OkDialog("Error opening file " + filename + ".\nNot a HDF5-type file.\n")
        exit(0)
    }

    // Scan file
    TagGroup content = h5_info(filename)
    TagGroup names = NewTagList()
    TagGroup infos = NewTagList()
    VeloxEMDContentWalker(content, names, infos)
	
	Number num = TagGroupCountTags(infos)
    if (num == 0) {
        OkDialog("No supported datasets in file " + filename + ".\n")
        exit(-1)
    } else if (num > 1) {
        Number index = H5IMPORT_select_dataset(filename, infos)
        //result(index)
		if (index < 0) {
            Number n
            for (n = 0; n < num; n++) {
                String name
                TagGroupGetIndexedTagAsString(names, n*3, name)
                VeloxImg := VeloxEMDGetOneImage(filename, name)
                VeloxImg.showimage()
            }
        } else {
            String name
            TagGroupGetIndexedTagAsString(names, index*3, name)
            result(name+"\n")
			VeloxImg := VeloxEMDGetOneImage(filename, name)
			VeloxImg.showimage()
        }
    } else {
        // Single image in file
        String name
        TagGroupGetIndexedTagAsString(names, 0, name)
        VeloxImg := VeloxEMDGetOneImage(filename, name)
        VeloxImg.showimage()
    }
}

//============================= Main ==========================
// Scope needed to avoid memory leaks
{
	String filename
	number count=0
	While (!ShiftDown())
	{
		if (!OpenDialog(Null, "Select Velox file", "*.emd", filename)) exit(0)
		count++
		Result("Import EMD file "+count+":\t"+filename+"\n")
		Main(filename)
			
	}
}