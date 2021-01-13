// Function to show all hidden images
// D.R.G. Mitchell, adminnospam@dmscripting.com (remove the nospam to make this work)
// version:20121229, v4.0, December 2012, www.dmscripting.com

// This script has been rewritten as previous versions relied upon an exception to halt execution.
// This made it less unsuitable for incorporation into a function - since the previous script halted at the
// exception. Here, no exceptions are generated and the original image stacking order is preserved.


// function to show all hidden images

void showallhiddenimages()
	{
		// variables

		number nodocs,shown, hidden, imgdocid, counter
		number i
		imagedocument imgdoc


		// count the number of image documents - this includes both shown and hidden ones
		// and delete any temporary tag information which may have been left over from last time.
		
		shown=Countdocumentwindowsoftype(5) // images currently shown
		nodocs=countimagedocuments() // all image documents, including hidden ones
		hidden=nodocs-shown // the number of hidden images
		
		deletepersistentnote("Show Hidden Temp")


		// Image documents are indexed (0,1,2 . . . ) in the following manner
		// Front-most shown image (0) then the image shown behind that (1) etc.
		// Hidden images appear at the end of the stacking sequence. The last image in the sequence is the last
		// image to be hidden. In the Window menu, hidden images have a 'h' next to them. 
		// Hidden images have indexes which run top to bottom in this menu. If there are only 3 image documents
		// open and all three are hidden, then the top one in the 'Window' menu will have an index of 0
		// next down 1, and the last will have an index of 2. If there are four image documents open, one is shown
		// and three are hidden, the hidden image documents are indexed 1,2 and 3 ie
		
		// Shown images
		// Front most - index = 0
		// behind front - index = 1
		// rear most - index = n = number of shown imagedocuments-1
		
		// Hidden images
		// First to be hidden - index n+1
		// Last to be hidden - index number of imagedocuments-1

		// Get the IDs of the hidden image documents
		
		for(i=shown; i<nodocs; i++)
			{
				imgdoc=getimagedocument(i)
				imgdocid=imgdoc.imagedocumentgetid()
				setpersistentnumbernote("Show Hidden Temp:Image "+counter,imgdocid)
				counter=counter+1
			}

		// loop to display all the hidden image documents in a stacked sequence. 
		
		number j

		for(j=hidden-1; j>-1; j--) // loop through the number of hidden image documents
			{
				getpersistentnumbernote("Show Hidden Temp:Image "+j,imgdocid)
				imgdoc=getimagedocumentbyid(imgdocid)
				imagedocumentshow(imgdoc)
			}


		// Remove the temporary tags

		deletepersistentnote("Show Hidden Temp")
	}

number nodocs = countdocumentwindowsoftype(5)

if(nodocs<1)
	{
		showalert("Ensure a 1D profile or spectrum is front-most.",2)
		exit(0)
	}

number angle

if(!getnumber("rotation angle", 90, angle)) exit(0)

angle = angle * Pi() / 180

number i
for(i=0;i<nodocs;i++){
	image front := getfrontimage()
	image rotated
	rotated = rotate(front, angle)
	showimage(rotated)
	hideimage(rotated)
	closeimage(front)
}

showallhiddenimages()