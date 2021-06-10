// 20210108
// JINSEOK RYU (jinseuk56@gmail.com)
// Electron Microscopy and Spectroscopy Lab.
// Dept. of Materials Science and Engineering
// Seoul National University
// rolling axes of a 4D image
// (a, b, c, d) -> (c, d, a, b)

image img
if(!getoneimagewithprompt("","select a 4D image", img)) exit(0)
if(img.imagegetnumdimensions() != 4){
	result("\n ensure that the image is 4-dimensional")
	exit(0)
}

number a, b, c, d
a = img.imagegetdimensionsize(0)
b = img.imagegetdimensionsize(1)
c = img.imagegetdimensionsize(2)
d = img.imagegetdimensionsize(3)

string file_name = img.getname()

number sx, sy, dsx, dsy
If (!GetNumber("pixel size of the 1st dimension of the original image", a, sx)) Exit(0)
If (!GetNumber("pixel size of the 2nd dimension of the original image", b, sy)) Exit(0)
If (!GetNumber("pixel size of the 3rd dimension of the original image", c, dsx)) Exit(0)
If (!GetNumber("pixel size of the 4th dimension of the original image", d, dsy)) Exit(0)

image rolled

rolled = img.sliceN(4, 4, 0, 0, 0, 0, 2, dsx, 1, 3, dsy, 1, 0, sx, 1, 1, sy, 1).imageclone()

rolled.setname(file_name+"_axes-rolled")
rolled.showimage()