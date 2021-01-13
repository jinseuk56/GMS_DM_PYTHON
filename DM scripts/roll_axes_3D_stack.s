// 20200519
// JINSEOK RYU
// Electron Microscopy and Spectroscopy Lab.
// Seoul National University
// rolling axes of a 3D image
// (a, b, c) -> (b, c, a)

image img
if(!getoneimagewithprompt("","select a 3D image", img)) exit(0)
if(img.imagegetnumdimensions() != 3){
	result("\n ensure that the image is 3-dimensional")
	exit(0)
}

number ix, iy, iz
img.get3dsize(ix, iy, iz)

number sx, sy, ez
If (!GetNumber("size of the 1st dimension", iy, sx)) Exit(0)
If (!GetNumber("size of the 2nd dimension", iz, sy)) Exit(0)
If (!GetNumber("size of the 3rd dimension", ix, ez)) Exit(0)

image rolled

rolled := img.sliceN(3, 3, 0, 0, 0, 1, sx, 1, 2, sy, 1, 0, ez, 1)

rolled.setname("axes rolled image")
rolled.showimage()