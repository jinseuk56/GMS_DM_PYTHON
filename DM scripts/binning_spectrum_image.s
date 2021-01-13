// 2019. 9. 26
// JINSEOK RYU, Electron Microscopy and Spectroscopy Lab. Seoul National University
// binning a 3D image
// naive version

image binning(image SI, number sx, number sy, number sz, number bx, number by){
	
	image binned

	//result("\ntemporary check for x size "+sx)
	//result("\ntemporary check for y size "+sy)

	if(sx % bx != 0) sx = sx - sx % bx
	if(sy % by != 0) sy = sy - sy % by
	
	//result("\ntemporary check for x size "+sx)
	//result("\ntemporary check for y size "+sy)
	
	
	number new_sx, new_sy
	new_sx = sx / bx
	new_sy = sy / by
	
	//result("\ntemporary check for new x size "+new_sx)
	//result("\ntemporary check for new y size "+new_sy)
	
	binned := realimage("binned SI", 4, new_sx, new_sy, sz)
	
	number i, j, k, l
	
	for(i=0;i<new_sx;i++){
		for(j=0;j<new_sy;j++){
			for(k=bx*i;k<(i+1)*bx;k++){
				for(l=by*j;l<(j+1)*by;l++){
					slice1(binned, i, j, 0, 2, sz, 1) += slice1(SI, k, l, 0, 2, sz, 1)
				}
			}
		}
	}
		
	return binned
	
}
	

image SI
if(!(getoneimage("select a spectrum image", SI))) exit(0)
if(SI.imagegetnumdimensions() != 3){
	showalert("the dimensions of the imported image must be 3", 0)
	exit(0)
}

number sx,sy,sz,xscale,yscale,zscale,xorigin,yorigin,zorigin
string xunit, yunit, zunit
number bx, by
image binned

SI.get3dsize(sx,sy,sz)

xscale = SI.imagegetdimensionscale(0)
yscale = SI.imagegetdimensionscale(1)
zscale = SI.imagegetdimensionscale(2)

xorigin = SI.imagegetdimensionorigin(0)
yorigin = SI.imagegetdimensionorigin(1)
zorigin = SI.imagegetdimensionorigin(2)

xunit = SI.imagegetdimensionunitstring(0)
yunit = SI.imagegetdimensionunitstring(1)
zunit = SI.imagegetdimensionunitstring(2)


if(!getnumber("width of binning", 2, bx)) exit(0)
if(!getnumber("height of binning", 2, by)) exit(0)

binned = binning(SI, sx, sy, sz, bx, by)
binned.setname(bx+"x"+by+" binned SI")

taggroupcopytagsfrom(binned.imagegettaggroup(), SI.imagegettaggroup())

imagecopycalibrationfrom(binned, SI)

//binned.ImageSetDimensionOrigin( 0, xorigin ) 
binned.ImageSetDimensionScale( 0, xscale*bx ) 
//binned.ImageSetDimensionUnitString(0, xunit ) 
//binned.ImageSetDimensionOrigin(1, yorigin ) 
binned.ImageSetDimensionScale( 1, yscale*by ) 
//binned.ImageSetDimensionUnitString(1, yunit )
//binned.ImageSetDimensionOrigin(2, zorigin ) 
//binned.ImageSetDimensionScale( 2, zscale ) 
//binned.ImageSetDimensionUnitString(2, zunit )

showimage(binned)