// JINSEOK RYU, Electron Microscopy and Spectroscopy Lab., Seoul National University
// 2019. 11. 4
// remove zlp for low loss EELS (using reflected tail model)

// 2019. 11. 5
// update for dealing with some spikes in spectra

// 2019. 11. 21
// small updates

image SI
getfrontimage(SI)

number sx,sy,sz,xscale,yscale,zscale,xorigin,yorigin,zorigin
string xunit, yunit, zunit

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

number offset
If (!GetNumber("origin of the spectra (= (index-1) * "+zscale+" "+zunit+")", 0, offset)) Exit(0)

number f_length
If (!GetNumber("final length of the spectra (= origin + final length * "+zscale+" "+zunit+")", 0, f_length)) Exit(0)

number offset_ev, last_ev

offset_ev = (offset - 1)*zscale
last_ev = offset_ev + f_length * zscale

result("\n *******************************************")
result("\n energy range of the ourput : [ "+offset_ev+", "+last_ev+" ]")
result("\n *******************************************")

If(!ContinueCancelDialog("you should check the final energy range. continue ?")){
	ShowAlert("the process is terminated by the user",1)
	exit(0)
}

number i, j
image zl_removed := realimage("zero loss removed", 4, sx, sy, f_length)

for(i=0;i<sx;i++){
	for(j=0;j<sy;j++){
		number maxid_x, maxid_y
		image temp, temp_right, temp_left
		temp = slice1(SI, i, j, 0, 2, sz, 1)
		max(temp, maxid_x, maxid_y)
		if((sz-maxid_x) >= maxid_x+1){
			image temp, temp_right, temp_left
			temp = slice1(SI, i, j, 0, 2, sz, 1)
			max(temp, maxid_x, maxid_y)
			temp_right = slice1(SI, i, j, maxid_x, 2, sz-maxid_x, 1)
			temp_left = slice1(SI, i, j, maxid_x, 2, maxid_x+1, -1)
			temp_right.slice1(0, 0, 0, 0, maxid_x+1, 1) -= temp_left
			zl_removed.slice1(i, j, 0, 2, f_length, 1) = temp_right.slice1(offset, 0, 0, 0, f_length, 1)
		}
		
		else{
			result("\n ********************Warning!!!**********************")
			result("\n an error occurs at x: "+i+" y: "+j)
			result("\n it might be caused by some spikes in the spectrum")
			result("\n so this pixel will be replaced by the previous pixel")
			result("\n ****************************************************")
			image temp, temp_right, temp_left
			temp = slice1(SI, i, j-1, 0, 2, sz, 1)
			max(temp, maxid_x, maxid_y)
			temp_right = slice1(SI, i, j-1, maxid_x, 2, sz-maxid_x, 1)
			temp_left = slice1(SI, i, j-1, maxid_x, 2, maxid_x+1, -1)
			temp_right.slice1(0, 0, 0, 0, maxid_x+1, 1) -= temp_left
			zl_removed.slice1(i, j, 0, 2, f_length, 1) = temp_right.slice1(offset, 0, 0, 0, f_length, 1)
		}
	}
}

zl_removed.ImageSetDimensionOrigin( 0, xorigin ) 
zl_removed.ImageSetDimensionScale( 0, xscale ) 
zl_removed.ImageSetDimensionUnitString(0, xunit ) 
zl_removed.ImageSetDimensionOrigin(1, yorigin ) 
zl_removed.ImageSetDimensionScale( 1, yscale ) 
zl_removed.ImageSetDimensionUnitString(1, yunit )
zl_removed.ImageSetDimensionOrigin(2, (offset-1) * zscale ) 
zl_removed.ImageSetDimensionScale(2, zscale ) 
zl_removed.ImageSetDimensionUnitString(2, zunit )

taggroup reftags=SI.imagegettaggroup()
taggroup fittedtags=zl_removed.imagegettaggroup()
taggroupcopytagsfrom(fittedtags, reftags)

zl_removed.showimage()