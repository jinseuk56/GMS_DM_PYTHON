// JINSEOK RYU, Electron Microscopy and Spectroscopy Lab., Seoul National University
// 2019. 11. 21
// crop all the spectra for a spectrum image (3D)

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

number start, end

if(!getnumber("the start point of the energy range (eV)", 0, start)) exit(0)
if(!getnumber("the end point of the energy range (eV)", 0, end)) exit(0)

number start_ind, end_ind

start_ind = round((start - zorigin) / zscale)
end_ind = round((end - zorigin) / zscale) - 1

image cropped
cropped := imageclone(SI[0, 0, start_ind, sx, sy, end_ind+1])

showimage(cropped)
cropped.setname(SI.getname()+"_cropped")