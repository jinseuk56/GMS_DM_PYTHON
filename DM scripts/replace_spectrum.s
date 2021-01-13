// 20200722
// JINSEOK RYU
// Electron Microscopy and Spectroscopy Lab.
// Seoul National University
// replace an abnormal signal with the previous signal

image front := getfrontimage()

number x, y
if(!getnumber("x position", 0, x)) exit(0)
if(!getnumber("y position", 0, y)) exit(0)

number sx, sy, sz
front.get3dsize(sx, sy, sz)

front.slice1(x, y, 0, 2, sz, 1) = front.slice1(x-1, y, 0, 2, sz, 1)
showimage(front)