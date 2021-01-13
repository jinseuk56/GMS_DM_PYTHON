// 20200625
// JINSEOK RYU
// Electron Microscopy and Spectroscopy Lab.
// Seoul National University
// dimensions (a, b, 1) -> (a, b)

image front := getfrontimage()
number sx, sy
front.getsize(sx, sy)

showimage(front.slice2(0, 0, 0, 0, sx, 1, 1, sy, 1))