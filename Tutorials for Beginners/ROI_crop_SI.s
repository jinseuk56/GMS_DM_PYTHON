// 20201207
// JINSEOK RYU (jinseuk56@gmail.com)
// Electron Microscopy and Spectroscopy Lab.
// Dept. of Materials Science and Engineering
// Seoul National University
// crop a spectrum image (STEM-EELS or 4D-STEM)
// before cropping, there must be a ROI in the spectrum image

image front := getfrontimage()
number top, left, bottom, right

front.getselection(top, left, bottom, right)
result("\n the ROI information")
result("\n top left bottom right")
result("\n "+top+" "+left+" "+bottom+" "+right)

image cropped
cropped := front[top, left, bottom, right]

showimage(cropped)