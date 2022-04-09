// 20210322
// JINSEOK RYU (jinseuk56@gmail.com)
// Electron Microscopy and Spectroscopy Lab.
// Dept. of Materials Science and Engineering
// Seoul National University
// binary data (raw data from EMPAD, Thermo Fisher) -> 4D stack image
// You must know the original shape of 4D-STEM data

Object file_stream
Image img
String filename
Number size_x, size_y, file_ID
If (!OpenDialog(filename)) Exit(0)
number sx, sy, dsx, dsy
If (!GetNumber("pixel size of the 1st dimension", 128, sx)) Exit(0)
If (!GetNumber("pixel size of the 2nd dimension", 128, sy)) Exit(0)
If (!GetNumber("pixel size of the 3rd dimension", 128, dsx)) Exit(0)
If (!GetNumber("pixel size of the 4th dimension", 128, dsy)) Exit(0)
file_ID = OpenFileForReading(filename)
file_stream = NewStreamFromFileReference(file_ID, 1)
img := realimage("raw 4D", 4, dsx, dsy+2, sx, sy)
ImageReadImageDataFromStream(img, file_stream, 2) // 2 -> little-endian
CloseFile(file_ID)

image rolled

rolled := img.sliceN(4, 4, 0, 0, 0, 0, 3, sx, 1, 2, sy, 1, 0, dsx, 1, 1, dsy, 1)

closeimage(img)
rolled.showimage()