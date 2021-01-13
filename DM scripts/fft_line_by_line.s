image front = getfrontimage()
number sx, sy
number n_row
number fy
number i, j, l

getsize(front, sx, sy)
if(!getnumber("number of rows integrated", 0, n_row)) exit(0)
fy = sy - n_row + 1

compleximage linefft := compleximage("fft line by line", 8, sx, fy)
compleximage fft_temp := compleximage("", 8, sx, 1)

for(i=0;i<fy;i++){
	realimage line_integrated := realimage("", 4, sx, 1)
	for(j=0;j<n_row;j++){
		line_integrated = line_integrated + front[i+j, 0, i+1+j, sx]
	}
	
	fft_temp = realFFT(line_integrated)
	
	for(l=0;l<sx;l++){
		setpixel(linefft, l, i, getpixel(fft_temp, l, 0))
		}
}

showimage(linefft)