compleximage front := getfrontimage()
number sx, sy

getsize(front, sx, sy)
result("\n"+sx+" "+sy)

image out := realimage("log10(abs(FFT))", 4, sx, sy)

number i, j

for(i=0;i<sx;i++){
	for(j=0;j<sy;j++){
		setpixel(out, i, j, log10(abs(getpixel(front, i, j))))
		//result("\n"+getpixel(front, i, j))
	}
}

imagecopycalibrationfrom(out, front)
taggroup reftags=front.imagegettaggroup()
taggroup fittedtags=out.imagegettaggroup()
taggroupcopytagsfrom(fittedtags, reftags)

showimage(out)