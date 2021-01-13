compleximage front := getfrontimage()
number sx, sy

getsize(front, sx, sy)
result("\n"+sx+" "+sy)

image out_tmp := realimage("abs(FFT)", 4, sx, sy)

number i, j

for(i=0;i<sx;i++){
	for(j=0;j<sy;j++){
		setpixel(out_tmp, i, j, abs(getpixel(front, i, j)))
		//result("\n"+getpixel(front, i, j))
	}
}

out_tmp = out_tmp / max(out_tmp)

image out := realimage("(normalize(abs(FFT)))**3", 4, sx, sy)
for(i=0;i<sx;i++){
	for(j=0;j<sy;j++){
		setpixel(out, i, j, getpixel(out_tmp, i, j)**3)
		//result("\n"+getpixel(front, i, j))
	}
}


imagecopycalibrationfrom(out, front)
taggroup reftags=front.imagegettaggroup()
taggroup fittedtags=out.imagegettaggroup()
taggroupcopytagsfrom(fittedtags, reftags)

showimage(out)