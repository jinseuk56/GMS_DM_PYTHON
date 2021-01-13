image front := getfrontimage()
number top, left, bottom, right

getselection(front, top, left, bottom, right)

number xsize, ysize

result(top+" "+left+" "+bottom+" "+right+"\n")

xsize = right - left
ysize = bottom - top

number i, j
number x, y

for(i=0;i<xsize;i++)
{
	for(j=0;j<ysize;j++)
	{	
		x = i + left
		y = j + top
		result(x+" "+y+"\n")
	}
}