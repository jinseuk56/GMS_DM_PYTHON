void swap(image &arr, number i, number j, number row){
	number temp = getpixel(arr, i, row)
	setpixel(arr, i, row, getpixel(arr, j, row))
	setpixel(arr, j, row, temp)
	return
}

void partition(image &arr, number low, number high, number row, number &pi){
	number pivot = getpixel(arr, high, row)
	
	number i = low - 1
	number j
	
	image temp = arr
	for(j=low;j<=high-1;j++)
	{
		if(getpixel(temp, j, row) <= pivot)
		{
			i = i + 1
			swap(temp, i, j, row)
			swap(temp, i, j, row-1)
		}
	}
	swap(temp, i+1, high, row)
	swap(temp, i+1, high, row-1)
	arr = temp
	pi = i + 1
	return
}

void quicksort(image &arr, number low, number high, number row)
{
	image temp = arr
	if(low < high)
	{
		number pi
		partition(temp, low, high, row, pi)
		//result("\n"+pi)
		quicksort(temp, low, pi-1, row)
		quicksort(temp, pi+1, high, row)
	}
	arr = temp
	return
}

image original
image peak_position
string prompt

if(!(getoneimagewithprompt(prompt, "select the original image", original))) exit(0)
if(!(getoneimagewithprompt(prompt, "select the peak position image", peak_position))) exit(0)

number xos, yos
getsize(original, xos, yos)

number xps, yps
getsize(peak_position, xps, yps)

number p, q, temp
for(p=0;p<xps;p++){
	for(q=0;q<yps;q++){
		temp = round(getpixel(peak_position, p, q))
		setpixel(peak_position, p, q, temp)
	}
}

number A_col, A_row
if(!getnumber("number of A-site atoms in one row", 0, A_col)) exit(0)
if(!getnumber("number of A-site atoms in one column", 0, A_row)) exit(0)

number B_col, B_row
B_col = A_col - 1
B_row = A_row - 1

image A_pos := realimage("A site position", 4, A_col, A_row*2)
image B_pos := realimage("B site position", 4, B_col, B_row*2)

number i, j, k
i = 2
j = 0
k = 0
while(1)
{

	for(j=0;j<A_col;j++)
	{
		//result("l = "+l+"\n")
		number x = getpixel(peak_position, 0, j+i)
		number y = getpixel(peak_position, 1, j+i)
		//result("(x, y) = "+"("+x+" ,"+y+")\n")
		setpixel(A_pos, j, k, x)
		setpixel(A_pos, j, k+1, y)
	}
	i = i + A_col
	
	if(i == yps) break

	for(j=0;j<B_col;j++)
	{
		number x = getpixel(peak_position, 0, j+i)
		number y = getpixel(peak_position, 1, j+i)
		setpixel(B_pos, j, k, x)
		setpixel(B_pos, j, k+1, y)
	}
	i = i + B_col
	k = k + 2
}
result("\n***********************************")

for(i=1;i<B_row*2;i=i+2){
	result("\ni = "+i)
	quicksort(B_pos, 0, B_col-1, i)
}
for(i=1;i<A_row*2;i=i+2){
	quicksort(A_pos, 0, A_col-1, i)
}

//showimage(A_pos)
//A_pos.setdisplaytype(7)
//showimage(B_pos)
//B_pos.setdisplaytype(7)

image peak_point := rgbimage("peak point image", 4, xos, yos)

ImageDisplay clutDisplay, imgDisp = original.ImageGetImageDisplay(0)
number exportMode = 7

RGBImage rgb_original:= ImgDisp.ImageDisplayGetExportImage(exportmode, clutDisplay)


for(i=0;i<A_row-1;i++){
	for(j=0;j<A_col-1;j++){		
		number x1 = getpixel(A_pos, j, i*2)
		//result("\ni*2 +2 = "+(i*2+2))
		number x2 = getpixel(A_pos, j, i*2+2)
		number x3 = getpixel(A_pos, j+1, i*2)
		number x4 = getpixel(A_pos, j+1, i*2+2)
	
		number y1 = getpixel(A_pos, j, i*2+1)
		number y2 = getpixel(A_pos, j, i*2+3)
		number y3 = getpixel(A_pos, j+1, i*2+1)
		number y4 = getpixel(A_pos, j+1, i*2+3)
		
		rgbnumber white = RGB(255, 255, 255)
		setpixel(peak_point, x1, y1, white)
		setpixel(peak_point, x2, y2, white)
		setpixel(peak_point, x3, y3, white)
		setpixel(peak_point, x4, y4, white)
		setpixel(rgb_original, x1, y1, white)
		setpixel(rgb_original, x2, y2, white)
		setpixel(rgb_original, x3, y3, white)
		setpixel(rgb_original, x4, y4, white)
		
		number mid_y
		mid_y = (y1+y2+y3+y4) / 4.0
		
		rgbnumber red = RGB(255, 0, 0)
		rgbnumber blue = RGB(0, 0, 255)
		rgbnumber green = RGB(0, 255, 0)
		
		number B_x = getpixel(B_pos, j, i*2)
		number B_y = getpixel(B_pos, j, i*2+1)
			
		if(B_y > mid_y){
			setpixel(peak_point, B_x, B_y, red)
			/*setpixel(peak_point, B_x-1, B_y, red)
			setpixel(peak_point, B_x+1, B_y, red)
			setpixel(peak_point, B_x, B_y-1, red)
			setpixel(peak_point, B_x-1, B_y-1, red)
			setpixel(peak_point, B_x+1, B_y-1, red)
			setpixel(peak_point, B_x, B_y+1, red)
			setpixel(peak_point, B_x-1, B_y+1, red)
			setpixel(peak_point, B_x+1, B_y+1, red)*/
			setpixel(rgb_original, B_x, B_y, red)
			setpixel(rgb_original, B_x-1, B_y, red)
			setpixel(rgb_original, B_x+1, B_y, red)
			setpixel(rgb_original, B_x, B_y-1, red)
			setpixel(rgb_original, B_x-1, B_y-1, red)
			setpixel(rgb_original, B_x+1, B_y-1, red)
			setpixel(rgb_original, B_x, B_y+1, red)
			setpixel(rgb_original, B_x-1, B_y+1, red)
			setpixel(rgb_original, B_x+1, B_y+1, red)
		}
		if(B_y < mid_y){
			setpixel(peak_point, B_x, B_y, blue)
			/*setpixel(peak_point, B_x-1, B_y, blue)
			setpixel(peak_point, B_x+1, B_y, blue)
			setpixel(peak_point, B_x, B_y-1, blue)
			setpixel(peak_point, B_x-1, B_y-1, blue)
			setpixel(peak_point, B_x+1, B_y-1, blue)
			setpixel(peak_point, B_x, B_y+1, blue)
			setpixel(peak_point, B_x-1, B_y+1, blue)
			setpixel(peak_point, B_x+1, B_y+1, blue)*/
			setpixel(rgb_original, B_x, B_y, blue)
			setpixel(rgb_original, B_x-1, B_y, blue)
			setpixel(rgb_original, B_x+1, B_y, blue)
			setpixel(rgb_original, B_x, B_y-1, blue)
			setpixel(rgb_original, B_x-1, B_y-1, blue)
			setpixel(rgb_original, B_x+1, B_y-1, blue)
			setpixel(rgb_original, B_x, B_y+1, blue)
			setpixel(rgb_original, B_x-1, B_y+1, blue)
			setpixel(rgb_original, B_x+1, B_y+1, blue)
		}
		if(B_y == mid_y){
			setpixel(peak_point, B_x, B_y, green)
			/*setpixel(peak_point, B_x-1, B_y, green)
			setpixel(peak_point, B_x+1, B_y, green)
			setpixel(peak_point, B_x, B_y-1, green)
			setpixel(peak_point, B_x-1, B_y-1, green)
			setpixel(peak_point, B_x+1, B_y-1, green)
			setpixel(peak_point, B_x, B_y+1, green)
			setpixel(peak_point, B_x-1, B_y+1, green)
			setpixel(peak_point, B_x+1, B_y+1, green)*/			
			setpixel(peak_point, B_x, B_y, green)
			setpixel(rgb_original, B_x, B_y, green)
			setpixel(rgb_original, B_x-1, B_y, green)
			setpixel(rgb_original, B_x+1, B_y, green)
			setpixel(rgb_original, B_x, B_y-1, green)
			setpixel(rgb_original, B_x-1, B_y-1, green)
			setpixel(rgb_original, B_x+1, B_y-1, green)
			setpixel(rgb_original, B_x, B_y+1, green)
			setpixel(rgb_original, B_x-1, B_y+1, green)
			setpixel(rgb_original, B_x+1, B_y+1, green)
		}
	}
}

showimage(peak_point)
showimage(rgb_original)
