// RYU, J., Material Science and Engineering, Seoul National University
// 2018.3.14
// showing amounts of B-site atom displacement (110 direction perovskite) by a color gradient
// it needs an original image and coordinates of atom peaks created by PPA (from a third-party application by HREM)


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
image peak_position_original
string prompt

// select the original image and peak position data
if(!(getoneimagewithprompt(prompt, "select the original image", original))) exit(0)
if(!(getoneimagewithprompt(prompt, "select the peak position image", peak_position_original))) exit(0)

number xos, yos
getsize(original, xos, yos)

number xps, yps
getsize(peak_position_original, xps, yps)

image peak_position_transpose := realimage("peak position transpose", 4, yps-2, xps-1)

number p, q

for(p=0;p<yps-2;p++){
	number x = getpixel(peak_position_original, 0, p+2)
	number y = getpixel(peak_position_original, 1, p+2)
	setpixel(peak_position_transpose, p, 0, x)
	setpixel(peak_position_transpose, p, 1, y)
}

quicksort(peak_position_transpose, 0, yps-3, 1)

image peak_position := realimage("peak position rearranged", 4, xps-1, yps-2)
for(p=0;p<yps-2;p++){
	setpixel(peak_position, 0, p, getpixel(peak_position_transpose, p, 0))
	setpixel(peak_position, 1, p, getpixel(peak_position_transpose, p, 1))
}

//showimage(peak_position)
//peak_position.setdisplaytype(7)

// get the size of A-site atom matrix
number col, row
if(!getnumber("number of atoms in one row", 0, col)) exit(0)

row = (yps-2) / col
result("\nrow = "+row+" col = "+col)
image pos := realimage("rearrangement", 4, col, row*2)

number i, j, k
i = 0
j = 0
k = 0
while(1)
{

	for(j=0;j<col;j++)
	{
		//result("l = "+l+"\n")
		number x = getpixel(peak_position, 0, j+i)
		number y = getpixel(peak_position, 1, j+i)
		//result("\n(x, y) = "+"("+x+" ,"+y+")")
		setpixel(pos, j, k, y)
		setpixel(pos, j, k+1, x)
	}
	i = i + col
	if(i == yps-2) break
	
	k = k + 2
}
result("\n***********************************")

for(i=1;i<row*2;i=i+2){
	//result("\ni = "+i)
	quicksort(pos, 0, col-1, i)
}

showimage(pos)
pos.setdisplaytype(7)

image peak_point := rgbimage("peak point image", 4, xos, yos)
image B_displacement := realimage("B site displacement", 4, col, round(row/2))

ImageDisplay clutDisplay, imgDisp = original.ImageGetImageDisplay(0)
number exportMode = 7

RGBImage rgb_original:= ImgDisp.ImageDisplayGetExportImage(exportmode, clutDisplay)


for(i=0;i<col;i++){
	for(j=0;j<row-2;j=j+2){
		result("\ni = "+i+" j = "+j)
		number x1 = getpixel(pos, i, 2*j+1)
		number x2 = getpixel(pos, i, 2*j+5)

		number y1 = getpixel(pos, i, 2*j)
		number y2 = getpixel(pos, i, 2*j+4)
		
		number mid_y
		mid_y = (y1+y2) / 2.0
		
		
		number red = 255
		number blue = 255
		rgbnumber white = RGB(255, 255, 255)
		
		setpixel(peak_point, round(x1), round(y1), white)
		setpixel(peak_point, round(x2), round(y2), white)

		setpixel(rgb_original, round(x1), round(y1), white)
		setpixel(rgb_original, round(x2), round(y2), white)

		
		number B_x = getpixel(pos, i, 2*j+3)
		number B_y = getpixel(pos, i, 2*j+2)
		
		
		result("\nx1 = "+x1+" y1 = "+y1)
		result("\nx2 = "+x2+" y2 = "+y2)
		result("\nb_x = "+B_x+" b_y = "+B_y)
		
		number displacement = mid_y-B_y
		setpixel(B_displacement, i, j/2.0, displacement)
		
		//asumming the range of B-site displacement -3.0 ~ 3.0
		number degree = ( displacement + 3 ) / 6
		
		// red -> above the center
		// blue -> below the center
		// including extra 8 pixels surrounding the center (corresponding to a B-site atom)
		setpixel(peak_point, round(B_x), round(B_y), RGB(red*degree, 0, blue*(1-degree)))
		setpixel(peak_point, round(B_x-1), round(B_y), RGB(red*degree, 0, blue*(1-degree)))
		setpixel(peak_point, round(B_x+1), round(B_y), RGB(red*degree, 0, blue*(1-degree)))
		setpixel(peak_point, round(B_x), round(B_y-1), RGB(red*degree, 0, blue*(1-degree)))
		setpixel(peak_point, round(B_x-1), round(B_y-1), RGB(red*degree, 0, blue*(1-degree)))
		setpixel(peak_point, round(B_x+1), round(B_y-1), RGB(red*degree, 0, blue*(1-degree)))
		setpixel(peak_point, round(B_x), round(B_y+1), RGB(red*degree, 0, blue*(1-degree)))
		setpixel(peak_point, round(B_x-1), round(B_y+1), RGB(red*degree, 0, blue*(1-degree)))
		setpixel(peak_point, round(B_x+1), round(B_y+1), RGB(red*degree, 0, blue*(1-degree)))
		setpixel(rgb_original, round(B_x), round(B_y), RGB(red*degree, 0, blue*(1-degree)))
		setpixel(rgb_original, round(B_x-1), round(B_y), RGB(red*degree, 0, blue*(1-degree)))
		setpixel(rgb_original, round(B_x+1), round(B_y), RGB(red*degree, 0, blue*(1-degree)))
		setpixel(rgb_original, round(B_x), round(B_y-1), RGB(red*degree, 0, blue*(1-degree)))
		setpixel(rgb_original, round(B_x-1), round(B_y-1), RGB(red*degree, 0, blue*(1-degree)))
		setpixel(rgb_original, round(B_x+1), round(B_y-1), RGB(red*degree, 0, blue*(1-degree)))
		setpixel(rgb_original, round(B_x), round(B_y+1), RGB(red*degree, 0, blue*(1-degree)))
		setpixel(rgb_original, round(B_x-1), round(B_y+1), RGB(red*degree, 0, blue*(1-degree)))
		setpixel(rgb_original, round(B_x+1), round(B_y+1), RGB(red*degree, 0, blue*(1-degree)))
	}
}

number max_dis, min_dis

max_dis = max(B_displacement)
min_dis = min(B_displacement)

result("\nmaximum displacement of B-site = "+max_dis)
result("\nminimum displacement of B-site = "+min_dis)

//showimage(B_displacement)
showimage(peak_point)
showimage(rgb_original)