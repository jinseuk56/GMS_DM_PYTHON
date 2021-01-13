// RYU, J., Material Science and Engineering, Seoul National University
// 2018.3.8
// showing amounts of B-site atom displacement (100 direction perovskite) by a color gradient
// it needs an original image and coordinates of atom peaks created by PPA (from a third-party application by HREM)
// you can change the spot size (2018.3.29)


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


// select the original image and peak position data
if(!(getoneimagewithprompt(prompt, "select the original image", original))) exit(0)
if(!(getoneimagewithprompt(prompt, "select the peak position image", peak_position))) exit(0)

number xos, yos
getsize(original, xos, yos)

number xps, yps
getsize(peak_position, xps, yps)


// get the size of A-site atom matrix
number A_col, A_row
if(!getnumber("number of A-site atoms in one row", 0, A_col)) exit(0)
if(!getnumber("number of A-site atoms in one column", 0, A_row)) exit(0)

number B_col, B_row
B_col = A_col - 1
B_row = A_row - 1

image A_pos := realimage("A site position", 4, A_col, A_row*2)
image B_pos := realimage("B site position", 4, B_col, B_row*2)

number i, j, k, l
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
	//result("\ni = "+i)
	quicksort(B_pos, 0, B_col-1, i)
}
for(i=1;i<A_row*2;i=i+2){
	quicksort(A_pos, 0, A_col-1, i)
}

//showimage(A_pos)
//A_pos.setdisplaytype(7)
//showimage(B_pos)
//B_pos.setdisplaytype(7)

image peak_point_y := rgbimage("peak point image (y-axis displacement)", 4, xos, yos)
image peak_point_x := rgbimage("peak point image (x-axis displacement)", 4, xos, yos)

// images where displacements of B-site atoms are saved
image By_displacement := realimage("B site y-axis displacement", 4, B_col, B_row)
image Bx_displacement := realimage("B site x-axis displacement", 4, B_col, B_row)

ImageDisplay clutDisplay, imgDisp = original.ImageGetImageDisplay(0)
number exportMode = 7

RGBImage rgb_original_y:= ImgDisp.ImageDisplayGetExportImage(exportmode, clutDisplay)
RGBImage rgb_original_x:= ImgDisp.ImageDisplayGetExportImage(exportmode, clutDisplay)

number spot_size
if(!getnumber("spot size", 5, spot_size)) exit(0)

for(i=0;i<A_row-1;i++){
	for(j=0;j<A_col-1;j++){		
		//get 4 A-site atoms
		number x1 = getpixel(A_pos, j, i*2)
		//result("\ni*2 +2 = "+(i*2+2))
		number x2 = getpixel(A_pos, j, i*2+2)
		number x3 = getpixel(A_pos, j+1, i*2)
		number x4 = getpixel(A_pos, j+1, i*2+2)
	
		number y1 = getpixel(A_pos, j, i*2+1)
		number y2 = getpixel(A_pos, j, i*2+3)
		number y3 = getpixel(A_pos, j+1, i*2+1)
		number y4 = getpixel(A_pos, j+1, i*2+3)
		
		number mid_y, mid_x
		mid_y = (y1+y2+y3+y4) / 4.0
		mid_x = (x1+x2+x3+x4) / 4.0
		
		
		number red = 255
		number blue = 255
		rgbnumber white = RGB(255, 255, 255)
		
		/*
		for(k=0;k<spot_size;k++){
			for(l=0;l<spot_size;l++){
				setpixel(peak_point_y, (round(x1)-spot_size+k), (round(y1)-spot_size+l), white)
				setpixel(peak_point_y, (round(x2)-spot_size+k), (round(y2)-spot_size+l), white)
				setpixel(peak_point_y, (round(x3)-spot_size+k), (round(y3)-spot_size+l), white)
				setpixel(peak_point_y, (round(x4)-spot_size+k), (round(y4)-spot_size+l), white)
				setpixel(rgb_original_y, (round(x1)-spot_size+k), (round(y1)-spot_size+l), white)
				setpixel(rgb_original_y, (round(x2)-spot_size+k), (round(y2)-spot_size+l), white)
				setpixel(rgb_original_y, (round(x3)-spot_size+k), (round(y3)-spot_size+l), white)
				setpixel(rgb_original_y, (round(x4)-spot_size+k), (round(y4)-spot_size+l), white)
				
				setpixel(peak_point_x, (round(x1)-spot_size+k), (round(y1)-spot_size+l), white)
				setpixel(peak_point_x, (round(x2)-spot_size+k), (round(y2)-spot_size+l), white)
				setpixel(peak_point_x, (round(x3)-spot_size+k), (round(y3)-spot_size+l), white)
				setpixel(peak_point_x, (round(x4)-spot_size+k), (round(y4)-spot_size+l), white)
				setpixel(rgb_original_x, (round(x1)-spot_size+k), (round(y1)-spot_size+l), white)
				setpixel(rgb_original_x, (round(x2)-spot_size+k), (round(y2)-spot_size+l), white)
				setpixel(rgb_original_x, (round(x3)-spot_size+k), (round(y3)-spot_size+l), white)
				setpixel(rgb_original_x, (round(x4)-spot_size+k), (round(y4)-spot_size+l), white)
			}
		}
		*/

		number B_x = getpixel(B_pos, j, i*2)
		number B_y = getpixel(B_pos, j, i*2+1)
		
		//result("\nx1, y1 = "+x1+", "+y1)
		//result("\nx2, y2 = "+x2+", "+y2)
		//result("\nx3, y3 = "+x3+", "+y3)
		//result("\nx4, y4 = "+x4+", "+y4)
		//result("\nbx, by = "+B_x+", "+B_y)
		
		number displacement_y = mid_y-B_y
		setpixel(By_displacement, j, i, displacement_y)
		
		number displacement_x = mid_x-B_x
		setpixel(Bx_displacement, j, i, displacement_x)
		
		
		//asumming the range of B-site displacement -0.5 ~ 0.5
		number degree_y = displacement_y + 0.5
		number degree_x = displacement_x + 0.5
		
		// red -> above the center, left of the center
		// blue -> below the center, right of the center
		// including extra 8 pixels surrounding the center (corresponding to a B-site atom)
		
		for(k=0;k<(spot_size*2+1);k++){
			for(l=0;l<(spot_size*2+1);l++){
				setpixel(peak_point_y, (round(B_x)-spot_size+k), (round(B_y)-spot_size+l), RGB(red*degree_y, 0, blue*(1-degree_y)))
				setpixel(rgb_original_y, (round(B_x)-spot_size+k), (round(B_y)-spot_size+l), RGB(red*degree_y, 0, blue*(1-degree_y)))
				setpixel(peak_point_x, (round(B_x)-spot_size+k), (round(B_y)-spot_size+l), RGB(red*degree_x, 0, blue*(1-degree_x)))
				setpixel(rgb_original_x, (round(B_x)-spot_size+k), (round(B_y)-spot_size+l), RGB(red*degree_x, 0, blue*(1-degree_x)))
			}
		}
		
	}
}

number max_dis_y, min_dis_y, max_dis_x, min_dis_x

max_dis_y = max(By_displacement)
min_dis_y = min(By_displacement)

max_dis_x = max(Bx_displacement)
min_dis_x = min(Bx_displacement)

result("\nmaximum y-axis displacement of B-site = "+max_dis_y)
result("\nminimum y-axis displacement of B-site = "+min_dis_y)

result("\nmaximum x-axis displacement of B-site = "+max_dis_x)
result("\nminimum x-axis displacement of B-site = "+min_dis_x)

//showimage(B_displacement)
showimage(peak_point_y)
showimage(peak_point_x)
rgb_original_y.setname("B site y-axis displacement")
showimage(rgb_original_y)
rgb_original_x.setname("B site x-axis displacement")
showimage(rgb_original_x)

//showimage(by_displacement)
//by_displacement.setdisplaytype(7)
//showimage(bx_displacement)
//bx_displacement.setdisplaytype(7)