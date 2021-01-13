// 2018. 02. 27
// JINSEOK RYU, Electron Microscopy and Spectroscopy Lab. Seoul National University
// quicksort funtion for one line

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
		}
	}
	swap(temp, i+1, high, row)
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
		result("\n"+pi)
		quicksort(temp, low, pi-1, row)
		quicksort(temp, pi+1, high, row)
	}
	arr = temp
	return
}