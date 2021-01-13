// J. S. RYU, Electron Microscopy and Spectroscopy Lab., Material Science and Engineering, Seoul National University
// 2018. 4. 30

string path
number file_ref
number value
string line

number width
number height

if(!getnumber("width", 0, width)) exit(0)
if(!getnumber("height", 0, height)) exit(0)

image imported := realimage("imported", 4, width, height)
number i, j

if(!opendialog(path)) exit(0)
file_ref = openfileforreading(path)

for(i=0;i<height;i++){
	for(j=0;j<width;j++){
		readfileline(file_ref, line)
		value = line.val()
		setpixel(imported, j, i, value)
	}
}

showimage(imported)