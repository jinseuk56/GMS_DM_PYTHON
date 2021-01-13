// EELS background fitting with chi square


// selecting a reference spectrum and observed spectrum
// selected spectra are saved to tmp_refspec and tmp_obsspec respectively
image tmp_refspec, tmp_obsspec
string prompt

if(!(getoneimagewithprompt(prompt, "Select a reference spectrum", tmp_refspec))) exit(0)
if(!(getoneimagewithprompt(prompt, "Select a observed sepctrum", tmp_obsspec))) exit(0)

// edge slices of each spectrum are saved to refspec and obsspec respectively
// each partial spectrum corresponding to the ROI area which is already assigned is saved to tmp_ref and tmp_obs respectively
image refspec, obsspec, tmp_ref, tmp_obs
number top, left, bottom, right

refspec := tmp_refspec{2}
obsspec := tmp_obsspec{2}
getselection(obsspec, top, left, bottom, right) // getting coordinates of ROI
//result(top+" "+left+" "+bottom+" "+right+"\n")
tmp_ref := refspec[top, left, bottom, right]
tmp_obs := obsspec[top, left, bottom, right]
//showimage(tmp_ref)
//showimage(tmp_obs)


// calculating the sumy of intensities in tmp_ref and obs_ref
// getting the intensity avgs of tmp_ref and obs_ref
// finally obtaining the avg ratio (= ref_avg / obs_avg)
number xsize, ysize
number sum_ref, i, avg_ref
number sum_obs, j, avg_obs
number rat_avg

getsize(tmp_ref, xsize, ysize)

i=0
while(i<xsize)
{
sum_ref += getpixel(tmp_ref, i, 0)
sum_obs += getpixel(tmp_obs, i, 0)
i++
}
avg_ref = sum_ref/xsize
avg_obs = sum_obs/xsize
//result(avg_ref+"\n"+avg_obs+"\n")
rat_avg = avg_ref/avg_obs
result("average ratio is "+rat_avg+"\n")

// obs_tmp variation(range) -> +- 10% of rat_avg & number of steps = (range/increment)*2
// calculating chi squares between ref_tmp and obs_tmp multiplied by avg_ratio +- range with increment value
// chi squares are saved to savechi
image savechi, ingr_obs
number rat_tmp, rat_var, k, l, chi_comp
number range, percent, increment, numofsteps
percent = 10
range = round(rat_avg/(100/percent))

result(range)
increment = 0.001
numofsteps = (range/increment)*2

savechi := realimage("chisquare vs ratio",4,numofsteps,1)
rat_var = -range

for(k=0;k<numofsteps;k++)
{
rat_tmp = rat_avg + rat_var
ingr_obs = tmp_obs * rat_tmp
rat_var += increment

chi_comp = 0
l = 0
while(l<xsize)
{
chi_comp += (getpixel(tmp_ref, l, 0)-getpixel(ingr_obs, l, 0))**2/getpixel(tmp_ref, l, 0)**2
l++
}

setpixel(savechi, k, 0, chi_comp)
}
showimage(savechi) //showing the value of chi square depending on the step


// finding the minimum chi square and its position within savechi
// by using the x-position of the minimum chi square, obtaining min_chi_ratio 
number min_x, min_y, min_chi_ratio

min(savechi, min_x, min_y)

min_chi_ratio = rat_avg -1 + increment*min_x
result("When the observed spectrum is multiplied by "+min_chi_ratio+", the minimum value of chi square is obtained\n")

// displaying the image and copy across the calibrations and tag groups
// adding the fitted oberved spectrum to the reference spectrum
image fittedspec
string resultname
fittedspec = refspec
getname(fittedspec, resultname)
setname(fittedspec, "Result")
showimage(fittedspec)
updateimage(fittedspec)

imagecopycalibrationfrom(fittedspec, tmp_refspec)
taggroup reftags=tmp_refspec.imagegettaggroup()
taggroup fittedtags=fittedspec.imagegettaggroup()
taggroupcopytagsfrom(fittedtags, reftags)

imagedisplay imgdisp = fittedspec.imagegetimagedisplay(0)
imgdisp.lineplotimagedisplaysetlegendshown(1)
object sliceid=imgdisp.imagedisplaygetsliceidbyindex(0)
imgdisp.imagedisplaysetslicelabelbyid(sliceid, "reference")
imgdisp.imagedisplayaddimage(obsspec*min_chi_ratio, "fitted to reference")
