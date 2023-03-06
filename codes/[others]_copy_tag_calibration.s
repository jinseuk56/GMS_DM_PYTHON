image ref
image out
string prompt

if(!(gettwoimageswithprompt("image 0: reference, image 1: target", "Select the reference and target", ref, out))) exit(0)

imagecopycalibrationfrom(out, ref)
taggroup reftags=ref.imagegettaggroup()
taggroup fittedtags=out.imagegettaggroup()
taggroupcopytagsfrom(fittedtags, reftags)

showimage(out)