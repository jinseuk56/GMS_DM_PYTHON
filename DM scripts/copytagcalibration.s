compleximage ref
image out
string prompt

if(!(getoneimagewithprompt(prompt, "Select the reference image", ref))) exit(0)
if(!(getoneimagewithprompt(prompt, "Select the target image", out))) exit(0)


imagecopycalibrationfrom(out, ref)
taggroup reftags=ref.imagegettaggroup()
taggroup fittedtags=out.imagegettaggroup()
taggroupcopytagsfrom(fittedtags, reftags)

showimage(out)