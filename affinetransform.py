import pyvips
import sys
import random
import os
import shutil
import time


#  python affinetransform.py ~/Documents/MATLAB/imageinformation_040_050.txt ~/Documents/acfer094.dzi

# size of output image
#WIDTH = 8000
#HEIGHT = 3000


#infofilepath=os.path.expanduser("~/Documents/MATLAB/imageinformation_054_052.txt")

pyvips.cache_set_max(0)
pyvips.leak_set(True)
infofilepath=sys.argv[1]
outputname=os.path.splitext(sys.argv[2])[0]
outputnamefolder=outputname+'_files'

print(outputname)

try:
	os.remove(outputname)
except OSError:
	pass

try:
	shutil.rmtree(outputnamefolder)
except OSError:
	pass


count = 0
for line in open(infofilepath).xreadlines(  ): count += 1
print count-1,' files'
file_object  = open(infofilepath, "r")
words = file_object.readline()
words = words.split(",")

for ii, line in enumerate(file_object):
#    print(ii)
    words = line.split(",")
    inputfname=os.path.expanduser(words[0])
    filenamebase, file_extension = os.path.splitext(inputfname)
    if not os.path.isfile(filenamebase+".v"):
        print(inputfname) 
    	image0 = pyvips.Image.new_from_file(inputfname)
    	image0 = image0.scaleimage()
    	image0 = image0.cast("uchar")
   	image0 = image0.bandjoin(1)
    	image0.write_to_file(filenamebase+".v")

file_object  = open(infofilepath, "r")
words = file_object.readline()
words = words.split(",")
xMin=float(words[0])
xMax=float(words[1])
yMin=float(words[2])
yMax=float(words[3])

WIDTH=int(round(xMax-xMin))
HEIGHT=int(round(yMax-yMin))

print("%d x %d pixels = %.2f gigapixels" %(WIDTH,HEIGHT,WIDTH*HEIGHT/1e9))

for ii, line in enumerate(file_object):
    print(ii)
    t0 = time.time()
    words = line.split(",")
    inputfname=os.path.expanduser(words[0])
    filenamebase, file_extension = os.path.splitext(inputfname)
    image0 = pyvips.Image.new_from_file(filenamebase+".v", access="sequential")
    odx_=float(words[5])-xMin
    ody_=float(words[6])-yMin
    image0 = image0.affine([float(words[1]),float(words[2]),float(words[3]),float(words[4])],odx=odx_,ody=ody_,oarea=[0,0,WIDTH,HEIGHT], background=[0,0])
    [gray1, alpha] = image0.bandsplit()
    if ii == 0:
        image = gray1.copy()
    else:
        image = alpha.ifthenelse(gray1, image)
    #image.write_to_file(os.path.expanduser("~/Documents/x.v"))
    #t1 = time.time()
    #print("{0} seconds".format(t1 - t0))
#print(image.bands)
print("mosaic construction complete")
#[gray1, alpha] = image.bandsplit()
#print("bandsplit complete")

# Image & thumbnail write:
image.write_to_file(os.path.expanduser("~/Documents/x.v"))
print("file write complete")
imageth=image.thumbnail(3000)
imageth.write_to_file("x_th.jpg")
print("thumbnail write complete")

#print("write buffer complete")
#gray1 = gray1.scaleimage()
#gray1 = gray1.cast("uchar")
#print("8-bit conversion complete")
#gray1.write_to_file("x.v")
#print("file write complete")
#image00 =  pyvips.Image.new_from_file("x.v")
#print("file read complete")
#print("bandsplit complete")
#gray1=gray1.cast("uchar")
#gray1 = gray1.resize(0.1, kernel = "linear")

#grayth=gray1.thumbnail_image(gray1.width/10)
#grayth.write_to_file(outputfilename)


gray1.dzsave(outputname)
print("dzsave complete")
