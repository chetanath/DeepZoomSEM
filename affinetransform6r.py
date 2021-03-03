import pyvips
import sys
import random
import os
import time
import logging
#import numpy as np

logging.basicConfig(level=logging.WARNING)

t0 = time.time()

os.environ["VIPS_PROGRESS"] = "FALSE"
os.environ["VIPS_DISC_THRESHOLD"] = "1m" # DEFAULT=100m
os.environ["VIPS_CONCURRENCY"] = "1" # DEFAULT=0 (16 threads)
#temporarydirectory = "/var/TMP/" # Set temporary directory
temporarydirectory = "/backup/tmp/" # Set temporary directory
os.environ["TMPDIR"] = temporarydirectory

pyvips.cache_set_max(0)
pyvips.leak_set(True)

infofilepath=sys.argv[1]
outputdir=sys.argv[2]

thisrow=int(sys.argv[3])

maxv8bit=2**8-1
maxv16bit=2**16-1

f  = open(infofilepath, "r+")
lines = [line for line in f.readlines()]
f.close()

vals = lines[0].split(",")
xMin=float(vals[0])
xMax=float(vals[1])
yMin=float(vals[2])
yMax=float(vals[3])

WIDTH=int(round(xMax-xMin))
HEIGHT=int(round(yMax-yMin))


print("%d x %d pixels = %.2f gigapixels" %(WIDTH,HEIGHT,WIDTH*HEIGHT/1e9))

nrows = int(lines[1])
nn=2
rowimagelist=[]

#deltay=[]
compositelist=[]
compositenumber=0

for jj in range(nrows):
    words = lines[nn].split(",")
    nimages = int(words[0])
    YL_l=int(words[1])
    YL_width=int(words[2])
    imagelist=[]
    
    for ii in range(nn,nn+nimages):
        words = lines[ii+1].split(",")
        offset=float(words[1])
        inputfname=os.path.expanduser(words[0])
        a=float(words[2])
        c=float(words[3])
        b=float(words[4])
        d=float(words[5])
        odx=float(words[6])
        ody=float(words[7])
        if jj==thisrow:
            print(inputfname)
            filenamebase, file_extension = os.path.splitext(inputfname)
            image0 = pyvips.Image.new_from_file(inputfname)
            image0 = image0.cast("float")
            image0 = image0+offset
            image0 = maxv8bit*image0/maxv16bit
            image0 = image0.cast("uchar")
            image0 = image0.bandjoin(255).copy(interpretation="b-w")
            image0 = image0.affine([a,b,c,d],odx=odx,ody=ody,oarea=[round(xMin),YL_l,WIDTH,YL_width],interpolate=pyvips.Interpolate.new("bicubic"),background=[0,0])
            imagelist.insert(0,image0)
    
    if jj==thisrow:
        RR=imagelist[0].composite(imagelist[1:], "over")
        imagename=temporarydirectory+outputdir+"/rowimage"+"{:03d}".format(jj)+".v"
        RR[0].write_to_file(imagename)
        print("Row = "+"{:03d}".format(jj))
    nn=nn+nimages+1

        
t1 = time.time()

print("")
print("")
print("{0} seconds".format(t1 - t0))
print("")
print("")
