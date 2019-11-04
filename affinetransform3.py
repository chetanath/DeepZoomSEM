import pyvips
import sys
import random
import os
import time

t0 = time.time()

dothumbnail=False
dodzsave=True
dochopimage=False

thumbnailsize=20000
choptilesize=1024

os.environ["VIPS_PROGRESS"] = "FALSE"
os.environ["VIPS_DISC_THRESHOLD"] = "100m" # DEFAULT=100m
os.environ["VIPS_CONCURRENCY"] = "0" # DEFAULT=0 (16 threads)
os.environ["TMPDIR"] = "/backup/tmp/" # Set temporary directory

# python ~/Documents/affinetransform3.py ~/Documents/MATLAB/imageinformation_010_135.txt ~/Documents/acfer094.dzi
# sudo sysctl -w  kernel.pid_max=4194304
# ulimit -n     : 65535


pyvips.cache_set_max(0)
pyvips.leak_set(True)
infofilepath=sys.argv[1]
outputname='/data/Data/deepzoom/'+os.path.splitext(sys.argv[2])[0]
outputnamefolder=outputname+'_files'


try:
	print('Removing '+outputname+'.dzi')
	os.remove(outputname+'.dzi')
except OSError:
	pass

try:
	print('Removing '+outputnamefolder)
	os.system('rm -fr "%s"' % outputnamefolder)
except OSError:
	pass

if dochopimage:
	try:
		print('Removing '+outputname+'_chop.dzi')
		os.remove(outputname+'.dzi')
	except OSError:
		pass
	try:
		print('Removing '+outputnamefolder)
		os.system('rm -fr "%s"' % outputname+'_chop_files')
	except OSError:
		pass


with open(infofilepath) as f:
   count=sum(1 for _ in f)

print(str(count-1)+" files")

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

xs=0
ys=0

bg = pyvips.Image.black(WIDTH, HEIGHT, bands=1)

maxv8bit=2**8-1
maxv16bit=2**16-1

for ii, line in enumerate(file_object):
    words = line.split(",")
    slope=float(words[1])
    offset=float(words[2])
    print(str(ii+1)+" / "+str(count-1)+" : "+str(slope))
    inputfname=os.path.expanduser(words[0])
    filenamebase, file_extension = os.path.splitext(inputfname)
    image0 = pyvips.Image.new_from_file(inputfname).cast("float")
    image0 = maxv8bit*(slope*(image0+offset))/maxv16bit # +offset
    image0 = image0.cast("uchar")
    #image0 = image0.bandjoin(2**8-1)
    a=float(words[3])
    c=float(words[4])
    b=float(words[5])
    d=float(words[6])
    odx_=float(words[7])
    ody_=float(words[8])
    bordercutpixels=min(int(round(0.5*image0.width)-1),int(round(max(abs((a+b-1)*image0.width),abs((c+d-1)*image0.width))))) #int(round(0.005*image0.width)); # pixels to cut from border before insert
    #print('Border cut pixels='+str(bordercutpixels))
    #xpe=int(image0.width*a+image0.height*b+odx_)
    #ype=int(image0.width*c+image0.height*d+ody_)
    image0 = image0.affine([a,b,c,d],oarea=[bordercutpixels,bordercutpixels,image0.width-2*bordercutpixels,image0.height-2*bordercutpixels],interpolate=pyvips.Interpolate.new("bicubic"))
    bg=bg.insert(image0, int(odx_-xMin+bordercutpixels), int(ody_-yMin+bordercutpixels), background=maxv8bit)
    del image0	

print("Mosaic construction complete")

if dothumbnail:

    bg.thumbnail_image(thumbnailsize).write_to_file(os.path.expanduser(outputname+".png"))

    print("Thumbnail file write complete")

if dodzsave:

    bg.dzsave(outputname,suffix='.jpg[Q=100]')

    print("dzsave complete")

if dochopimage:
        
    bg.dzsave(outputname+'_chop',suffix='.jpg[Q=100]',depth='one',tile_size=choptilesize,overlap=0)
        
del bg
t1 = time.time()

print("{0} seconds".format(t1 - t0))
