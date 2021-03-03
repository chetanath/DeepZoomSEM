import pyvips
import sys
import random
import os
import time
import logging
import glob
import shutil
#import numpy as np

logging.basicConfig(level=logging.WARNING)

t0 = time.time()

dothumbnail=True
dodzsave=True
dopngsave=False
dochopimage=False

#baseoutputdir='/backup/'
#baseoutputdir='/home/ogliore/Data/deepzoom/'
baseoutputdir=os.path.expanduser('~/Data/deepzoom/')

thumbnailsize=5000
choptilesize=1024

os.environ["VIPS_PROGRESS"] = "FALSE"
os.environ["VIPS_DISC_THRESHOLD"] = "1m" # DEFAULT=100m
os.environ["VIPS_CONCURRENCY"] = "0" # DEFAULT=0 (16 threads)
#temporarydirectory = "/var/TMP/" # Set temporary directory
temporarydirectory = "/backup/tmp/" # Set temporary directory
os.environ["TMPDIR"] = temporarydirectory



# python ~/Documents/affinetransform3.py ~/Documents/MATLAB/imageinformation_010_135.txt ~/Documents/acfer094.dzi
# sudo sysctl -w  kernel.pid_max=4194304
# ulimit -n     : 65535

pyvips.cache_set_max(0)
pyvips.leak_set(True)

outputname=baseoutputdir+os.path.splitext(sys.argv[1])[0]
outputext=os.path.splitext(sys.argv[1])[1]
outputnamefolder=outputname+'_files'

if dodzsave:
        try:
            print('Removing '+baseoutputdir+sys.argv[1]+'.dzi')
            os.remove(baseoutputdir+sys.argv[1]+'.dzi')
        except OSError:
            pass
 
        try:
            print('Removing '+baseoutputdir+sys.argv[1]+'_files')
            shutil.rmtree(baseoutputdir+sys.argv[1]+'_files')
        except OSError as e:
            print("Error: %s : %s" % (baseoutputdir+sys.argv[1]+'_files', e.strerror))
 
#        try:
#            print('Removing '+baseoutputdir+sys.argv[1]+'_files')
#            os.system('rm -fr "%s"' % baseoutputdir+sys.argv[1]+'_files')
#        except OSError:
#            pass

if dopngsave or dothumbnail:
        try:
            print('Removing '+baseoutputdir+sys.argv[1]+".png")
            os.remove(baseoutputdir+sys.argv[1]+".png")
        except OSError:
            pass

if dochopimage:
    try:
        print('Removing '+baseoutputdir+sys.argv[1]+'_chop.dzi')
        os.remove(outputname+'.dzi')
    except OSError:
        pass
    try:
        print('Removing '+baseoutputdir+sys.argv[1]+'_chop_files')
        os.system('rm -fr "%s"' % baseoutputdir+sys.argv[1]+'_chop_files')
    except OSError:
        pass

import glob

print(temporarydirectory+sys.argv[1])
rowimagefolder=sorted(glob.glob(temporarydirectory+sys.argv[1]+'/*.v'))

rowimagelist=[]
for f in rowimagefolder:
    print(f)
    rowimagelist.append(f)
    

images = [pyvips.Image.new_from_file(filename, access="sequential") for filename in rowimagelist]

TT=pyvips.Image.arrayjoin(images,across=1,background=[0]) # ,hspacing=0,vspacing=0,shim=0,halign="low"


if dothumbnail:

    TT.thumbnail_image(thumbnailsize).write_to_file(os.path.expanduser(baseoutputdir+sys.argv[1]+".png"))
    print(baseoutputdir+sys.argv[1]+".png thumbnail saved")
    

if dopngsave:

        TT.write_to_file(os.path.expanduser(baseoutputdir+sys.argv[1]+".png"))
        print(baseoutputdir+sys.argv[1]+".png saved")

if dodzsave:

        TT.dzsave(os.path.expanduser(baseoutputdir+sys.argv[1]),suffix='.jpg[Q=100]')
        print(baseoutputdir+sys.argv[1]+".dzi saved")

if dochopimage:
        
        TT.dzsave(os.path.expanduser(baseoutputdir+sys.argv[1])+'_chop',suffix='.jpg[Q=100]',depth='one',tile_size=choptilesize,overlap=0)
        print(baseoutputdir+sys.argv[1]+"_chop.dzi saved")
        
t1 = time.time()

print("")
print("")
print("{0} seconds".format(t1 - t0))
print("")
print("")

