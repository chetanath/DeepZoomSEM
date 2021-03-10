#!/bin/sh                                                                                             

# $1 = directory                                                                                      

CURRENTDIR=$PWD

cd $1

imgext=$(find . -name "*.hdr" -print -quit | head -1 | cut -d'-' -f2- | cut -d"." -f1)

echo image extension = $imgext

find -regex '\./[0-9_]+' -type d -exec mv -n -- {}/Snap_1.$imgext {}.$imgext 2>/dev/null \;
find -regex '\./[0-9_]+' -type d -exec mv -n -- {}/Snap_1-$imgext.hdr {}-$imgext.hdr 2>/dev/null \; -empty -delete

find . -type d -empty -delete

cd $CURRENTDIR
