#!/bin/bash
if [[ $# -le 1 ]] ; then
    echo 'Argument 1: infofile, Argument 2: outname. Exiting...'
    exit 0
fi
#cwd=$(pwd)
echo Creating row images...
infofile=$1
outname=$2
tmpdir=/backup/tmp/
#tmpdir=/var/TMP/
mkdir -p $tmpdir$outname
rm -f $tmpdir$outname/*.v
nrows=$(awk 'NR==2' $infofile)
echo $nrows rows
for (( c=1; c<=$nrows; c++ ))
do
   echo Row $c of $nrows
   python ~/Documents/affinetransform6r.py $infofile $outname $c &
   pids[${c}-1]=$!
done
for pid in ${pids[*]}; do
    wait $pid
done
python ~/Documents/affinetransform6f.py $outname
#cd $(echo ~/Data/deepzoom/"$outname"_files/14/)
#v=$(ls *.jpg | sort -V | tail -1 | cut -d'_' -f1)
#acrossnumber=$((v+1))
#vips arrayjoin "$(ls *.jpg | sort -t_ -k2g -k1g)" ~/Data/deepzoom/$outname.jpg[Q=99] --across $acrossnumber
#cd $(echo "$cwd")
