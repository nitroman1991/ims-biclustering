#!/bin/bash

d=$1
suffix=$(date '+%y%m%d-%H%M%S')
if [ -z "$2" ]
  then
    M=10
  else
  	M=$2
fi
if [ -z "$3" ]
  then
    N=10
  else
  	N=$3
fi
python python/visualize.py -x data/$d -y $mattype -e data/$d.val.csv -v data/$d.vec.csv -c data/$d.coords.csv -N $N -M $M -t $d
cd reports/latex && pdflatex tmp >/dev/null && pdflatex tmp >/dev/null && mv tmp.pdf ../$d-$suffix.pdf && cd ../..
echo $d-$suffix
exit