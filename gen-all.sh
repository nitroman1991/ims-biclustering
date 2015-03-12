#!/bin/bash

d=$1
if [ -z "$2" ]
  then
    param_r=0.25
  else
  	param_r=$2
fi

if [ -z "$3" ]
  then
    param_mf=120
  else
  	param_mf=$3
fi

if [ -z "$4" ]
  then
    param_thr=1.5
  else
  	param_thr=$4
fi

echo "Running ims-bicluster for $d..."
echo "./bin/ims-bicluster --input=data/$d --eigens=20 --r=$param_r --maxforce=$param_mf --threshold=$param_thr"
./bin/ims-bicluster --input=data/$d --eigens=20 --r=$param_r --maxforce=$param_mf --threshold=$param_thr
suffix=$(date '+%y%m%d-%H%M%S')
echo "Making report for $d..."
reportname=$(bash gen-report.sh $d)
echo "Report written to reports/$reportname.pdf"
exit