#!/bin/bash

d=$1
mode=$2

if [ -z "$3" ]
  then
    param_r=0.25
  else
  	param_r=$3
fi

if [ -z "$4" ]
  then
    param_mf=120
  else
  	param_mf=$4
fi

if [ -z "$5" ]
  then
    param_thr=1.5
  else
  	param_thr=$5
fi

echo "Running ims-bicluster for $d..."
echo "./bin/ims-bicluster --mat$mode --input=data/$d --eigens=20 --r=$param_r --maxforce=$param_mf --threshold=$param_thr"
./bin/ims-bicluster --mat$mode --input=data/$d --eigens=20 --r=$param_r --maxforce=$param_mf --threshold=$param_thr
suffix=$(date '+%y%m%d-%H%M%S')
echo "Making report for $d..."
reportname=$(bash gen-report.sh $d $mode)
echo "Report written to reports/$reportname.pdf"
exit
