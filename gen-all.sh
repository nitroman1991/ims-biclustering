#!/bin/bash

d=$1
mode=$2

if [ -z "$3" ]
  then
    param_r=$3
  else
  	param_r=0.2
fi

if [ -z "$4" ]
  then
    param_mf=$4
  else
  	param_mf=100
fi

echo "Running ims-bicluster for $d..."
echo "./bin/ims-bicluster --mat$mode --input=data/$d --eigens=20 --r=$param_r --maxforce=$param_mf"
./bin/ims-bicluster --mat$mode --input=data/$d --eigens=20 --r=$param_r --maxforce=$param_mf
suffix=$(date '+%y%m%d-%H%M%S')
echo "Making report for $d..."
reportname=$(bash gen-report.sh $d $mode)
echo "Report written to reports/$reportname.pdf"
exit
