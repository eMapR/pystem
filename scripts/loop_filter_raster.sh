#!/bin/bash
for file in $( find $1 -type f -name "*.txt" )
do
  python /vol/v2/stem/scripts/filter_raster.py $file
  rm $file
done
