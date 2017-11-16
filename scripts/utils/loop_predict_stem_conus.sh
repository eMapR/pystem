#!/bin/bash
for file in $( find $1 -type f -name "predict_stem_params_*.txt" )
do
  python /vol/v2/stem/stem-git/scripts/predict_stem_conus.py $file
  rm $file
done
