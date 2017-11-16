#!/bin/bash
for file in $( find $1 -type f -name "predict_rf_params_*.txt" )
do
  python /vol/v2/stem/scripts/predict_rf.py $file
  rm $file
done
