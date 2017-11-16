#!/bin/bash
for file in $( find $1 -type f -name "confusion_params_*.txt" )
do
  python /vol/v2/stem/stem-git/scripts/confusion_matrix.py $file
  rm $file
done
