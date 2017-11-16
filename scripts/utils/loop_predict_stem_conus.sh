#!/bin/bash

# Make new dir to store processed files if it doesn't already exist
mv_dir=$1"processed"
if [ ! -d "$mv_dir" ]; then
	mkdir $mv_dir
fi

# Loop through and run the script for each param file found
for file in $( find $1 -type f -name "predict_stem_params_*.txt" )
do
	python /vol/v2/stem/stem-git/scripts/predict_stem_conus.py $file
	mv "$file" "$mv_dir"
done
