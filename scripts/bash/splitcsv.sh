#!/bin/bash

FILE1=$1

tail -n +2 $FILE1 | split -l 1000000 - split_
for file in split_*
do
    head -n 1 $FILE1 > tmp_file
    cat $file >> tmp_file
    mv -f tmp_file "$file".csv
    rm $file
done