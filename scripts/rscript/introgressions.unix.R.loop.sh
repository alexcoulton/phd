#!/bin/bash

array=( Synthetic_39 Synthetic_46 Savannah Opata )
for i  in "${array[@]}"
do
	ls /home/ac14037/project.phd.main/rotation1scripts_v4/original_data/genotype.data/820k.genotypes/ordered.split/ | parallel /home/ac14037/project.phd.main/rotation1scripts_v4/scripts/rscript/introgressions.unix.R $i
done



