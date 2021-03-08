#!/bin/bash
# request resources:
#PBS -l nodes=1:ppn=16
cd $PBS_O_WORKDIR
parallel --joblog ./log_parallelordering.log "./p_ord_quad2.R {}" ::: ./chromosomes_preordering_quad2/*.csv
