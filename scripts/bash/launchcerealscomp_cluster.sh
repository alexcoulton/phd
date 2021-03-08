#!/bin/bash
# request resources:
#PBS -l nodes=1:ppn=16
cd $PBS_O_WORKDIR
parallel --joblog ./log_cerealscomp.log "./cerealscomp_cluster.R {}" ::: *.param
