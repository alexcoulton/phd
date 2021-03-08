#!/bin/bash

#create array of values from 1 to 50, 1 step
vals=($(seq 1 1 50))

parallel -j50 ~/project.phd.main/rotation1scripts_v4/scripts/rscript/recomb.simulation2.rscript.R ::: ${vals[@]}
