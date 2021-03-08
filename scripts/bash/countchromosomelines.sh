#!/bin/bash
# request resources:
#PBS -l nodes=1:ppn=16
cd $PBS_O_WORKDIR

#probename=$(cat /home/ac14037/wheat/genome/probename.txt)
#echo $probename

#grep -A 1 $probename /home/ac14037/wheat/genome/35kprobes.fa > ./probe.extracted$probename.fa


LISTOFFILES=($(ls -l /home/ac14037/wheat/genome/fullchromosomesequences/ | awk '{print $9}' | awk -v ORS=" " 1 | sed 's/ //'))

for i in "${LISTOFFILES[@]}" 
do
	echo "the number of lines in '$i' is" >> chromosomecounts.txt
	wc -l /home/ac14037/wheat/genome/fullchromosomesequences/$i >> chromosomecounts.txt
done
