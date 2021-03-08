#!/bin/bash
# request resources:
#PBS -l nodes=1:ppn=16
#cd $PBS_O_WORKDIR

#probename=$(cat /home/ac14037/wheat/genome/probename.txt)
#echo $probename

#grep -A 1 $probename /home/ac14037/wheat/genome/35kprobes.fa > ./probe.extracted$probename.fa


LISTOFFILES=($(ls -l ./fasta_files_chromosome_probes/ | awk '{print $9}' | awk -v ORS=" " 1 | sed 's/ //'))

for i in "${LISTOFFILES[@]}" 
do
	echo "performing blast on $i"
	echo blastn -db wheat_genome_blast_db -query ./fasta_files_chromosome_probes/$i -outfmt 6 -out blast-files-LG-vs-genome/$i.blast -num_threads 16
done
