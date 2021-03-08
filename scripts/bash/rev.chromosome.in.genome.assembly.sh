#!/bin/bash
#Returns a genome assembly with the reverse complement of a particular chromosome
genome=$1
ext=".oneline.fa"
echo $genome$ext

~/bin/conv.2.oneline.fasta.sh $genome $genome$ext

sed -n 2p $genome$exti
#this script is incomplete
