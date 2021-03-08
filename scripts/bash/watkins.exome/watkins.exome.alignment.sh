#!/bin/bash

folder1=$1

cd $folder1
for i in ./*
do
	echo "unzipping $i"
	gunzip -k $i
done

read1=$(ls ./*_R1.fastq)
read2=$(ls ./*_R2.fastq)
read3=$(ls ./*_R0.fastq)



echo "performing alignment"

bwa mem -M -t 60 ~/project.phd.main/genome_assemblies/iwgsc/161010_Chinese_Spring_v1.0_pseudomolecules.fasta $read1 $read2 > alignment.sam

rm $read1
rm $read2
rm $read3

echo "converting to bam"

samtools view -@ 60 -bS alignment.sam > alignment.bam

echo "sorting bam"

samtools sort -O bam -o alignment.sorted.bam alignment.bam -@ 60

echo "indexing bam"

samtools index -c alignment.sorted.bam

echo "extracting unique reads"

samtools view -q 10 -b alignment.sorted.bam > alignment.uniq.sorted.bam -@ 60

echo "extracting only mapped reads"

samtools view -b -F 4 alignment.uniq.sorted.bam -@ 60 > alignment.uniq.nounmap.sort.bam

echo "sort by name"

samtools sort -n -o a5.bam alignment.uniq.nounmap.sort.bam -@ 60

echo "add mate coordinates (this is required for duplicate removal)"

samtools fixmate -m a5.bam a6.bam

echo "sort by coordinate again"

samtools sort -O bam -o a7.bam a6.bam -@ 60

echo "remove duplicates"

samtools markdup -r -s a7.bam a8.bam

echo "index again"

samtools index -c a8.bam

mkdir bams

mv a8.bam bams/
mv a8.bam.csi bams/

rm ./*.bam
rm ./*.sam

mkdir chromosomes

echo "calculate genotype likelihoods using mpileup"

parallel samtools mpileup -g -u -E -D --skip-indels -f ~/project.phd.main/genome_assemblies/iwgsc/161010_Chinese_Spring_v1.0_pseudomolecules.fasta -r {} bams/a8.bam -o chromosomes/{}.bcf ::: chr1A chr1B chr1D chr2A chr2B chr2D chr3A chr3B chr3D chr4A chr4B chr4D chr5A chr5B chr5D chr6A chr6B chr6D chr7A chr7B chr7D chrUn

cd chromosomes

"call SNPs / make vcf files"

ls | parallel ~/bin/bcftools call -v -c -N {} -o {.}.vcf

cat ./*.vcf > all.vcf

mv all.vcf ../
rm ./*.vcf
rm ./*.bcf

cd ..

bcftools view -i '%QUAL>=20' all.vcf > all2.vcf

rm all.vcf

cd bams

samtools depth a8.bam > a8.coverage

cat a8.coverage | awk '{ if ($3 > 19) { print } }' > a8.cov

rm a8.coverage