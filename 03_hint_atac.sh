#!/bin/bash
#
#### Job name
#SBATCH -J 30
#SBATCH -e logs/30.txt
#SBATCH -o logs/30.txt
#SBATCH -t 520:00:00
#SBATCH --mem=180G --cpus-per-task=10


## 00 please change the path to your own path
## 00 please change the organism to your own organism


PDIR=`pwd`

bam_loc=$PDIR/DiffFootprinting/BAM
peaks_loc=$PDIR/DiffFootprinting/Peaks
footprint_loc=$PDIR/DiffFootprinting/Footprints
motifmatching_loc=$PDIR/DiffFootprinting/MotifMatching
diff_loc=$PDIR/DiffFootprinting/Diff
bigwig_loc=$PDIR/DiffFootprinting/Diff


mkdir -p ${peaks_loc}
mkdir -p ${footprint_loc}
mkdir -p ${bigwig_loc}
mkdir -p ${motifmatching_loc}
mkdir -p ${bam_loc}
mkdir -p ${diff_loc}

echo ${footprint_loc}


## 01 merge all bam files that belongs to the same cluster
## C1~C6 stands for the six clusters
samtools merge -f --threads 10 ${bam_loc}/C1.bam  ../BAM/*_C1.bam
samtools merge -f --threads 10 ${bam_loc}/C2.bam  ../BAM/*_C2.bam
samtools merge -f --threads 10 ${bam_loc}/C3.bam  ../BAM/*_C3.bam
samtools merge -f --threads 10 ${bam_loc}/C4.bam  ../BAM/*_C4.bam
samtools merge -f --threads 10 ${bam_loc}/C5.bam  ../BAM/*_C5.bam
samtools merge -f --threads 10 ${bam_loc}/C6.bam  ../BAM/*_C6.bam

samtools index  ../BAM/*_C1.bam
samtools index  ../BAM/*_C2.bam
samtools index  ../BAM/*_C3.bam
samtools index  ../BAM/*_C4.bam
samtools index  ../BAM/*_C5.bam
samtools index  ../BAM/*_C6.bam





## 02 call peaks for each cluster & perform footprinting
footprinting(){
    	filename=$1
    	samtools index $1
    	macs2 callpeak -t $filename -n ${filename%.bam} --outdir $2 -g mm --nomodel -f BAMPE -q 0.01 --keep-dup all
    	sed -i '/random\|chrUn\|chrM/d' $2/${filename%.bam}_peaks.narrowPeak
    	sed -i '/random\|chrUn\|chrM/d' $2/${filename%.bam}_summits.bed
	rgt-hint footprinting --organism=mm10 --atac-seq --paired-end --output-location=$4 --output-prefix=${filename%.bam} $1 $2/${filename%.bam}_peaks.narrowPeak
	rgt-motifanalysis matching --organism=mm10 --output-location=$5 --input-files $4/${filename%.bam}.bed
}

cd ${bam_loc}
for filename in *.bam;
do
	footprinting $filename ${peaks_loc} ${bigwig_loc} ${footprint_loc} ${motifmatching_loc} &
done

wait


## 03 calculate the differential footprinting
cd ${diff_loc}
rgt-hint differential --window-size=1000 --organism=mm10 --bc --nc 64 \
--mpbs-files=${motifmatching_loc}/C1_mpbs.bed,\
${motifmatching_loc}/C2_mpbs.bed,\
${motifmatching_loc}/C3_mpbs.bed,\
${motifmatching_loc}/C4_mpbs.bed,\
${motifmatching_loc}/C5_mpbs.bed,\
${motifmatching_loc}/C6_mpbs.bed \
--reads-files=${bam_loc}/C1.bam,\
${bam_loc}/C2.bam,\
${bam_loc}/C3.bam,\
${bam_loc}/C4.bam,\
${bam_loc}/C5.bam,\
${bam_loc}/C6.bam \
--conditions=C1,C2,C3,C4,C5,C6 \
--output-location=${diff_loc} \
--output-prefix=All

