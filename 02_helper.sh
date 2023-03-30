#!/bin/sh
### Job name
#SBATCH -J D10_2
#SBATCH -e logs/D10_2.txt
#SBATCH -o logs/D10_2.txt

### Time your job needs to execute, e. g. 15 min 30 sec
#SBATCH -t 100:00:00

### Memory your job needs per node, e. g. 1 GB
#SBATCH --mem=200G



source ~/miniconda3/bin/activate
conda activate Seurat3


aaaaa="Healthy_Female"
bbbbb="prameter-b"
ccccc="prameter-c"
ddddd="prameter-d"
eeeee="prameter-e"


cluster_file=./data/clusters/${aaaaa}.txt
bam_file=../data/bamfiles/${aaaaa}/outs/possorted_bam.bam
output_location=../BAM

mkdir ../BAM

echo "python split_bam_by_cluster.py $cluster_file $bam_file $output_location $aaaaa"

python split_bam_by_cluster.py $cluster_file $bam_file $output_location $aaaaa
