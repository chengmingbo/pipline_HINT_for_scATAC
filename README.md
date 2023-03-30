# pipline_HINT_for_scATAC

## 0. put `sbatch_param` to your local execution path (e.g. `$HOME/.local/bin/`)
  sbatch_param enables you to set the parameters for sbatch command
  e.g. `sbatch_param -f yourscript.sh -n yourjobname -a your1stparam -b your2ndparam -c your3rdparam`
  where `$aaaaa` is the first parameter, `$bbbbb` is the second parameter, and `$ccccc` is the third parameter in your script

  you can also run without `sbatch_param` directly using `sbatch 02_helper.sh` by setting all the parameters in the script

## 1. run script to split barcodes by clusters and samples
   `Rscript 01_getCluster.R`

   this will generate barcodes text file recording the cluster information in each sample named file.


## 2. run script to split bamfiles by the split barcode produced in step 1
   `02_all_run.sh` will submit all jobs to slurm nodes by passing paramters to `02_helper.sh`,
   where `02_helper.sh` call `split_bam_by_cluster.py` to split bam files by clusters and samples

   This will generate a lot of bam files, "`$sampleName`\_`$clusterName`.bam" in the `BAM` folder

## 3. run HINT-ATAC on bam files
  `sh 03_hint_atac.sh` includes a bunch of commands to run HINT-ATAC on bam files
    01. merge all bam files that belong to the same cluster
    02. call peaks for each cluster & perform footprinting
    03. calculate the differential footprinting

