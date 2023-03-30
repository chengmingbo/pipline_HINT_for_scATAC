from __future__ import print_function

import os

samples = ['Sham_Female', 'Sham_Male', 'Day3_Female', 'Day3_Male', 'Day10_Female', 'Day10_Male']
for sample in samples:
    job_name = "split_{}".format(sample)
    command = "sbatch -J " + job_name + " -o " + "./cluster_out/" + job_name + "_out.txt -e " + \
              "./cluster_err/" + job_name + "_err.txt -t 120:00:00 --mem 180G ./run.zsh"
    os.system(command + " " + sample)
