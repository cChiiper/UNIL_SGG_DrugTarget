#!/bin/bash

#SBATCH --job-name='Data_extract'
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --partition urblauna
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=500MB

snakemake -j 1000 --latency-wait=10 --cluster "sbatch -p urblauna --job-name='Protein_Exposure' --time={resources.time} --nodes={resources.nodes} --cpus-per-task={resources.cpus} --mem={resources.mem}"

