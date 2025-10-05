#!/bin/bash -l

#SBATCH --job-name='MR_PPP_forward'
#SBATCH --nodes=1
#SBATCH --time=23:00:00
#SBATCH --partition=cpu
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=8GB

source /.../bin/activate
cd /.../

snakemake -j 1000 --latency-wait=40 --keep-going --cluster "sbatch -p urblauna --job-name='Protein_Exposure' --time={resources.time} --nodes={resources.nodes} --cpus-per-task={resources.cpus} --mem={resources.mem}"