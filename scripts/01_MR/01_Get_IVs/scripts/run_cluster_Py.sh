#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=Add_relative_path # Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          			# Define number of nodes required
#SBATCH --time=0-00:05:00   		# Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=urblauna  		# Define the partition on which teh job shall run. May be omitted
#SBATCH --mem=4GB           		# Memory required per node

#SBATCH --output=/.../%x_%j.out

#################
#   JOB INFO    #
#################

python /.../01_relative_path.py

echo "Job done"
