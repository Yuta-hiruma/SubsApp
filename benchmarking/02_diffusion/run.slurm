#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --threads-per-core=1
#SBATCH --partition=allq
#SBATCH --nodes=1
#SBATCH -o job.%N.%j.out
#SBATCH -e job.%N.%j.err
#SBATCH --job-name jstmk
#SBATCH --mem-per-cpu=5G
#SBATCH --exclusive

make
