#!/bin/bash
#SBATCH --partition=allq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=diff
#SBATCH --mem=2G
#SBATCH -o job.%N.%j.out
#SBATCH -e job.%N.%j.err

make
