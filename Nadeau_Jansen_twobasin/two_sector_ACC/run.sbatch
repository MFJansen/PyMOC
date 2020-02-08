#!/bin/bash
#SBATCH --job-name=PyMOC
#SBATCH --output=pymoc_sbatch.out
#SBATCH --error=pymoc_sbatch.err
#SBATCH --time=96:00:00
#SBATCH --partition=mfj
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --exclusive
#SBATCH --mem-per-cpu=2000

module load python
python twobasin_NadeauJansen_kprofile_series.py 
