#!/bin/bash
#SBATCH --job-name=ref
#SBATCH --output=1DAMOC.out
#SBATCH --error=1DAMOC.err
#SBATCH --time=2:00:00
#SBATCH --partition=mfj
#SBATCH --nodes=1
#SBATCH --exclusive

module load python

python run_JansenNadeau_2018.py  --pickup='save_ref.npz'  --diagfile='diags_ref' 

