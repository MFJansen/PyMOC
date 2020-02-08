#!/bin/bash
#SBATCH --job-name=adiabatic_ref
#SBATCH --output=1DAMOC.out
#SBATCH --error=1DAMOC.err
#SBATCH --time=2:00:00
#SBATCH --partition=mfj
#SBATCH --nodes=1
#SBATCH --exclusive

module load python

python run_JansenNadeau_2018.py  --adiabatic --pickup='save_adiabatic.npz' --diagfile='diags_adiabatic_ref' 

