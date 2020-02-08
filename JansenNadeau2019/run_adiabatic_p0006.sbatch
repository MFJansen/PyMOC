#!/bin/bash
#SBATCH --job-name=adiabatic_p0006
#SBATCH --output=1DAMOC.out
#SBATCH --error=1DAMOC.err
#SBATCH --time=2:00:00
#SBATCH --partition=mfj
#SBATCH --nodes=1
#SBATCH --exclusive

module load python

python run_JansenNadeau_2018.py  --adiabatic --db=0.0006 --pickup='save_adiabatic.npz' --diagfile='diags_adiabatic_p0006'

