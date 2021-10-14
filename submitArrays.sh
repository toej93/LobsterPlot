#!/bin/bash

#SBATCH --job-name=run_Arrays
#SBATCH --ntasks=2 --nodes=1 --mem-per-cpu 5g --partition scavenge --requeue
#SBATCH --mail-type=FAIL
#SBATCH --time=02:35:00
#SBATCH --output=logs/Arrays.out    # Standard output and error log

#source /gpfs/loomis/home.grace/jat89/.bash_cuoreMC
module load miniconda
conda activate ROOT6
python getArrays.py 0 &
python getArrays.py 1 &
wait

