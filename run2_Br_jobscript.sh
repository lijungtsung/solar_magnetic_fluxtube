#!/bin/bash
#SBATCH --account=XXXXXX
#SBATCH --job-name="Br_network"
#SBATCH --output="out.log"
#SBATCH --nodes=5                # node count
#SBATCH --ntasks=200             # total number of tasks
#SBATCH --cpus-per-task=1        # cpu-cores per task
#SBATCH --mem-per-cpu=4G         # memory per cpu-core
#SBATCH --export=ALL
#SBATCH --time=3:00:00
#SBATCH --no-requeue
#SBATCH --mail-type=ALL


module load miniconda3/4.10.3-py37
module load openmpi/4.0.3
source activate fast-mpi4py

srun -m NoPack python -m mpi4py.futures run2_Br.py & wait