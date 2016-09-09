#!/bin/bash
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -p sched_mit_hill,newnodes
#SBATCH --time=12:00:00
#SBATCH --exclusive
#SBATCH -J 4m-la@@SEED@@
#SBATCH -o slurm.%j.%N.out # STDOUT
#SBATCH -e slurm.%j.%N.err # STDERR

export CROSS_SECTIONS=/pool001/bforget/openmc/data/nndc/cross_sections.xml
export GFORTRAN_UNBUFFERED_ALL=y

../../../../bin/openmc-loo2-mapped