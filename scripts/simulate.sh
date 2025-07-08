#!/bin/bash
#SBATCH --job-name marc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:4
#SBATCH --time=24:00:00
#SBATCH --partition=boost_usr_prod

export OMP_NUM_THREADS=8

export GMX_GPU_DD_COMMS=true
export GMX_GPU_PME_PP_COMMS=true
export GMX_FORCE_UPDATE_DEFAULT_GPU=true

set -e

module purge
module load profile/lifesc gromacs

main() {
  if [[ -f "md_prod1.cpt" ]]; then
    gmx mdrun -v -s ./production.tpr -deffnm ./md_prod2 -cpi ./md_prod1.cpt -ntmpi 4 -ntomp 8 -pin on -pinstride 1 -nb gpu -pme gpu -npme 1 -pmefft gpu -update gpu
  else
    gmx mdrun -v -s ./production.tpr -deffnm ./md_prod1 -ntmpi 4 -ntomp 8 -pin on -pinstride 1 -nb gpu -pme gpu -npme 1 -pmefft gpu -update gpu
  fi
}

main "$@"
