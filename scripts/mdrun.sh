#!/bin/bash
#SBATCH --job-name marc
#SBATCH -N1 --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --account=ELIX6_castrign
#SBATCH --partition=g100_usr_prod

module load profile/lifesc
module load autoload gromacs
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

gmx=$(which gmx_mpi)
config_dir="$(realpath "$(cd "$(dirname "$0")" && pwd)/../config")"

gmx grompp -f "$config_dir/minimization.mdp" -c ./ions.gro -p ./topol.top -o ./minimization.tpr
srun "$gmx" mdrun -s topol.tpr -ntomp 1 -v -maxh 24.0 -deffnm ./minimization

gmx make_ndx -f ./co2.gro -o ./co2.ndx << EOF
0 & ! a H*
q
EOF
echo "CO2" | gmx genrestr -f ./co2.gro -n ./co2.ndx -o ./posre_co2.itp -fc 1000 1000 1000
sed "/#include \"co2.itp\"/a\\
\\
; Ligand position restraints\
\\
#ifdef POSRES\
\\
#include \"posre_co2.itp\"\
\\
#endif\
\\
" ./topol.top > ./topol.top.tmp && mv ./topol.top.tmp ./topol.top

gmx make_ndx -f ./minimization.gro -o ./index.ndx << EOF
1 | 13
q
EOF

gmx grompp -f "$config_dir/nvt.mdp" -c ./minimization.gro -r ./minimization.gro -p ./topol.top -n ./index.ndx -o ./nvt.tpr
srun "$gmx" mdrun -s topol.tpr -ntomp 1 -v -maxh 24.0 -deffnm ./nvt

gmx grompp -f "$config_dir/npt.mdp" -c ./nvt.gro -t ./nvt.cpt -r ./nvt.gro -p ./topol.top -n ./index.ndx -o ./npt.tpr
srun "$gmx" mdrun -s topol.tpr -ntomp 1 -v -maxh 24.0 -deffnm ./npt

gmx grompp -f "$config_dir/production.mdp" -c ./npt.gro -t ./npt.cpt -p ./topol.top -n ./index.ndx -o ./production.tpr
srun "$gmx" mdrun -s topol.tpr -ntomp 1 -v -maxh 24.0 -deffnm ./production
