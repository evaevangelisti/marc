#!/bin/bash
#SBATCH --job-name marc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:4
#SBATCH --time=01:00:00
#SBATCH --partition=boost_usr_prod

export OMP_NUM_THREADS=8

export GMX_GPU_DD_COMMS=true
export GMX_GPU_PME_PP_COMMS=true
export GMX_FORCE_UPDATE_DEFAULT_GPU=true

module purge
module load profile/lifesc gromacs

set -e

GAS="co2"
MOLECULES=1

main() {
  source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../config" && pwd)/config.sh"

  if [ ! -f "$GASES_DIR/$GAS/$GAS.pdb" ] || [ ! -f "$GASES_DIR/$GAS/$GAS.itp" ]; then
    echo "error: gas '$GAS' not found in $GASES_DIR/" >&2; exit 1
  fi

  gas_pdb="$GASES_DIR/$GAS.pdb"
  gas_itp="$GASES_DIR/$GAS.itp"

  if [[ ! "$MOLECULES" =~ ^[0-9]+$ || "$MOLECULES" -lt 1 ]]; then
    echo "error: MOLECULES must be a positive integer" >&2; exit 1
  fi

  gmx editconf -f "$gas_pdb" -o "./$GAS.gro"

  cp "$gas_itp" "./$GAS.itp"
  sed "/#include \"oplsaa.ff\/forcefield.itp\"/a\\
\\
; Include topology for $GAS\
\\
#include \"$GAS.itp\"\
\\
" ./topol.top > ./topol.top.tmp && mv ./topol.top.tmp ./topol.top
  printf "%s%-$((21 - ${#GAS} - ${#MOLECULES}))s%d\n" "$GAS" "" "$MOLECULES" >> ./topol.top

  gmx insert-molecules -f ./processed.gro -ci "./$GAS.gro" -o ./complex.gro -nmol "$MOLECULES"

  gmx editconf -f ./complex.gro -o ./box.gro -bt dodecahedron -d 1.0

  gmx solvate -cp ./box.gro -cs tip4p.gro -p ./topol.top -o ./solvated.gro

  gmx grompp -f "$MDP_IONS" -c ./solvated.gro -p ./topol.top -o ./ions.tpr
  echo SOL | gmx genion -s ./ions.tpr -o ./ions.gro -p ./topol.top -pname NA -nname CL -neutral -conc 0.100

  gmx grompp -f "$MDP_MINIMIZATION" -c ./ions.gro -p ./topol.top -o ./minimization.tpr
  gmx mdrun -v -s ./minimization.tpr -deffnm ./minimization -ntmpi 4 -ntomp 8 -pin on -pinstride 1 -nb gpu -pme gpu -npme 1 -pmefft gpu -update gpu

  gmx make_ndx -f "./$GAS.gro" -o "./$GAS.ndx" << EOF
0 & ! a H*
q
EOF
  echo "${GAS^^}" | gmx genrestr -f "./$GAS.gro" -n "./$GAS.ndx" -o "./posre_$GAS.itp" -fc 1000 1000 1000
  sed "/#include \"$GAS.itp\"/a\\
\\
; Ligand position restraints\
\\
#ifdef POSRES\
\\
#include \"posre_$GAS.itp\"\
\\
#endif\
\\
  " ./topol.top > ./topol.top.tmp && mv ./topol.top.tmp ./topol.top

  gmx make_ndx -f ./minimization.gro -o ./index.ndx << EOF
1 | 13
q
EOF

  gmx grompp -f "$MDP_NPT" -c ./minimization.gro -r ./minimization.gro -p ./topol.top -n ./index.ndx -o ./nvt.tpr
  gmx mdrun -v -s ./nvt.tpr -deffnm ./nvt -ntmpi 4 -ntomp 8 -pin on -pinstride 1 -nb gpu -pme gpu -npme 1 -pmefft gpu -update gpu

  gmx grompp -f "$MDP_NVT" -c ./nvt.gro -t ./nvt.cpt -r ./nvt.gro -p ./topol.top -n ./index.ndx -o ./npt.tpr
  gmx mdrun -v -s ./npt.tpr -deffnm ./npt -ntmpi 4 -ntomp 8 -pin on -pinstride 1 -nb gpu -pme gpu -npme 1 -pmefft gpu -update gpu

  gmx grompp -f "$MDP_PRODUCTION" -c ./npt.gro -t ./npt.cpt -p ./topol.top -n ./index.ndx -o ./production.tpr
}

main "$@"
