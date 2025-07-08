#!/bin/bash
#SBATCH --job-name marc
#SBATCH -N1 --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --account=ELIX6_castrign
#SBATCH --partition=g100_usr_prod

set -e

module load profile/lifesc
module load autoload gromacs

directory="."
gas="co2"

usage() {
  echo "usage: $0 [-h] [-d <directory>] [-g <gas>]"
}

main() {
  while getopts ":hd:g:" opt; do
    case $opt in
      h) usage; exit 0 ;;
      d) directory="$OPTARG" ;;
      g) gas="${OPTARG,,}" ;;
      :) echo "error: missing argument for option '-$OPTARG'" >&2; usage; exit 1 ;;
      ?) echo "error: unknown option '-$OPTARG'" >&2; usage; exit 1 ;;
    esac
  done

  config_dir="$(realpath "$(cd "$(dirname "$0")" && pwd)/../config")"

  if [[ ! -d "$config_dir" ]]; then
    echo "error: '$config_dir' does not exist" >&2
    exit 1
  fi

  required_files=(
    minimization.mdp
    npt.mdp
    nvt.mdp
    production.mdp
  )

  for file in "${required_files[@]}"; do
    if [[ ! -f "$config_dir/$file" ]]; then
      echo "error: missing file $config_dir/$file" >&2
      exit 1
    fi
  done

  pushd "$(pwd)" > /dev/null
  cd "$directory"
  trap 'popd > /dev/null' EXIT

  if [[ ! -f "./topol.top" ]]; then
    echo "error: 'topol.top' not found in $directory" >&2
    exit 1
  fi

  if [[ ! -f "./ions.gro" ]]; then
    echo "error: 'ions.gro' not found in $directory" >&2
    exit 1
  fi

  if [[ ! -f "./$gas.gro" ]]; then
    echo "error: '$gas.gro' not found in $directory" >&2
    exit 1
  fi

  gmx_mpi grompp -f "$config_dir/minimization.mdp" -c ./ions.gro -p ./topol.top -o ./minimization.tpr
  srun gmx_mpi mdrun -s topol.tpr -ntomp 1 -v -maxh 24.0 -deffnm ./minimization

  gmx_mpi make_ndx -f "./$gas.gro" -o "./$gas.ndx" << EOF
0 & ! a H*
q
EOF
  echo "${gas^^}" | gmx_mpi genrestr -f "./$gas.gro" -n "./$gas.ndx" -o "./posre_$gas.itp" -fc 1000 1000 1000
  sed "/#include \"$gas.itp\"/a\\
\\
; Ligand position restraints\
\\
#ifdef POSRES\
\\
#include \"posre_$gas.itp\"\
\\
#endif\
\\
  " ./topol.top > ./topol.top.tmp && mv ./topol.top.tmp ./topol.top

  gmx_mpi make_ndx -f ./minimization.gro -o ./index.ndx << EOF
1 | 13
q
EOF

  gmx_mpi grompp -f "$config_dir/nvt.mdp" -c ./minimization.gro -r ./minimization.gro -p ./topol.top -n ./index.ndx -o ./nvt.tpr
  srun gmx_mpi mdrun -s topol.tpr -ntomp 1 -v -maxh 24.0 -deffnm ./nvt

  gmx_mpi grompp -f "$config_dir/npt.mdp" -c ./nvt.gro -t ./nvt.cpt -r ./nvt.gro -p ./topol.top -n ./index.ndx -o ./npt.tpr
  srun gmx_mpi mdrun -s topol.tpr -ntomp 1 -v -maxh 24.0 -deffnm ./npt

  gmx_mpi grompp -f "$config_dir/production.mdp" -c ./npt.gro -t ./npt.cpt -p ./topol.top -n ./index.ndx -o ./production.tpr
  srun gmx_mpi mdrun -s topol.tpr -ntomp 1 -v -maxh 24.0 -deffnm ./production
}

main "$@"
