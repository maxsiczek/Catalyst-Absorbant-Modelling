#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --tasks-per-node=36
#SBATCH --cpus-per-task=1
#SBATCH --mem=200G
#SBATCH --job-name=job
###SBATCH --partition=standard
#SBATCH --partition=standard
###SBATCH --partition=ccei_biomass
#SBATCH --time=48:00:00
#SBATCH --export=NONE
#SBATCH --no-step-tmpdir
#SBATCH --mail-type=ALL
#SBATCH --mail-user=salmanak@udel.edu
export VALET_PATH=/work/ccei_biomass/sw/valet
vpkg_require vasp/5.4.1:VTST ase/3.16.2:python3
. /opt/shared/slurm/templates/libexec/openmpi.sh
. /work/ccei_biomass/sw/vasp/vasp-runtime.sh
export VASP_COMMAND="mpiexec vasp_std"
export VASP_PP_PATH=/work/ccei_biomass/sw/vasp/vasp_psp/v54


ulimit -c unlimited
export f77_dump_flag=TRUE
exec_vasp
mpi_rc=$?
exit $mpi_rc
