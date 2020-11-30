#!/bin/zsh
#$ -S /bin/zsh
#$ -cwd
#$ -V
#$ -j y
#$ -N ace
#$ -o std.log
#$ -pe all_pe* 36
#============ Shell Script ============

mpirun -np 36 /storage/common_new/src/vasp.5.4.4/bin/vasp_gam
