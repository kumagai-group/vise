#!/bin/zsh
#$ -S /bin/zsh
#$ -cwd
#$ -V
#$ -j y
#$ -N test
#$ -o std.log
#$ -pe all_pe* 36
#$ -l h_rt=10:00:00
#============ Shell Script ============
# A script for sequential calculations.
#   Written by Yu Kumagai
#   Ver. 0.3.1  Time: Mon May  9 19:32:20 JST 2016
python ~/my_bin/vise/vise/main.py vr -vc "/home/kuma/bin/openmpi-1.8.1-intel16.0.2/bin/mpirun -np $NSLOTS /home/kuma/bin/vasp.5.4.4/bin/vasp_std"
