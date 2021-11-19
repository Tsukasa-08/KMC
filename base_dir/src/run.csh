#!/bin/csh
#$ -cwd
#$ -q kod1
#$ -q kod2
#$ -V
#$ -N NAME
#$ -S /bin/zsh
#$ -pe smp 1

./test_kmc_test  >& log_cout
