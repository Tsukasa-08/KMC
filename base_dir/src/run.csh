#!/bin/csh
#$ -cwd
#$ -V
#$ -N NAME
#$ -S /bin/zsh
#$ -pe smp 1

./test_kmc_test  >& log_cout
