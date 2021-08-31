#!/bin/csh
#$ -cwd
#$ -V
#$ -N NAME
#$ -S /bin/zsh
#$ -pe smp 40

./test_kmc_test > /dev/null
