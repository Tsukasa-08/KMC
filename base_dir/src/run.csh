#!/bin/csh
#$ -cwd
#$ -V
#$ -N NAME
#$ -S /bin/zsh
#$ -pe smp 40

rm DiffusionCoefficient
rm ElectoricalConductivity
rm log
./test_kmc_test  > dev/null
