#!/bin/csh
#$ -cwd
#$ -V
#$ -N NAME
#$ -S /bin/zsh
#$ -pe smp 40

rm DiffusionCoefficient
rm ElectoricalConductivity
./test_kmc_test  > dev/null
