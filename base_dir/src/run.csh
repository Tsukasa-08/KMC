#!/bin/csh
#$ -cwd
#$ -V
#$ -N NAME
#$ -S /bin/zsh
#$ -pe smp 40

rm DiffusionCoefficient
rm ElectricalConductivity
rm log
rm sigma_log
./test_kmc_test  > dev/null
