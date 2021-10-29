rm ./DiffusionCoefficient
rm ./ElectricalConductivity
rm ./log_cout
rm ./mean_displacement.csv
rm ./OUTPUT

make all

sed -e "s/NAME/TEST/g" ./run.csh | qsub
