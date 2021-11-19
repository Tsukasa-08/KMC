rm ./DiffusionCoefficient
rm ./ElectricalConductivity
rm ./log_cout
rm ./mean_displacement.csv
rm ./OUTPUT

sed -e "s/NAME/TEST/g" ./run.csh | qsub
