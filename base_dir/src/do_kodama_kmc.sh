rm -f ./DiffusionCoefficient
rm -f ./IonicConductivity
rm -f ./log_cout
rm -f ./mean_displacement.csv
rm -f ./OUTPUT

sed -e "s/NAME/TEST/g" ./run.csh | qsub
