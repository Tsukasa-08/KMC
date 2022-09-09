rm -f DiffusionCoefficient
rm -f IonicConductivity
rm -f mean_displacement.csv
rm -f log_cout
rm -f OUTPUT
./test_kmc_test 1> log_cout 2>&1
