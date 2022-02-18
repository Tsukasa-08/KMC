#!/bin/bash

cd "$(dirname "$0")"
cd ../..
./tools/goodtools/make_toml_Efield_0.sh
./tools/goodtools/cat_mean_displacement.sh Efield_0
./tools/goodtools/calc_tracer_D_Sigma.py -d cat_mean_displacement.csv -p Efield_0_parameter.toml


