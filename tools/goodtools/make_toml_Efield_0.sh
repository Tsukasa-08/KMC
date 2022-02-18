#!/bin/bash

find . -name INPUT | head -n 1 | xargs -I{} cp {} .
find . -name OUTPUT | head -n 1 | xargs -I{} cp {} .
exec >> Efield_0_parameter.toml
#echo "[Efield_0]"
echo -n "temperture = " ; cat INPUT OUTPUT | grep TEMP | sed -e 's/[^0-9]//g'
echo -n "concentration = " ; cat INPUT OUTPUT | grep concentration | awk '{print $4}'
echo -n "total_time = " ; cat INPUT OUTPUT | grep total_time | awk '{print $4}'
echo -n "ion_charge = 1"
exec >> /dev/tty
