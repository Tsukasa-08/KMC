#!/bin/bash

find . -name INPUT | head -n 1 | xargs -I{} cp {} .
find . -name OUTPUT | head -n 1 | xargs -I{} cp {} .
exec >> Efield_0_parameter.toml
#echo "[Efield_0]"
echo -n "temperture = " ; cat INPUT OUTPUT | grep TEMP | head -n 1 | sed -e 's/[^0-9]//g'
echo -n "concentration = " ; cat INPUT OUTPUT | grep concentration | head -n 1| awk '{print $4}'
echo -n "total_time = " ; cat INPUT OUTPUT | grep total_time | head -n 1 | awk '{print $4}'
echo -n "ion_charge = 1"
exec >> /dev/tty
