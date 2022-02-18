#!/bin/bash

awk -F, -v 'OFS=,' '{print $1,"2", $2, "1",$NF}' jmpdata*.csv > JMPDATA
sed -i /^\#/d JMPDATA
