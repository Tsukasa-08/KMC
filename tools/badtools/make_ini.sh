#! /bin/bash

SCRIPT_DIR=$(cd $(dirname $0); pwd)
#Efield_i/calc_jを抽出
EFIELD_CALC_DIR=${SCRIPT_DIR##*K/}
#Efield_iを抽出
EFIELD_DIR=${EFIELD_CALC_DIR%%/*}
#for_sigma_calc_.iniを配置するディレクトリを指定(Efieldの手前、jmp*Kディレクトリに配置)
PUT_DIR=${SCRIPT_DIR%/E*}


exec > $PUT_DIR/for_sigma_calc_"$EFIELD_DIR".ini

echo "[$EFIELD_DIR]"

echo -n "Concentration = " ; grep concentration OUTPUT | uniq | awk '{print $4}'  
echo  "#単位 : [/Å^3]"

if [[ "$EFIELD_DIR" != *Efield_0* ]]  ; then
	echo -n "Efield = " ; grep Efield_strength OUTPUT | head -n 1 | awk '{print $3}'
	echo  "#単位 : [V/Å]"
fi


echo -n "Total_time = " ; grep total OUTPUT | awk '{print $4}' | awk -f ~/awk_file/standard_deviation.awk | awk -F"," '{print $1}' 
echo  "#単位 : [s]"

echo -n "Temperture = " ; grep TEMP INPUT | sed -e "s/[^0-9]//g"
echo  "#単位 : [K]"

echo -n "Mean_x = " ; grep mean_displacement_x OUTPUT | awk '{print $3}' | awk -f ~/awk_file/standard_deviation.awk | awk -F"," '{print $1}' 
echo  "#単位 : [Å]"

echo -n "Mean_y = " ; grep mean_displacement_y OUTPUT | awk '{print $3}' | awk -f ~/awk_file/standard_deviation.awk | awk -F"," '{print $1}' 
echo  "#単位 : [Å]"

echo -n "Mean_z = " ; grep mean_displacement_z OUTPUT | awk '{print $3}' | awk -f ~/awk_file/standard_deviation.awk | awk -F"," '{print $1}' 
echo  "#単位 : [Å]"

echo -n "Mean_square_x = " ; cat mean_displacement.csv | grep -v "[^0-9]$" | sed '/^$/d' | awk -F"," '{print $2}' | awk -f ~/awk_file/mean_square.awk
echo  "#単位 : [Å^2]"

echo -n "Mean_square_y = " ; cat mean_displacement.csv | grep -v "[^0-9]$" | sed '/^$/d' | awk -F"," '{print $3}' | awk -f ~/awk_file/mean_square.awk
echo  "#単位 : [Å^2]"

echo -n "Mean_square_z = " ; cat mean_displacement.csv | grep -v "[^0-9]$" | sed '/^$/d' | awk -F"," '{print $4}' | awk -f ~/awk_file/mean_square.awk
echo  "#単位 : [Å^2]"
