awk -F"," '{print $5}' mean_displacement.csv | sed '/^$/d' | awk -f ~/awk_file/standard_deviation.awk
