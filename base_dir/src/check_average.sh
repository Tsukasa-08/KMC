awk -F"," '{print $5}' mean_displacement.csv | head -n -1 | sed -n '3,$p' | awk -f ~/awk_file/standard_deviation.awk
