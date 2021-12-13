#! /bin/bash

#4つのiniファイルを1つにまとめる
for i in for_sigma_calc_E* ; do cat $i ; echo "" ; done > for_sigma_calc_sum.ini
