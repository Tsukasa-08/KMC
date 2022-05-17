#!/bin/bash

#MASTEQ形式のjmpdata.csvをtakahashiのKMC形式のJMPDATAに変換する
#拡散機構が空孔機構で、拡散種がプロトン(と空孔)のみの場合に有効
#takahashiのKMC内では空孔を1, プロトンを2に割り当てている
rm JMPDATA
awk -F, -v 'OFS=,' '{print $1,"2", $2, "1",$NF}' jmpdata*.csv > JMPDATA
sed -i /^\#/d JMPDATA
