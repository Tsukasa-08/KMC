#!/bin/bash

#Efield_0ディレクトリから変位一覧をまとめたcsvファイルを作成する
#引数 : 対象となるEfield_0ディレクトリ

if [ $# != 1 ]; then
	echo 引数エラー: 対象となるEfield_0ディレクトリを引数に設定してください
	exit 1
fi

exec > cat_mean_displacement.csv

echo "#mean_x,mean_y,mean_z"

Efield_0_dir=$1

#mean_displacement.csvを読み込む
grep -h -r --include="mean*" '[0-9]' $1 |
	#平均変位以外の行を削除
	grep -v '^[^0-9]' |
	#平均変位x,y,zをcsv形式で出力
	awk -F"," 'BEGIN{OFS=","}{print $2,$3,$4}'

exec > /dev/tty
