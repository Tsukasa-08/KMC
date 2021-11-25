#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#importまとめ
import pymatgen as mg
import pprint
import csv
import sys

#おおもとの全て含まれたPOSCAR(vasp形式)を読み込む
args = sys.argv
if len(args) == 2 : 
    structure = mg.Structure.from_file(args[1])
else:
    print("引数にはPOSCAR(vasp形式)を1つ入力してください")

#酸素イオンのindexをリスト化
O_sites_index = [i for i, site in enumerate(structure) if site.specie == mg.Element("O") ]

#どのindexからプロトンが始まるかをproton_start_numberに格納
for i, site in enumerate(structure) :
    if site.specie == mg.Element("H") : 
        proton_start_number = i
        break

#同時に専有できないプロトンサイトを1行ずつcsv形式で出力
blocking_list = []
for i in O_sites_index :
    list_proton_site = []
    #酸素イオンから距離2Å以下のプロトンサイトを抜き出し、
    #プロトンサイトのみのPOSCARにしたときのindexを記録する
    for proton_site in structure.get_neighbors(structure[i], 2):
        list_proton_site.append(proton_site[2]-proton_start_number)
    blocking_list.append(list_proton_site)

#csv形式のファイルを開き、出力する
with open("blocking_list.csv", "w") as f:
    writer = csv.writer(f)
    writer.writerows(blocking_list)



# In[116]:



    


# In[ ]:




