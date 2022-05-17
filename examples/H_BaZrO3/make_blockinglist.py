#!/usr/bin/env python
# coding: utf-8

# In[2]:


#importまとめ
import pymatgen as mg
import csv
import argparse


# In[6]:


#おおもとの全て含まれたPOSCAR(vasp形式)を読み込む
parser = argparse.ArgumentParser(description='find blocking groups. print by csv format.')

parser.add_argument('-f', '--poscar', type=str, required=True,                         help='Structure file name in POSCAR format containing all site of                                atom including diffusion species.')

args = parser.parse_args()

structure = mg.Structure.from_file(args.poscar)


# In[ ]:


#酸素イオンのindexをリスト化
O_sites_index = [i for i, site in enumerate(structure) if site.specie == mg.Element("O") ]

#どのindexからプロトンが始まるかをproton_start_numberに格納
for i, site in enumerate(structure) :
    if site.specie == mg.Element("H") : 
        proton_start_number = i
        break


# In[ ]:


#同時に専有できないプロトンサイトを1行ずつcsv形式で出力
blocking_list = []
for i in O_sites_index :
    list_proton_site = []
    #酸素イオンから距離2Å以下のプロトンサイトを抜き出し、
    #プロトンサイトのみのPOSCARにしたときのindexを記録する、indexは1から始まることに注意
    for proton_site in structure.get_neighbors(structure[i], 2):
        list_proton_site.append(proton_site[2]-proton_start_number+1)
    blocking_list.append(list_proton_site)


# In[ ]:


#csv形式のファイルを開き、出力する
with open("blocking_list.csv", "w") as f:
    writer = csv.writer(f)
    writer.writerows(blocking_list)

print("indexは1から始まっていることに注意してください")

