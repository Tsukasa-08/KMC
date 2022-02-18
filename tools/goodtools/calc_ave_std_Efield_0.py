#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import argparse
import toml


# In[58]:


#parserを作成
parser = argparse.ArgumentParser(description="calculate mean displacement and standard deviation under NO Efield.                                                Output file is 'Efield_0_ave_std.toml' .")

parser.add_argument('-d', '--displacement_csv', required=True,                    help="csv file containing all diffusion atoms displacement.")

args = parser.parse_args()
#変位一覧が入力されたcsvファイルを読み込む
a = np.loadtxt(args.displacement_csv, delimiter=",", skiprows=1)


# In[59]:


#x,y,z軸方向の平均変位と標準偏差を求める
b = np.average(a, axis=0)
c = np.std(a, axis=0)


# In[67]:


#tomlファイルに平均変位と標準偏差を出力する
toml_string = """

['Efield_0']
    mean_displacement =  [{0}, {1}, {2}]
    stddev_displacement =  [{3}, {4}, {5}]
   
""".format(str(b[0]), str(b[1]), str(b[2]), str(c[0]), str(c[1]), str(c[2]),)

parsed_toml = toml.loads(toml_string)

with open('Efield_0_ave_std.toml', 'wt') as f:
    toml.dump(parsed_toml, f)

