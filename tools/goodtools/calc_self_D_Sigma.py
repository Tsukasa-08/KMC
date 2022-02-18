#!/usr/bin/env python
# coding: utf-8

# In[20]:


import numpy as np
import pprint
import configparser
import argparse
import toml
from scipy import constants


# In[21]:



#parserを作成
parser = argparse.ArgumentParser(description="calculate D_j : self-diffusion coefficient , and Sigma_j : electronic conductivity .                                                 Output file is 'self_D_Sigma.toml' .")

parser.add_argument('-d', '--displacement_csv', required=True, default='cat_mean_displacement.csv',                    help="csv file containing all diffusion atoms displacement.")

parser.add_argument('-p', '--parameter', required=True, default='OUTPUT',                    help="toml file containing temperture, total_time, concentration, ion_charge, and Efield_strength.")

args = parser.parse_args()


# In[22]:


#変位一覧のcsvファイルを読み込む
cat_mean_csv = np.loadtxt(args.displacement_csv, delimiter=",", skiprows=1)

#計算に必要なパラメータを読み込む
Efield_0_parameter = toml.load(args.parameter)


# In[26]:


#x,y,z方向の変位を全て合計し全変位を算出する
sum_mean_csv = cat_mean_csv.sum(axis=0)


# In[31]:


#x,y,z軸方向の合計変位を2乗して拡散種数で割る
mean_square_displacement = np.power(sum_mean_csv, 2) / cat_mean_csv.shape[0]


# In[32]:


#パラメータを格納する
temperture = Efield_0_parameter["temperture"]
concentration = Efield_0_parameter["concentration"]
total_time = Efield_0_parameter["total_time"]
ion_charge = Efield_0_parameter["ion_charge"]


# In[33]:


#拡散係数を計算する
self_D_Ang = []  #[Ang.^2/s]
self_D_cm = []  #[cm^2/s]

for i in range(3):
    self_D_Ang.append(mean_square_displacement[i] / (2 * total_time))
    self_D_cm.append(self_D_Ang[i] * pow(10,-16))
    
#伝導度を計算する
self_Sigma = [] #[S/cm]

for i in range(3):
    self_Sigma.append(pow(ion_charge * constants.e , 2) * concentration * self_D_Ang[i] / (constants.k * temperture) * pow(10,8))


# In[34]:


#一旦toml形式で出力する
toml_string = """

self_D = [{0[0]},{0[1]},{0[2]}]
self_Sigma = [{1[0]},{1[1]},{1[2]}]

""".format(self_D_cm, self_Sigma)

parsed_toml = toml.loads(toml_string)

with open('self_D_Sigma.toml', 'w') as f:
    toml.dump(parsed_toml, f)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




