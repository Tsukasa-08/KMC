#!/usr/bin/env python
# coding: utf-8

# In[24]:


import numpy as np
import pprint
import configparser
import sys
import toml
from scipy import constants


# In[25]:


#変位一覧のcsvファイルを読み込む
cat_mean_csv = np.loadtxt("cat_mean_displacement.csv", delimiter=",", skiprows=1)

#計算に必要なパラメータを読み込む
Efield_0_parameter = toml.load("Efield_0_parameter.toml")


# In[26]:


#変位を全て2乗する
cat_mean_square_csv = np.power(cat_mean_csv,2)


# In[27]:


#x,y,z軸方向の平均二乗変位と標準偏差を求める
mean_square_displacement = np.average(cat_mean_square_csv, axis=0)
stddev_square_displacement = np.std(cat_mean_square_csv, axis=0)


# In[28]:


#パラメータを格納する
temperture = Efield_0_parameter["Efield_0"]["temperture"]
concentration = Efield_0_parameter["Efield_0"]["concentration"]
total_time = Efield_0_parameter["Efield_0"]["total_time"]
ion_charge = Efield_0_parameter["Efield_0"]["ion_charge"]


# In[29]:


#拡散係数を計算する
tracer_D_Ang = []  #[Ang.^2/s]
tracer_D_cm = []  #[cm^2/s]

for i in range(3):
    tracer_D_Ang.append(mean_square_displacement[i] / (2 * total_time))
    tracer_D_cm.append(tracer_D_Ang[i] * pow(10,-16))
    print(tracer_D_cm[i] )
    
#伝導度を計算する
tracer_Sigma = [] #[S/cm]

for i in range(3):
    tracer_Sigma.append(pow(ion_charge * constants.e , 2) * concentration * tracer_D_Ang[i] / (constants.k * temperture) * pow(10,8))
    print(tracer_Sigma[i])


# In[30]:


#一旦toml形式で出力する
toml_string = """

tracer_D = [{0[0]},{0[1]},{0[2]}]
tracer_Sigma = [{1[0]},{1[1]},{1[2]}]

""".format(tracer_D_cm, tracer_Sigma)

parsed_toml = toml.loads(toml_string)

with open('tracer_D_Sigma.toml', 'w') as f:
    toml.dump(parsed_toml, f)


# In[ ]:




