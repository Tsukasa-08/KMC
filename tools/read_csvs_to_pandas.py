#!/usr/bin/env python
# coding: utf-8

# In[244]:


import numpy as np
import pandas as pd
import argparse
import toml
import glob
import re
import os
from scipy import constants


# In[245]:


#parserを作成
#parser = argparse.ArgumentParser(description="calculate D_tag : tracer-diffusion coefficient , and Sigma_tag : electronic conductivity .                                                 Output file is 'tracer_D_Sigma.toml' .")

#parser.add_argument('-d', '--displacement_csv', required=True, default='all_mean_displacement.csv',                    help="csv file containing all diffusion atoms displacement.")

#parser.add_argument('-p', '--parameter', required=True, default='OUTPUT',                    help="toml file containing temperture, total_time, concentration, ion_charge, and Efield_strength.")

#args = parser.parse_args()


# In[246]:


#calc[n]が配下にある場合

filepaths = glob.glob('./calc*/**/mean_*.csv', recursive=True)

#extract n : calc[n] for using n as MultiIndex : calc
list_of_calc_files = [int("".join(re.findall(r'\d+', filepath))) for filepath in filepaths]


# In[247]:


#print(filepaths)
#print(list_of_calc_files)


# In[248]:


data_list = []
for fileitem in filepaths:
    tmp = pd.read_csv(fileitem, index_col=[0,1], comment="#")
    data_list.append(tmp)

dataset = pd.concat(data_list, keys=list_of_calc_files, ignore_index= False)


# In[249]:


dataset.index.names = ["calc", "KMC_times", "diffusion_id"]


# In[251]:


parameter = toml.load("OUTPUT")


# In[250]:


dataset


# In[ ]:




