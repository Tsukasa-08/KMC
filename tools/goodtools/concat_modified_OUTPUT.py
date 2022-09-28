#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import argparse
import glob
import re
import os
import shutil
import toml
from scipy import constants

#split path into a list : [dirname, subdirname, ... ,and filename] 
def recursive_path_split(path_list, path):
    if os.path.abspath(path) == os.getcwd() : 
        return 0

    path_list.insert(0,os.path.split(path)[-1])
    recursive_path_split(path_list, os.path.split(path)[0])
    


# In[39]:


#parserを作成
parser = argparse.ArgumentParser(description="concat all modified_OUTPUT under search directory.")

parser.add_argument('-p', '--search_dir_path', required=True, default='.',                    help="Directory path which you wanna search for modified_OUTPUT. 1866K or Random is recommended. ")


args = parser.parse_args()

search_dir_path = args.search_dir_path


#search_dir_path = "."
#search for all "mean_displacement.csv"
filepaths = glob.glob(os.path.join(search_dir_path, '**/modified_OUTPUT'), recursive=True)


data_list = []
#list_of_calc_files = [int("".join(re.findall(r'\d+', filepath))) for filepath in filepaths]
list_of_calc_files = []

#modified_mean_displacement.csvに対して修正操作を行う
for filepath in filepaths :
    print(filepath)
    
    if ('1866K' in filepath or 'Random' in filepath): 
        #with open(os.path.join(os.path.dirname(filepath),'concat_modified_mean_displacement.csv'), mode='w') as t:

            #with open(os.path.join(os.path.dirname(filepath),'modified_mean_displacement.csv')) as f:


        path_split_list = []
        recursive_path_split(path_split_list, filepath)

        config_str = [l for l in path_split_list if '1866K' in l or 'Random' in l][0]
        #calc_str = [l for l in path_split_list if 'calc' in l][0]
        Efield_str = [l for l in path_split_list if 'Efield' in l][0]
        Efield_str = re.sub(r'\D',"",Efield_str)
        #list_of_calc_files.append(int(re.sub(r'\D',"", Efield_str)))
        Y_concentration_str = [l for l in path_split_list if 'Y' in l][0]

        calc_jmp_dirname_str = [l for l in path_split_list if 'jmp_' in l][0]

        #config_n_strを求める、"4"など
        config_n_str = re.findall('_([0-9]*)_', calc_jmp_dirname_str)[0]

        #temperture_strを求める、右から_を探してそれより右側を抽出、"600K"など
        target = '_'
        idx = calc_jmp_dirname_str.rfind(target)
        temperture_str = calc_jmp_dirname_str[idx+len(target):]

        tmp = toml.load(filepath)
        if ('Efield_strength' in tmp.keys()) : 
            tmp_df = pd.DataFrame({'total_time':tmp['total_time'], 'concentration' : tmp['concentration'], 'Efield_strength':tmp['Efield_strength']['Ang']} , index=[Efield_str])
        else : 
            tmp_df = pd.DataFrame({'total_time':tmp['total_time'], 'concentration' : tmp['concentration']} , index=[Efield_str])
        
        tmp_df.index.names = ['Efield']
        tmp_df = pd.concat([tmp_df],keys=[config_n_str],names=["config_n"])
        tmp_df = pd.concat([tmp_df],keys=[temperture_str],names=["temperture"])
        tmp_df = pd.concat([tmp_df],keys=[Y_concentration_str],names=["xY"])
        tmp_df = pd.concat([tmp_df],keys=[config_str],names=["config"])
        data_list.append(tmp_df)


# In[41]:


#concatnate data_list, adding a Index "calc"
if len(data_list) == 0:
    print("found no mean_displacement.csv")
else:
    dataset = pd.concat(data_list, ignore_index= False)

    if (dataset.to_csv(os.path.join(search_dir_path, "concat_modified_OUTPUT"))) :
        print("file saved at ", os.path.join(search_dir_path, "concat_modified_OUTPUT"))

