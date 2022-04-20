#!/usr/bin/env python
# coding: utf-8

# In[28]:


import numpy as np
import pandas as pd
import argparse
import glob
import re
import os
import shutil
from scipy import constants

#split path into a list : [dirname, subdirname, ... ,and filename] 
def recursive_path_split(path_list, path):
    if os.path.abspath(path) == os.getcwd() : 
        return 0

    path_list.insert(0,os.path.split(path)[-1])
    recursive_path_split(path_list, os.path.split(path)[0])
    


# In[29]:


#parserを作成
"""parser = argparse.ArgumentParser(description="concatnate all mean_displacement.csv under calc[n], with MultiIndex[calc, KMC_times, diffusion_id]. n is the number of calc directories.")

parser.add_argument('-p', '--search_dir_path', required=True, default='.',                    help="Directory path which you wanna search for mean_displacement.csv ")


args = parser.parse_args()

search_dir_path = args.search_dir_path
"""
search_dir_path = "."


# In[30]:


#search paths to calc[n] directory
calc_paths = glob.glob(os.path.join(search_dir_path,'**/calc1'), recursive=True)


for calc_path in calc_paths:
    calc_dir = os.path.dirname(calc_path)
    
    #以下concat_mean_displacement_under_dir.pyを転記
    
    #calc[n]が配下にある場合

    #search for all "mean_displacement.csv"
    filepaths = glob.glob(os.path.join(calc_dir, 'calc*/**/mean_*.csv'), recursive=True)
    
    #search for a "OUTPUT"
    #output_path = glob.glob(os.path.join(calc_path, '**', "OUTPUT"), recursive=True)[0]
    output_path = glob.glob(os.path.join(calc_path, '**', "OUTPUT"), recursive=True)[0]
    
    data_list = []
    
    
    #extract n : calc"$n" for using n as MultiIndex,"calc"
    #list_of_calc_files = [int("".join(re.findall(r'\d+', filepath))) for filepath in filepaths]
    list_of_calc_files = []

    for fileitem in filepaths:
        path_split_list = []
        recursive_path_split(path_split_list, fileitem)
        calc_str = [l for l in path_split_list if 'calc' in l][0]
        list_of_calc_files.append(int(re.sub(r'\D',"", calc_str)))
        #mcsp_str = [l for l in path_split_list if 'mcsp' in l][0]
        #naverage_str = [l for l in path_split_list if 'average' in l][0]
    
        tmp = pd.read_csv(fileitem, index_col=[0,1], comment="#")
        data_list.append(tmp)

    #concatnate data_list, adding a Index "calc"
    if len(data_list) == 0:
        print("found no mean_displacement.csv")
    else:
        dataset = pd.concat(data_list, keys=list_of_calc_files, ignore_index= False)
    
        dataset.index.names = ["calc", "KMC_times", "diffusion_id"]
        
        #export and copy
    
        if (dataset.to_csv(os.path.join(calc_dir, "concat_mean_displacement.csv"))) :
            print("file saved at ", os.path.join(calc_dir, "concat_mean_displacement.csv"))
            
        if (shutil.copy(output_path, calc_dir)):
            print("file saved at ", os.path.join(calc_dir))
         
        


# In[79]:


ex_file_path = './mcsp_pow_10_1/n_average_pow_10_1/calc1/base_dir/src/mean_displacement.csv'


# In[65]:


os.path.exists(ex_file_path)


# In[80]:


type(os.path.split(ex_file_path))


# In[81]:


path_split_list = []
def recursive_path_split(path_list, path):
    if os.path.abspath(path) == os.getcwd() : 
        return 0
    
    path_list.insert(0,os.path.split(path)[-1])
    recursive_path_split(path_list, os.path.split(path)[0])
    print("if do", datetime.datetime.now())

#    else:
#        print("else do", datetime.datetime.now())
recursive_path_split(path_split_list, ex_file_path)


# In[82]:


print(path_split_list)


# In[85]:


calc_str = [l for l in path_split_list if 'calc' in l][0]
mcsp_str = [l for l in path_split_list if 'mcsp' in l][0]
naverage_str = [l for l in path_split_list if 'average' in l][0]


# In[60]:


os.path.abspath(ex_file_path)


# In[23]:


1 + 1


# In[24]:


print(_)


# In[25]:


1 == 1


# In[26]:


print(_)


# In[27]:


1 == 1
if _ :
    print("it is True!")


# In[33]:


print("what")
print(1 + 1)
print(_)


# In[34]:


1 + 1


# In[35]:


print(_)


# In[36]:


1 + 1


# In[37]:


print("hello")


# In[38]:


print(_)


# In[ ]:




