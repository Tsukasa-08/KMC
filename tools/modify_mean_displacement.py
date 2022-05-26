#!/usr/bin/env python
# coding: utf-8

# In[10]:


import argparse
import glob
import os


# In[11]:


#parserを作成
parser = argparse.ArgumentParser(description="modify all mean_displacement.csv under search directory to modified_mean_displacement.csv.")

parser.add_argument('-p', '--search_dir_path', required=True, default='.',                    help="Directory path which you wanna search for mean_displacement.csv ")


args = parser.parse_args()

search_dir_path = args.search_dir_path

#search_dir_path = "."
        
#search for all "mean_displacement.csv"
filepaths = glob.glob(os.path.join(search_dir_path, '**/mean_*.csv'), recursive=True)

#mean_displacement.csvに対して修正操作を行う
for filepath in filepaths :
    
    print(filepath)
    counter = 0
    with open(os.path.join(os.path.dirname(filepath),'modified_mean_displacement.csv'), mode='w') as t:

        t.write("KMC_times,diffusion_id,dx,dy,dz,counter\n")

        with open(os.path.join(os.path.dirname(filepath),'mean_displacement.csv')) as f:
            for c,line in enumerate(f):
                if (("#" in line) or (not line.replace('\n',''))) : 
                    continue
                else :
                    if ("KMC" in line) :
                        counter+=1
                    else : 
                        t.write(",".join([str(counter),line]))


# In[ ]:




