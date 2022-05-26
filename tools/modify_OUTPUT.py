#!/usr/bin/env python
# coding: utf-8

# In[10]:


import argparse
import glob
import os


# In[33]:


#parserを作成
parser = argparse.ArgumentParser(description="modify all OUTPUT under search directory to modified_OUTPUT.")

parser.add_argument('-p', '--search_dir_path', required=True, default='.',                    help="Directory path which you wanna search for OUTPUT ")


args = parser.parse_args()

search_dir_path = args.search_dir_path
#search_dir_path = "."
        
#search for all "mean_displacement.csv"
filepaths = glob.glob(os.path.join(search_dir_path, '**/OUTPUT'), recursive=True)

#mean_displacement.csvに対して修正操作を行う
for filepath in filepaths :
    
    print(filepath)
    counter = 0
    with open(os.path.join(os.path.dirname(filepath),'modified_OUTPUT'), mode='w') as t:

        with open(os.path.join(os.path.dirname(filepath),'OUTPUT')) as f:
            for c,line in enumerate(f):
                if ("KMC" in line):
                    counter+=1
                    if (counter >= 2) :
                        break
                else : 
                 if ("total_time" in line):
                     split_list = line.split()
                     del split_list[-1]
                     split_list.remove('t')
                     modified_str = " ".join(split_list)
                     t.write(modified_str+'\n')
                 if ("concentration" in line):
                     split_list = line.split()
                     del split_list[-1]
                     split_list.remove('c')
                     modified_str = " ".join(split_list)
                     t.write(modified_str+'\n')
                 if (("Efield_strength" in line) and (not "cm" in line)):
                     split_list = line.split()
                     del split_list[-1]
                     split_list[0] = split_list[0] + ".Ang"
                     modified_str = " ".join(split_list)
                     t.write(modified_str+'\n')

# In[ ]:




