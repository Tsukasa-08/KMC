#!/usr/bin/env python
# coding: utf-8

# In[28]:


import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns 
from scipy import constants
idx = pd.IndexSlice


# In[29]:


#parserを作成
"""parser = argparse.ArgumentParser(description="concatnate all mean_displacement.csv under calc[n], with MultiIndex[calc, KMC_times, diffusion_id]. n is the number of calc directories.")

parser.add_argument('-p', '--search_dir_path', required=True, default='.',                    help="Directory path which you wanna search for mean_displacement.csv ")


args = parser.parse_args()

search_dir_path = args.search_dir_path
"""
search_dir_path = "."


# In[30]:


concat_1866K_df = pd.read_csv('./concat_modified_mean_displacement_1.csv', index_col=[i for i in range(0,8)])
#concat_mean_df = concat_mean_df.append(pd.read_csv('random_2to8_concat_modified_mean_displacement.csv', index_col=[i for i in range(0,7)]))


# In[31]:


concat_output_df = pd.read_csv('./concat_modified_OUTPUT', index_col=[i for i in range(0,5)])
#concat_output_df = concat_output_df.append(pd.read_csv('Random_2to8_concat_modified_OUTPUT', index_col=[i for i in range(0,5)]))


# In[32]:


concat_1866K_df = concat_1866K_df.join(concat_output_df) 


# In[33]:




# In[34]:


dx_sum_df = concat_1866K_df.loc[:,idx['dx':'dz']].groupby(level=['config','xY','temperture','config_n','Efield','calc']).sum()


# In[35]:


dx_sum_df


# In[36]:


mean_drift_df = dx_sum_df.groupby(level=['config','xY','temperture','config_n','Efield']).mean()


# In[37]:



#元データからdrift分を差し引いてセンタリングする、dx,dy,dz以外は0埋めして影響なし
concat_centering_df = dx_sum_df.subtract(mean_drift_df,fill_value=0)
concat_centering_df


# In[54]:


del concat_1866K_df


# In[55]:


concat_centering_df.std(ddof=0)


# In[56]:


concat_centering_df = concat_centering_df.join(concat_output_df)


# In[57]:


concat_centering_df['squared_dx'] = concat_centering_df['dx']**2
concat_centering_df['squared_dy'] = concat_centering_df['dy']**2
concat_centering_df['squared_dz'] = concat_centering_df['dz']**2


# In[58]:


target_Efield_0 = (concat_centering_df.index.get_level_values('Efield') == 0 )
#print(target_Efield_0)


# In[59]:


#calc_Dcollective_all_df = concat_centering_df[target_Efield_0].groupby(level=['config','xY','temperture','config_n']).mean().loc[:,idx['total_time':'squared_dz']]
#tempertureのKを除去(600K→600)
calc_Dcollective_all_df = concat_centering_df[target_Efield_0].rename(index=lambda s:int(s[:-1]), level='temperture').sort_index(axis=0)
#xYのYを消去(Y10→10)
calc_Dcollective_all_df = calc_Dcollective_all_df.rename(index=lambda s:int(s[1:]), level='xY').sort_index(axis=0)

#Dcollectiveを計算する、最後に拡散粒子数で割るのを忘れない
calc_Dcollective_all_df["Dcollective_x"] = calc_Dcollective_all_df["squared_dx"] / (2*calc_Dcollective_all_df["total_time"] * 27 * calc_Dcollective_all_df.index.get_level_values('xY')) * pow(10,-16)
calc_Dcollective_all_df["Dcollective_y"] = calc_Dcollective_all_df["squared_dy"] / (2*calc_Dcollective_all_df["total_time"] * 27 * calc_Dcollective_all_df.index.get_level_values('xY')) * pow(10,-16)
calc_Dcollective_all_df["Dcollective_z"] = calc_Dcollective_all_df["squared_dz"] / (2*calc_Dcollective_all_df["total_time"] * 27 * calc_Dcollective_all_df.index.get_level_values('xY')) * pow(10,-16)


# In[60]:


calc_Dcollective_all_df['Dcollective_average'] = calc_Dcollective_all_df.loc[:,idx['Dcollective_x':'Dcollective_z']].mean(axis=1)
calc_Dcollective_all_df


# In[61]:


analysis_df = calc_Dcollective_all_df.loc[:,'Dcollective_x':'Dcollective_average']


# In[62]:


calc_Dcollective_df = analysis_df.groupby(level=['config','xY','temperture','config_n','Efield']).mean()


# In[63]:


calc_Dcollective_all_df.to_pickle('calc_Dcollective_df_by_confign_and_calc_Random.pickle')
#calc_Dsigma_ave_df.to_pickle('calc_Dsigma_ave_df_Random.pickle')
calc_Dcollective_all_df.to_csv('calc_Dcollective_df_by_confign_and_calc_Random.csv')
#calc_Dsigma_ave_df.to_csv('calc_Dsigma_ave_df_Random_1296atom.csv')


# In[64]:


calc_Dcollective_df.to_pickle('calc_Dcollective_df_by_confign_Random.pickle')
#calc_Dsigma_ave_df.to_pickle('calc_Dsigma_ave_df_Random.pickle')
calc_Dcollective_df.to_csv('calc_Dcollective_df_by_confign_Random.csv')
#calc_Dsigma_ave_df.to_csv('calc_Dsigma_ave_df_Random_1296atom.csv')


# In[ ]:




