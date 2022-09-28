#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd
import numpy as np
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#import seaborn as sns 
from scipy import constants
idx = pd.IndexSlice


# In[3]:


#parserを作成
"""parser = argparse.ArgumentParser(description="concatnate all mean_displacement.csv under calc[n], with MultiIndex[calc, KMC_times, diffusion_id]. n is the number of calc directories.")

parser.add_argument('-p', '--search_dir_path', required=True, default='.',                    help="Directory path which you wanna search for mean_displacement.csv ")


args = parser.parse_args()

search_dir_path = args.search_dir_path
"""
search_dir_path = "."


# In[2]:


concat_1866K_df = pd.read_csv('./concat_modified_mean_displacement_1.csv', index_col=[i for i in range(0,8)])
#concat_mean_df = concat_mean_df.append(pd.read_csv('random_2to8_concat_modified_mean_displacement.csv', index_col=[i for i in range(0,7)]))
concat_1866K_df


# In[3]:


concat_output_df = pd.read_csv('./concat_modified_OUTPUT', index_col=[i for i in range(0,5)])
#concat_output_df = concat_output_df.append(pd.read_csv('Random_2to8_concat_modified_OUTPUT', index_col=[i for i in range(0,5)]))
concat_output_df


# In[4]:


concat_1866K_df = concat_1866K_df.join(concat_output_df) 


# In[ ]:


#concat_1866K_df = concat_1866K_df.rename(index=lambda s:int(s[:-1]), level='temperture').sort_index(axis=0)


# In[ ]:


#concat_1866K_df = concat_1866K_df.rename(index=lambda s:int(s[1:]), level='xY').sort_index(axis=0)


# In[7]:


#concat_1866K_df = concat_1866K_df.join(concat_output_df) 
#tempertureのKを除去(600K→600)
#concat_1866K_df = concat_1866K_df.rename(index=lambda s:int(s[:-1]), level='temperture').sort_index(axis=0)
#xYのYを消去(Y10→10)
#concat_1866K_df = concat_1866K_df.rename(index=lambda s:int(s[1:]), level='xY').sort_index(axis=0)
concat_1866K_df


# In[12]:


#Efield=0の平均変位を算出して、level"Efield"をdropすることで、config[n]内のEfield=0,1,2,3に対してbroadcastされてセンタリングできる
mean_drift_df = concat_1866K_df.groupby(level=["config",'xY','temperture','config_n','Efield']).mean().loc[idx[:,:,:,:,0],idx['dx':'dz']].reset_index(level='Efield', drop=True)
mean_drift_df


# In[13]:


#元データからdrift分を差し引いてセンタリングする、dx,dy,dz以外は0埋めして影響なし
concat_centering_df = concat_1866K_df.subtract(mean_drift_df,fill_value=0)
concat_centering_df


# In[14]:


del concat_1866K_df


# In[17]:


concat_centering_df['squared_dx'] = concat_centering_df['dx']**2
concat_centering_df['squared_dy'] = concat_centering_df['dy']**2
concat_centering_df['squared_dz'] = concat_centering_df['dz']**2
concat_centering_df.groupby(level=["config",'xY','temperture','config_n','Efield']).mean()


# In[18]:


target_Efield_0 = (concat_centering_df.index.get_level_values('Efield') == 0 )
#print(target_Efield_0)


# In[19]:


calc_Dtracer_df = concat_centering_df[target_Efield_0].groupby(level=['config','xY','temperture','config_n']).mean().loc[:,idx['total_time':'squared_dz']]
#calc_Dtracer_df


# In[20]:


calc_Dtracer_df["Dtracer_x"] = calc_Dtracer_df["squared_dx"] / (2*calc_Dtracer_df["total_time"]) * pow(10,-16)
calc_Dtracer_df["Dtracer_y"] = calc_Dtracer_df["squared_dy"] / (2*calc_Dtracer_df["total_time"]) * pow(10,-16)
calc_Dtracer_df["Dtracer_z"] = calc_Dtracer_df["squared_dz"] / (2*calc_Dtracer_df["total_time"]) * pow(10,-16)
#tempertureのKを除去(600K→600)
calc_Dtracer_df = calc_Dtracer_df.rename(index=lambda s:int(s[:-1]), level='temperture').sort_index(axis=0)
#xYのYを消去(Y10→10)
calc_Dtracer_df = calc_Dtracer_df.rename(index=lambda s:int(s[1:]), level='xY').sort_index(axis=0)

pd.options.display.precision = 4
#calc_Dtracer_df


# In[21]:


calc_Dtracer_df['Dtracer_average'] = calc_Dtracer_df.loc[:,idx['Dtracer_x':'Dtracer_z']].mean(axis=1)
#calc_Dtracer_df


# In[22]:

"""
target_Efield_123 = (concat_centering_df.index.get_level_values('Efield') != 0 )
#print(target_Efield_123)


# In[23]:


calc_Dsigma_df = concat_centering_df[target_Efield_123].groupby(level=['config','xY','temperture','config_n','Efield']).mean()
#tempertureのKを除去(600K→600)
calc_Dsigma_df = calc_Dsigma_df.rename(index=lambda s:int(s[:-1]), level='temperture').sort_index(axis=0)
#xYのYを消去(Y10→10)
calc_Dsigma_df = calc_Dsigma_df.rename(index=lambda s:int(s[1:]), level='xY').sort_index(axis=0)
calc_Dsigma_df['Dsigma'] = calc_Dsigma_df.loc[:,idx['dx':'dz']].max(axis=1) * constants.k * calc_Dsigma_df.index.get_level_values('temperture') / (constants.elementary_charge * calc_Dsigma_df['Efield_strength'] * calc_Dsigma_df['total_time']) * pow(10,-16)
#calc_Dsigma_df


# In[98]:



haven_x = calc_Dtracer_df['Dtracer_x'] / calc_Dsigma_df.loc[idx[:,:,:,:,1],:]['Dsigma']
#haven_y = calc_Dtracer_df['Dtracer_y'] / calc_Dsigma_df.loc[idx[:,:,:,:,2],:]['Dsigma']
#haven_z = calc_Dtracer_df['Dtracer_z'] / calc_Dsigma_df.loc[idx[:,:,:,:,3],:]['Dsigma']
#haven_x.groupby(level=['temperture','config_n']).plot()
#plt.legend()
haven_x_df = pd.DataFrame(haven_x)
haven_x_df.columns = ['haven_x']
haven_x_df['Dtracer'] = calc_Dtracer_df['Dtracer_x']
haven_x_df['Dsigma'] = calc_Dsigma_df.loc[idx[:,:,:,:,1],:]['Dsigma'].reset_index(level='Efield',drop=True)



# In[104]:


#sns.relplot(x='xY' ,y='haven_x' , row='temperture',data=haven_x_df)

"""
# In[1]:


calc_Dtracer_df.to_pickle('calc_Dtracer_df_Random.pickle')
#calc_Dsigma_df.to_pickle('calc_Dsigma_df_Random.pickle')
#haven_x_df.to_pickle('haven_x_df.pickle')

calc_Dtracer_df.to_csv('calc_Dtracer_df_Random.csv')
#calc_Dsigma_df.to_csv('calc_Dsigma_df_Random_1296atom.csv')
#haven_x_df.to_csv('haven_x_df.csv')

