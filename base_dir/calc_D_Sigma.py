#!/usr/bin/env python
# coding: utf-8

# In[1]:


import configparser
from decimal import Decimal
from scipy import constants
from icecream import ic



#繰り返しで読み込んでいく
for_sigma_calc_ini = configparser.ConfigParser()
for_sigma_calc_ini.read("for_sigma_calc_sum.ini")

#Efield_0のmean_displacementはあとで使うのでリスト化しとく
mean_displacement_E0 = [""] * 3

#記録用のリスト作成
D_star = [""] * 3
Sigma_star = [""] * 3
Sigma = [[""] * 3 for i in range(3)]

for axis_n in range(4):
    
    #ic(axis_n)
    #データをiniファイルから読み込む
    read_data = for_sigma_calc_ini['Efield_{}'.format(str(axis_n))]
    
    
    #Efieldなしのとき
    if (axis_n == 0) :
        
        #平均変位を取得しておく
        temp = float(read_data.get('Temperture'))
        time = float(read_data.get('Total_time'))
        c = float(read_data.get('Concentration'))
        mean_x = float(read_data.get('Mean_x'))
        mean_y = float(read_data.get('Mean_y'))
        mean_z = float(read_data.get('Mean_z'))
        mean_displacement_E0[0] = mean_x
        mean_displacement_E0[1] = mean_y
        mean_displacement_E0[2] = mean_z
        
        #拡散係数を導出する
        mean_square_x = float(read_data.get('Mean_square_x'))
        mean_square_y = float(read_data.get('Mean_square_y'))
        mean_square_z = float(read_data.get('Mean_square_z'))
        mean_square_displacement = [mean_square_x, mean_square_y, mean_square_z]
        
        D = [""] * 3
        for i in range(3):
            D[i] = mean_square_displacement[i] / (2 * time)
            D_star[i] = D[i] * pow(10,-16)
        #    ic(i)
        #    ic(D_star[i])
        
        #伝導度を算出する
        for i in range(3):
            Sigma_star[i] = pow(constants.e, 2) * c * D[i] / (constants.k * temp) *pow(10,8)
        #    ic(i)
        #    ic(Sigma_star[i])
    
    #Efieldありのとき
    else : 
        temp = float(read_data.get('Temperture'))
        time = float(read_data.get('Total_time'))
        c = float(read_data.get('Concentration'))
        E = float(read_data.get('Efield'))
        mean_x = float(read_data.get('Mean_x'))
        mean_y = float(read_data.get('Mean_y'))
        mean_z = float(read_data.get('Mean_z'))
        mean_displacement = [mean_x, mean_y, mean_z]
        #Efield_0の変位を引いて、真の変位を出力しておく
        mean_displacement_real = [mean_displacement[i]-mean_displacement_E0[i] for i in range(3)]
        
        #Sigmaを計算する
        for i in range(3):
            Sigma[axis_n-1][i] = constants.e * c * mean_displacement_real[i]/ (E * time) *pow(10,8)
        
    
#確認用
"""
for i in range(3):
    ic(i)
    for j in range(3):
        ic(j)
        ic(Sigma[i][j])
"""
        


# In[17]:


#出力用strのlistを作成する
D_star_str = [str(i) for i in D_star]
Sigma_star_str = [str(i) for i in Sigma_star]
Sigma_str = [[str(j) for j in i] for i in Sigma]
#出力用の対応辞書を作成
xyz_dict = {0 : "x", 1 : "y", 2 : "z"}

#ファイルに出力する
with open("D_Sigma.ini", "w") as f:
    #
    #拡散係数を出力
    for i in range(3):
        axis_1 = xyz_dict[i]
        f.write('D_star_{}'.format(axis_1) + " = " + D_star_str[i] + '\n')
    #伝導度(電場なし)を出力
    for i in range(3):
        axis_1 = xyz_dict[i]
        f.write('Sigma_star_{}'.format(axis_1) + " = " + Sigma_star_str[i] + '\n')
    #伝導度(電場あり)を出力
    for i in range(3):
        axis_1 = xyz_dict[i]
        for j in range(3):
            axis_2 = xyz_dict[j]
            f.write('Sigma_{0}{1}'.format(axis_1, axis_2) + " = " + Sigma_str[i][j] + '\n')
            
        


# In[ ]:




