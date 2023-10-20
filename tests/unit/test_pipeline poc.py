# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 08:18:58 2023

@author: rkerrid
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pipeline import filtering_diann_report as fdia
from pipeline import format_silac_precursors as pdia
from pipeline import extract_protein_level_ratios as rdia
from pipeline import calculate_protein_intensities as idia
from icecream import ic


test_data_no_spikein = 'C:/data/silac_dia_tools_files/data/no spikein data/'
test_data_spikein = 'C:/data/silac_dia_tools_files/data/spikein data/'

# # #test filter
# df_si = fdia.import_and_filter(test_data_spikein, update=True)
# df_nosi = fdia.import_and_filter(test_data_no_spikein, update=True)

# # #test format precursors
# df_pre_si = pdia.format_silac_channels(test_data_spikein)
# df_pre_nosi = pdia.format_silac_channels(test_data_no_spikein)

# # #extract data
# rdia.calculate_protein_level_ratios(test_data_no_spikein)
# rdia.calculate_protein_level_ratios(test_data_spikein)

# ##calculate intensities
# # no spike in
# idia.output_dlfq(test_data_no_spikein)
# idia.output_unnorm(test_data_no_spikein, False)

# # spike in
# idia.output_href(test_data_spikein)
# idia.output_unnorm(test_data_spikein, True)

# nsp_href_df = pd.read_csv(test_data_spikein + 'protein intensities/nsp_href.csv', sep=',')
# total_href_df = pd.read_csv(test_data_spikein + 'protein intensities/total_href.csv', sep=',')

# androgen_df = total_href_df[total_href_df['Protein.Group']== 'O60885' ]#ER P11474 

# nsp_feratins = nsp_href_df[nsp_href_df['Protein.Group']== 'P02794' ] #nsp_href_df['Protein.Group']== 'P02792' 

# nsp_feratins = nsp_href_df[nsp_href_df['Protein.Group']== 'P02792' ]

# total_cdk = total_href_df[total_href_df['Protein.Group']== 'P11802' ]#CDK4 (CDK 6 not in any) Q00534

df = pd.read_csv(test_data_spikein + 'protein intensities/nsp_href.csv', sep=',')

def getog2_FC(df, col1, col2, fc_col):
    df[fc_col] = np.log2((df[col1] ) / (df[col2] ))
    return df

treatment = 'FAC'
ax_number = 10
col1 = 'Janice_20231005_AJW_HSdia_px07_Hspike_control_1'
col2 = f'Janice_20231005_AJW_HSdia_px07_Hspike_{treatment}_1'
fc_col = 'fc rep 1'
df = getog2_FC(df, col1, col2, fc_col)
col1 = 'Janice_20231005_AJW_HSdia_px07_Hspike_control_2'
col2 = f'Janice_20231005_AJW_HSdia_px07_Hspike_{treatment}_2'
fc_col = 'fc rep 2'
df = getog2_FC(df, col1, col2, fc_col)

plt.scatter(df['fc rep 1'], df['fc rep 2'])
plt.grid(True)
plt.ylim(-ax_number, ax_number)
plt.xlim(-ax_number, ax_number)


for i, row in df.iterrows():
    if row['fc rep 1'] > 2 or row['fc rep 1'] < -2 or row['fc rep 2'] > 2 or row['fc rep 2'] < -2:
        
        plt.annotate(row['Protein.Group'], (row['fc rep 1'], row['fc rep 2']), textcoords="offset points", xytext=(0,10), ha='center')
        print('done')
plt.show()