# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 18:20:42 2023

@author: rkerrid

Step 4: Calculate protein level intensities


"""

from silac_dia_tools.pipeline.utils import dlfq_functions as dlfq
from silac_dia_tools.pipeline.report import protein_intensities_report
import os
import pandas as pd
import numpy as np


def create_protein_intensity_directory(path):
    # Combine the paths
    new_folder_path = os.path.join(path, 'protein intensities')
    
    # Create the new folder
    if not os.path.exists(new_folder_path):
        os.makedirs(new_folder_path)
        print(f"Folder preprocessing created successfully at {new_folder_path}")
    else:
        print(f"Folder preprocessing already exists at {new_folder_path}")
        
def output_href(path):
    create_protein_intensity_directory(path)
    
    # Import protein level ratios
    ratios = pd.read_csv(path + 'preprocessing/protein_ratios.csv')
    
    # Generate href df
    h_ref = ratios.groupby('Protein.Group')['H intensity'].median()
    h_ref = h_ref.reset_index().rename(columns={'H intensity': 'h_ref'})
    
    # merge href onto ratios
    merged_df_h = ratios.merge(h_ref, on='Protein.Group', how='inner')

    merged_df_h['H normalized total intensity'] = merged_df_h['h_ref'] / merged_df_h['H to stack ratio']
    merged_df_h['H intensity'] = merged_df_h['H normalized total intensity'] * merged_df_h['H to stack ratio']
    merged_df_h['M intensity'] = merged_df_h['H normalized total intensity'] * merged_df_h['M to stack ratio']
    merged_df_h['L intensity'] = merged_df_h['H normalized total intensity'] * merged_df_h['L to stack ratio']
        
    # Assign intensities to relevant columns
    merged_df_h['Total intensity'] =  merged_df_h['L intensity'] + merged_df_h['M intensity']
    merged_df_h['NSP intensity'] = merged_df_h['M intensity']
    
    # Generate subsetted dfs based on channels
    light_hnorm =  merged_df_h[['Run', 'Protein.Group', 'L intensity']]
    nsp_hnorm = merged_df_h[['Run', 'Protein.Group', 'NSP intensity']]
    total_hnorm = merged_df_h[['Run', 'Protein.Group', 'Total intensity']]
    reference = merged_df_h[['Run', 'Protein.Group', 'H intensity']]
    
    # Pivot tables to output format
    light_hnorm = light_hnorm.pivot(index='Protein.Group', columns='Run', values='L intensity')
    nsp_hnorm = nsp_hnorm.pivot(index='Protein.Group', columns='Run', values='NSP intensity')
    total_hnorm = total_hnorm.pivot(index='Protein.Group', columns='Run', values='Total intensity')
    reference = reference.pivot(index='Protein.Group', columns='Run', values='H intensity')
    
    #remove empty rows
    light_hnorm = remove_empyt_rows(light_hnorm)
    nsp_hnorm = remove_empyt_rows(nsp_hnorm)
    total_hnorm = remove_empyt_rows(total_hnorm)
    
    # nsp to light ratio and light to H ratio
    
    light_hnorm.to_csv(path + 'protein intensities/light_href.csv', sep=',')  
    nsp_hnorm.to_csv(path + 'protein intensities/nsp_href.csv', sep=',')
    total_hnorm.to_csv(path + 'protein intensities/total_href.csv', sep=',')
    reference.to_csv(path + 'protein intensities/reference_href.csv', sep=',')
    
    protein_intensities_report.create_intensities_report(total_hnorm, path)
    
    print('Saved H reference normalized protein intensities')
    return total_hnorm, nsp_hnorm, merged_df_h
        
def output_unnorm(path, contains_reference, pulse_channel='M'):
    create_protein_intensity_directory(path)
    # Import ratios table
    ratios = pd.read_csv(path + 'preprocessing/protein_ratios.csv')
    
    # Assign intensities to relevant columns
    nsp_channel = pulse_channel + ' intensity'
    ratios['Total intensity'] = ratios['L intensity'] +  ratios[nsp_channel]
    ratios['NSP intensity'] = ratios[nsp_channel]
    
    # Generate subsetted dfs based on channels
    light_unnorm = ratios[['Run', 'Protein.Group', 'L intensity']]
    nsp_unnorm = ratios[['Run', 'Protein.Group', 'NSP intensity']]
    total_unnorm = ratios[['Run', 'Protein.Group', 'Total intensity']]
    
    # Pivot tables to output format
    light_unnorm = light_unnorm.pivot(index='Protein.Group', columns='Run', values='L intensity')
    nsp_unnorm = nsp_unnorm.pivot(index='Protein.Group', columns='Run', values='NSP intensity')
    total_unnorm = total_unnorm.pivot(index='Protein.Group', columns='Run', values='Total intensity')
    
    #remove empty rows
    light_unnorm = remove_empyt_rows(light_unnorm)
    nsp_unnorm = remove_empyt_rows(nsp_unnorm)
    total_unnorm = remove_empyt_rows(total_unnorm)
    
    # nsp to light ratio and light to H ratio
    
    # Save tables to CSV
    light_unnorm.to_csv(path+'protein intensities/light_unnorm.csv', sep=',')
    nsp_unnorm.to_csv(path+'protein intensities/nsp_unnorm.csv', sep=',')
    total_unnorm.to_csv(path+'protein intensities/total_unnorm.csv', sep=',')
    if contains_reference:
        reference = ratios[['Run', 'Protein.Group', 'H intensity']]
        reference = reference.pivot(index='Protein.Group', columns='Run', values='H intensity')
        reference.to_csv(path+'protein intensities/reference_unnorm.csv', sep=',')
    print('Saved unnormalized protein intensities')

    return nsp_unnorm, total_unnorm

    
def output_dlfq(path, pulse_channel='M'):
    create_protein_intensity_directory(path)
    
    # Import ratios
    ratios = pd.read_csv(f'{path}preprocessing/protein_ratios.csv')
    nsp_ratio = pulse_channel + ' to stack ratio'
  
    # Read in precursors and save ms1 translated only file as input into directLFQ
    silac_precursors = pd.read_csv(path + 'preprocessing/silac_precursors.tsv', sep='\t')
    silac_precursors = silac_precursors[silac_precursors['quantity type']== 'Ms1.Translated']
    silac_precursors.to_csv(path+'preprocessing/silac_precursors_dlfq_in.tsv', sep='\t')
    dlfq.run_lfq( f'{path}preprocessing/silac_precursors_dlfq_in.tsv', 
              file =  f'{path}preprocessing/dlfq_protein_intensities.tsv',
              num_cores=1)
    
    # After directLFQ finishes running, read in results, format, and merge onto ratios df
    lfq_df = pd.read_csv(path+'preprocessing/dlfq_protein_intensities.tsv', sep='\t')
    lfq_df = lfq_df.melt(id_vars=['protein'], var_name = 'Run', value_name = 'Intensity')
    lfq_df.rename(columns = {'protein': 'Protein.Group'}, inplace=True)
    merged_df = ratios.merge(lfq_df, on=['Protein.Group','Run'], how='inner')
    
    # Assign intensities to relevant columns
    merged_df['L intensity'] = merged_df['L to stack ratio'] * merged_df['Intensity']
    merged_df['Total intensity'] = (merged_df['L to stack ratio'] * merged_df['Intensity']) + (merged_df[nsp_ratio] * merged_df['Intensity'])
    merged_df['NSP intensity'] = (merged_df[nsp_ratio] * merged_df['Intensity'])
    
    # Generate subsetted dfs based on channels
    total_lfq = merged_df[['Run', 'Protein.Group', 'Total intensity']]
    nsp_lfq = merged_df[['Run', 'Protein.Group', 'NSP intensity']]
    light_lfq = merged_df[['Run', 'Protein.Group', 'L intensity']]
    
    # Pivot tables to output format
    total_lfq = total_lfq.pivot(index='Protein.Group', columns='Run', values='Total intensity')
    nsp_lfq = nsp_lfq.pivot(index='Protein.Group', columns='Run', values='NSP intensity')
    light_lfq = light_lfq.pivot(index='Protein.Group', columns='Run', values='L intensity')
    
    #remove empty rows
    total_lfq = remove_empyt_rows(total_lfq)
    nsp_lfq = remove_empyt_rows(nsp_lfq)
    light_lfq = remove_empyt_rows(light_lfq)
    
    # nsp to light ratio and light to H ratio
    
    # Save tables to CSV
    total_lfq.to_csv(path+'protein intensities/total_dlfq.csv', sep=',')
    nsp_lfq.to_csv(path+'protein intensities/nsp_dlfq.csv', sep=',')
    light_lfq.to_csv(path+'protein intensities/light_dlfq.csv', sep=',')
    print('Saved directLFQ normalized protein intensities')
    
    protein_intensities_report.create_intensities_report(total_lfq, path)
    return total_lfq, nsp_lfq, light_lfq


def remove_empyt_rows(df):
    df.replace([np.nan, np.inf, -np.inf], 0, inplace=True) # replace infinity and nan values with 0
    cols = df.columns.values.tolist()[1:] # get a list of all sample columns
    df = df[df[cols].sum(axis=1) != 0] # remove all rows containing 0 values
    return df
