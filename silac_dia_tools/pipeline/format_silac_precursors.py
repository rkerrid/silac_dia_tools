# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 12:06:59 2023

@author: rkerrid

Step 2: Module for formatting report_filtered.tsv output from step one

"""
import pandas as pd
import numpy as np
import sys
sys.path.append('D:/Projects phd/General scripts for proteomics/SILAC DIA tools/')
from pipeline.report import precursor_report
from icecream import ic


#computing SILAC intensities for each precursor (and quantification type)
def format_silac_channels(path):
    new_path = f'{path}preprocessing/'
    df = import_filtered(new_path)
    ic(df)
    df = parse_data_for_channel_info(df) 
    ic(df)
    df = combine_modified_precursors(df)
    df = stack_intensities(df)
    df = drop_nan(df)
    print('Generating report...')
    precursor_report.silac_precursor_qc(df, path)
    print('Saving precursors...')
    df.to_csv(new_path+'silac_precursors.tsv', sep='\t')
    print('Done!')
    return df

##import/export file
def import_filtered(path, update=False):
    # Define chunk size (number of rows to load at a time)
    chunk_size = 100000
    chunks = []
    count = 1
    print("Importing report_filtered.tsv")
    # Iterate through the file in chunks and apply preprocessing functions
    for chunk in pd.read_table( f'{path}report_filtered.tsv',sep="\t", chunksize=chunk_size):
        chunks.append(chunk)
        if update:
            print('Importing chunk ', count)
            count+=1
        # if count == 2:
        #     break
    # Concatenate all chunks into a single DataFrame
    df = pd.concat(chunks, ignore_index=True)
    return df

def parse_data_for_channel_info(df):
    print('Parsing data for SILAC intensities')
    # Extracting SILAC labels from Precursor.ID and adding to separate column
    df['Label'] = df['Precursor.Id'].str.extract(r'\(SILAC-(K|R)-([HML])\)')[1]    
    df['Precursor'] = df['Stripped.Sequence'] + df['Precursor.Charge'].astype(str)
    
    # Splitting the data into two, one for 'Ms1.Translated' and another for 'Precursor.Translated'
    ms1_translated_df = df.copy()
    ms1_translated_df['intensity'] = df['Ms1.Translated']
    ms1_translated_df['quantity type'] = 'Ms1.Translated'

    precursor_translated_df = df.copy()
    precursor_translated_df['intensity'] = df['Precursor.Translated']
    precursor_translated_df['quantity type'] = 'Precursor.Translated'
    
    # Concatenate the two dataframes to create the parsed dataframe
    parsed_df = pd.concat([ms1_translated_df, precursor_translated_df], ignore_index=True)
    return parsed_df[['Run', 'Protein.Group', 'Precursor.Id', 'Precursor', 'Label', 'intensity', 'quantity type', 'Precursor.Quantity', 'Lib.PG.Q.Value']]


def combine_modified_precursors(df):
    print('combining modified precursors')
    # Define aggregation functions for columns
    agg_functions = {
        'intensity': 'first',
        'Lib.PG.Q.Value': 'first',
    }
    # Aggregate data using groupby and agg function
    agg_df = df.groupby(['Run', 'Protein.Group',  'Precursor', 'quantity type', 'Label']).agg(agg_functions).reset_index()

    # Pivot the table to get 'H', 'M', 'L' intensities in separate columns
    pivoted_df = agg_df.pivot_table(index=['Run', 'Protein.Group', 'Precursor', 'quantity type'], columns='Label', values='intensity', fill_value=0).reset_index()

    # Rename columns based on the label
    pivoted_df.columns.name = None
    pivoted_df = pivoted_df.rename(columns={'H': 'H intensity', 'M': 'M intensity', 'L': 'L intensity'})

    # Merge the pivoted_df with the original df to get other columns back
    df = pd.merge(agg_df.drop(columns=['Label', 'intensity']), pivoted_df, on=['Run', 'Protein.Group',  'Precursor', 'quantity type'])

    df = df.drop_duplicates()
    return df


#Stacking intensities and calculating intensity to stack ratio
def stack_intensities(df):
    columns_to_check = ['H intensity', 'M intensity', 'L intensity']
    for col in columns_to_check:
        if col not in df.columns:
            df[col] = 0
    df['Precursor.Quantity'] = df['H intensity'] + df['M intensity'] + df['L intensity']
    df['L to stack ratio'] = df['L intensity']/df['Precursor.Quantity']
    df['M to stack ratio'] = df['M intensity']/df['Precursor.Quantity']
    df['H to stack ratio'] = df['H intensity']/df['Precursor.Quantity']
    return df

def drop_nan(df):
    cols_to_check = [
    'Precursor.Quantity', 'H intensity', 'L intensity', 'M intensity', 
    'L to stack ratio', 'M to stack ratio', 'H to stack ratio'
    ]
    # Create a mask where all values in the specified columns are either 0 or NaN
    mask = df[cols_to_check].isin([0, np.nan]).all(axis=1)
    
    # Use the mask to drop rows from the DataFrame
    df = df[~mask]
    
    #rename precursor for directLFQ
    df = df.rename(columns={'Precursor':'Precursor.Id'})
    return df[['Run', 'Protein.Group', 'Precursor.Id','Precursor.Quantity','Lib.PG.Q.Value','quantity type','H intensity','M intensity','L intensity','H to stack ratio', 'M to stack ratio', 'L to stack ratio']]






    
