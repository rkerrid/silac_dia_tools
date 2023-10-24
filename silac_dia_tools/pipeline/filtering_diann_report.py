# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 11:43:50 2023

@author: rkerrid

Step 1: Module for filtering report.tsv output (DIA-NN version 1.8.1) with
SILAC settings as described in the README.md

Note: This script filters for contaminants by looking for the 'cont_' substring
in Protein.Groups so make sure your report.tsv is annotated in the same way or 
edit the remove_contaminants() funciton.

"""
import pandas as pd
import os
import json
import operator
from silac_dia_tools.pipeline.report import filtering_report
from pkg_resources import resource_filename


# Defining the relative path to configs directory 
CONFIG_DIR = os.path.join(os.path.dirname(__file__), '..', 'configs')
# json_path = 'filtering_parameters.json'

def import_and_filter(path, update=False):
    # Define chunk size (number of rows to load at a time)
    chunk_size = 100000
    chunks = []
    contams = []
    filtered_set = []
    
    # Load filtering parameters from JSON
    print('Loading filtering parameters')
    json_path = os.path.join(CONFIG_DIR, 'filtering_parameters.json')
    with open(json_path, 'r') as f:
        params = json.load(f)
        
    # Iterate through the file in chunks and apply preprocessing functions
    print('Beggining filtering')
    count = 1
    with open(f"{path}report.tsv", 'r', encoding='utf-8') as file:
        for chunk in pd.read_table(file,sep="\t", chunksize=chunk_size):
          
            # Apply filtering to each chunk
            chunk, contam = remove_contaminants(chunk)
            chunk, filtered_out = apply_filters(chunk, params)
            chunk = drop_cols(chunk) 
            
            # Append chunks from respective filtering steps
            filtered_set.append(filtered_out)
            contams.append(contam)
            chunks.append(chunk)
            # Update progress (optional)
            if update:
                print('chunk ', count,' processed')
            count+=1
        
    # Concatenate all chunks into a DataFrames
    filtered_set = pd.concat(filtered_set, ignore_index=True)
    contams = pd.concat(contams, ignore_index=True)
    df = pd.concat(chunks, ignore_index=True)
    
    # Pass filtering information to reports
    print('Generating filtering report')
    filtering_report.filtering_qc(df, contams, filtered_set, path, params)
    create_preprocessing_directory(path)
    print('Saving filtered_report.tsv')
    df.to_csv(f'{path}preprocessing/report_filtered.tsv',sep='\t')
    print('Filtering complete')
    return df, contams, filtered_set

#create preprocessing directory for new files 
def create_preprocessing_directory(path):
    # Combine the paths
    new_folder_path = os.path.join(path, 'preprocessing')
    
    # Create the new folder
    if not os.path.exists(new_folder_path):
        os.makedirs(new_folder_path)
        print(f"Folder preprocessing created successfully at {new_folder_path}")
    else:
        print(f"Folder preprocessing already exists at {new_folder_path}")
        
##Filtering
def remove_contaminants(chunk):
    # Create a contaminants mask based on the cont_ string and make sure all values are boolean
    contams_mask = chunk['Protein.Group'].str.contains('cont_', case=False, na=False)
    if not all(isinstance(x, bool) for x in contams_mask):
        print("contams_mask contains non-boolean values:", contams_mask[~contams_mask.isin([True, False])])

    contams_df = chunk[contams_mask]  # Dataframe with only contaminants
    cleaned_chunk = chunk[~contams_mask]  # Dataframe without contaminants
    return cleaned_chunk, contams_df
    
def apply_filters(chunk, params):
    # Initialize operator dict
    ops = {
        "==": operator.eq,
        "<": operator.lt,
        "<=": operator.le,
        ">": operator.gt,
        ">=": operator.ge
    }
    # Assign filtering parameter values and opperators to filtering conditions
    filtering_condition = (
        ops[params['apply_filters']['Global.PG.Q.Value']["op"]](chunk['Global.PG.Q.Value'], params['apply_filters']['Global.PG.Q.Value']["value"]) &
        ops[params['apply_filters']['Global.Q.Value']["op"]](chunk['Global.Q.Value'], params['apply_filters']['Global.Q.Value']["value"]) &
        ops[params['apply_filters']['Precursor.Charge']["op"]](chunk['Precursor.Charge'], params['apply_filters']['Precursor.Charge']["value"]) &
        ops[params['apply_filters']['Channel.Q.Value']["op"]](chunk['Channel.Q.Value'], params['apply_filters']['Channel.Q.Value']["value"])
    )
    # Filter chunk and return both filtered and filtered out dfs
    chunk_filtered = chunk[filtering_condition]
    chunk_filtered_out = chunk[~filtering_condition]

    return chunk_filtered, chunk_filtered_out

def  drop_cols(chunk):
    chunk['Genes'] = chunk['Genes'].fillna('')
    chunk['Protein.Group'] = chunk['Protein.Group'].str.cat(chunk['Genes'], sep='-')
    cols_to_keep = [ 'Run',
                     'Protein.Group',
                     'Stripped.Sequence',
                     'Precursor.Id', 
                     'Precursor.Charge',
                     'Lib.PG.Q.Value',
                     'Precursor.Quantity',
                     'Precursor.Translated',
                     'Ms1.Translated'
                     ]
    chunk = chunk[cols_to_keep]
    return chunk

