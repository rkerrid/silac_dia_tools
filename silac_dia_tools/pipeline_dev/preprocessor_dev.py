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
import json
import os
import operator

from icecream import ic


class Preprocessor:
    def __init__(self, path, parameter_file, meta=None, config_dir=None):
        self.path = path
        self.parameter_file = parameter_file
        self.meta = meta
        self.config_dir = config_dir or os.path.join(os.path.dirname(__file__), '..', 'configs')
        self.params = self._load_params()
        self.chunk_size = 10000
        self.relable_with_meta = self._confirm_metadata()
        if self.relable_with_meta:
            self.meta_data = pd.read_csv(f'{self.path}{self.meta}', sep=',')
            print('Will relabel runs with metadata sample column')
        self.update = True
        
    def _load_params(self):
        json_path = os.path.join(self.config_dir, self.parameter_file)
        with open(json_path, 'r') as f:
            params = json.load(f)
            return params 
    
    def _confirm_metadata(self):
        if self.meta is None:
            print("No metadata added, filtering will continue without relabeling")
            return False
        elif isinstance(self.meta, str):
            print("Metadata added, looking for the following file:", self.meta)
            meta_exists = self._check_directory()
            if meta_exists:
                return True
            else:
                print(f"Cannot find {self.meta} in {self.path}, check file name and location") 
                print("Filtering will continue without relabeling")
        else:
            print("File name is not a string, filering will continue without relabeling")
            return False
    
    def _check_directory(self):
        file_list = os.listdir(self.path)
        # Iterate through the list of filenames and check for a match
        found = False
        for filename in file_list:
            if filename == self.meta:
                found = True
                print(f"CSV file '{self.meta}' found in {self.path}")
                return True
        if not found:
            print(f"CSV file '{self.meta}' not found in the directory.")
            return False
    
    def import_no_filter(self, filter_cols):
        print('Beggining import no filter')
        if self.meta:
            print('reading in meta')
            metadata = pd.read_csv(f'{self.path}{self.meta}', sep=',')
        count = 1
        with open(f"{self.path}report.tsv", 'r', encoding='utf-8') as file:
            chunks = []
    
    
            for chunk in pd.read_table(file,sep="\t", chunksize=self.chunk_size):
                # print(chunk['Run'])
                if self.meta:
                    chunk = self.drop_non_meta_samples(chunk, metadata)
                    chunk['Genes'] = chunk['Genes'].fillna('')
                    chunk['Protein.Group'] = chunk['Protein.Group'].str.cat(chunk['Genes'], sep='-')
                # chunk = self.drop_cols(chunk, filter_cols)
                chunks.append(chunk)
                
                # Update progress (optional)
                if self.update:
                    print('chunk ', count,' processed')
                count+=1
            
            # Concatenate all chunks into a DataFrames
  
            df = pd.concat(chunks, ignore_index=True)
            print('Finished import no filter')
        ic(df)
        print(df.columns.values.tolist())
        return df
    
    def drop_non_meta_samples(self, chunk, meta):
        filtered_chunk = chunk[chunk['Run'].isin(meta['Run'])]
        return filtered_chunk
        
    
    def filter_formatted(self, formatted_precursors):
        precursors = self.relable_run(formatted_precursors)
        precursors, contam = self.remove_contaminants(precursors)
        precursors, filtered_out = self.apply_filters(precursors)
     
        return precursors, contam, filtered_out
        
    def import_and_filter(self):
        # Iterate through the file in chunks and apply preprocessing functions
        print('Beggining filtering')
        count = 1
        with open(f"{self.path}report.tsv", 'r', encoding='utf-8') as file:
            chunks = []
            contams = []
            filtered_set = []
    
            for chunk in pd.read_table(file,sep="\t", chunksize=self.chunk_size):
                # if applicable, relable chunk Run colum with metadata
                if self.relable_with_meta:
                    chunk = self.relable_run(chunk)
                
                # Apply filtering to each chunk
                chunk, contam = self.remove_contaminants(chunk)
                chunk, filtered_out = self.apply_filters(chunk)
                
                # Drop unnecesary columns
                chunk = self.drop_cols(chunk) 
                
                # Append chunks from respective filtering steps
                filtered_set.append(filtered_out)
                contams.append(contam)
                chunks.append(chunk)
                
                # Update progress (optional)
                if self.update:
                    print('chunk ', count,' processed')
                count+=1
            
            # Concatenate all chunks into a DataFrames
            filtered_set = pd.concat(filtered_set, ignore_index=True)
            contams = pd.concat(contams, ignore_index=True)
            df = pd.concat(chunks, ignore_index=True)
            
        
        # # Apply filtering to each chunk
        # df, contam = self.remove_contaminants(df)
        # df, filtered_out = self.apply_filters(df)
        
        # # Drop unnecesary columns
        # df = self.drop_cols(df) 
        
        print('Filtering complete')
        return df, contams, filtered_set

    #Filtering
    def remove_contaminants(self, chunk): # is self needed?
        # Create a contaminants mask based on the cont_ string and make sure all values are boolean
        contams_mask = chunk['Protein.Group'].str.contains('cont_', case=False, na=False)
        if not all(isinstance(x, bool) for x in contams_mask):
            print("contams_mask contains non-boolean values:", contams_mask[~contams_mask.isin([True, False])])

        contams_df = chunk[contams_mask]  # Dataframe with only contaminants
        cleaned_chunk = chunk[~contams_mask]  # Dataframe without contaminants
        return cleaned_chunk, contams_df
        
    def apply_filters(self, chunk):
        # Initialize operator dict
        ops = {
            "==": operator.eq,
            "<": operator.lt,
            "<=": operator.le,
            ">": operator.gt,
            ">=": operator.ge
        }

         # Create a boolean Series with all True values and explicitly set its index
        filtering_condition = pd.Series([True] * len(chunk), index=chunk.index)
        
        # Iterating over each filter condition in params['apply_filters']
        for column, condition in self.params['apply_filters'].items():
            op = ops[condition['op']]
            value = condition['value']
            
            # Updating filtering_condition by applying each condition
            filtering_condition &= op(chunk[column], value)

        # Filter chunk and return both filtered and filtered out dfs
        chunk_filtered = chunk[filtering_condition]
        chunk_filtered_out = chunk[~filtering_condition]

        return chunk_filtered, chunk_filtered_out

    def  drop_cols(self, chunk, filter_cols = []): # is self needed?
        chunk['Genes'] = chunk['Genes'].fillna('')
        chunk['Protein.Group'] = chunk['Protein.Group'].str.cat(chunk['Genes'], sep='-')
        cols_to_keep = [ 'Run',
                          'Protein.Group',
                          'Stripped.Sequence',
                          'Precursor.Id', 
                          'Precursor.Charge',
            
                          'Precursor.Quantity',
                          'Precursor.Translated',
                          'Ms1.Translated'
                          ] + filter_cols
        chunk = chunk[cols_to_keep]
        return chunk
    
    def relable_run(self, chunk):
        run_to_sample = dict(zip(self.meta_data['Run'], self.meta_data['Sample']))

        # Apply the mapping to df2['Run'] and raise an error if a 'Run' value doesn't exist in df1
        chunk['Run'] = chunk['Run'].map(run_to_sample)
        if chunk['Run'].isna().any():
            raise ValueError("Some Run values in the report.tsv are not found in the metadata, please ensure metadata is correct.")
            
        return chunk
 