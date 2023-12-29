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

from silac_dia_tools.pipeline_new_version.silac_formatter_dev import SilacFormatter


from icecream import ic


class Preprocessor:
    def __init__(self, path, params, meta=None):
        self.path = path
        self.params = params
        self.meta = meta
        self.filter_cols = list(self.params.keys())
        self.combined_precursors = None
        # self.config_dir = config_dir or os.path.join(os.path.dirname(__file__), '..', 'configs')
        # self.params = self._load_params()
        self.chunk_size = 10000
        self.relable_with_meta = self._confirm_metadata()
        if self.relable_with_meta:
            self.meta_data = pd.read_csv(f'{self.path}{self.meta}', sep=',')
            print('Will relabel runs with metadata sample column')
        self.update = True
        self.formatter = SilacFormatter(self.path, self.filter_cols)
        
        # filtered dataframes
        self.filter_precursors = None
        self.filtered_out_stringent = None
        self.filtered_out_loose = None
    
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
    
    def import_and_format(self):
        df = self.import_tsv()
        self.combined_precursors = self.formatter.format_silac_channels(df)
        
        return self.combined_precursors
        
    def import_tsv(self):
        # import tsv from path, drop all rows not in the metadata, drop all cols that are not needed, return imported data
        print('Beggining import')
        if self.meta:
            print('reading in meta')
            metadata = pd.read_csv(f'{self.path}{self.meta}', sep=',')
        count = 1
        with open(f"{self.path}report.tsv", 'r', encoding='utf-8') as file:
            chunks = []
    
            for chunk in pd.read_table(file,sep="\t", chunksize=self.chunk_size):
                if self.meta:
                    chunk = self.drop_non_meta_samples(chunk, metadata)
                # chunk = self.drop_cols(chunk, self.filter_cols)
                print(chunk['Run'])
                chunks.append(chunk)
                
                # Update progress (optional)
                if self.update:
                    print('chunk ', count,' processed')
                count+=1
            
            # Concatenate all chunks into a DataFrames
  
            df = pd.concat(chunks, ignore_index=True)
            print('Finished import ')

        return df
   
    
    def drop_non_meta_samples(self, chunk, meta):
        filtered_chunk = chunk[chunk['Run'].isin(meta['Run'])]
        # Create a mapping from 'Run' to 'Sample' using the meta DataFrame
        run_to_sample_map = dict(zip(meta['Run'], meta['Sample']))
    
           # Map the 'Run' column of the chunk to the 'Sample' values from meta
        filtered_chunk['Run'] = filtered_chunk['Run'].map(run_to_sample_map)

        return filtered_chunk
        
    
    def filter_formatted(self):
        precursors, contam = self.remove_contaminants(self.combined_precursors)
        self.filtered_precursors, filtered_out = self.apply_filters(precursors)
        self.filtered_precursors = self.format_table(self.filtered_precursors)
        return self.filtered_precursors, contam, filtered_out
    
    #Filtering
    def remove_contaminants(self, combined_precursors): # is self needed?
        # Create a contaminants mask based on the cont_ string and make sure all values are boolean
        contams_mask = combined_precursors['Protein.Group'].str.contains('cont_', case=False, na=False)
        if not all(isinstance(x, bool) for x in contams_mask):
            print("contams_mask contains non-boolean values:", contams_mask[~contams_mask.isin([True, False])])

        contams_df = combined_precursors[contams_mask]  # Dataframe with only contaminants
        filtered_precursors = combined_precursors[~contams_mask]  # Dataframe without contaminants
        return filtered_precursors, contams_df
    
    def format_table(self, precursors):
        print('before formatting')
        ic(precursors)
        precursors = self.drop_filter_cols(precursors)
        precursors = self.stack_intensities(precursors)
        print('after formatting')
        ic(precursors)
        return precursors
        
    def drop_filter_cols(self, precursors):
        # Base columns to retain
        base_columns = ['Run', 'Protein.Group', 'Precursor', 'quantity_type','passed_stringent']
        
        # Find all columns that contain 'intensity' with any suffix
        intensity_columns = [col for col in precursors.columns if 'intensity_' in col]
     
        # Combine base columns and intensity columns
        columns_to_keep = base_columns + intensity_columns
     
        # Filter the DataFrame to retain only the specified columns
        precursors = precursors[columns_to_keep]
        return precursors
        
    def stack_intensities(self, df):
        # Check for each intensity column and create/fill with 0 if not present
        for col in ['intensity_L', 'intensity_H', 'intensity_M']:
            if col not in df.columns:
                df[col] = 0
        # Replace missing values with 0 and sum up the intensity columns
        df['Precursor.Quantity'] = df[['intensity_L', 'intensity_H', 'intensity_M']].fillna(0).sum(axis=1)
        
        # Create ratio columns, handling division by zero by replacing NaNs with 0
        df['H_to_stack_ratio'] = (df['intensity_H'] / df['Precursor.Quantity']).fillna(0)
        df['L_to_stack_ratio'] = (df['intensity_L'] / df['Precursor.Quantity']).fillna(0)
        df['M_to_stack_ratio'] = (df['intensity_M'] / df['Precursor.Quantity']).fillna(0)

        return df
    
    # def filter_formatted(self, formatted_precursors):
    #     precursors = self.relable_run(formatted_precursors)
    #     precursors, contam = self.remove_contaminants(precursors)
    #     precursors, filtered_out = self.apply_filters(precursors)
     
    #     return precursors, contam, filtered_out
        
     
    # def apply_filters(self, chunk):
    #     # Initialize operator dict
    #     ops = {
    #         "==": operator.eq,
    #         "<": operator.lt,
    #         "<=": operator.le,
    #         ">": operator.gt,
    #         ">=": operator.ge
    #     }

    #      # Create a boolean Series with all True values and explicitly set its index
    #     filtering_condition = pd.Series([True] * len(chunk), index=chunk.index)
        
    #     # Iterating over each filter condition in params['apply_filters']
    #     for column, condition in self.params.items():
    #         op = ops[condition['op']]
    #         value = condition['value']
            
    #         # Updating filtering_condition by applying each condition
    #         filtering_condition &= op(chunk[column], value)

    #     # Filter chunk and return both filtered and filtered out dfs
    #     chunk_filtered = chunk[filtering_condition]
    #     chunk_filtered_out = chunk[~filtering_condition]
        
    #     #filter for each channel and add boolean column to say whether all channels passed
        
    #     return chunk_filtered, chunk_filtered_out
    
    def apply_filters(self, chunk):
        # Initialize operator dict
        ops = {
            "==": operator.eq,
            "<": operator.lt,
            "<=": operator.le,
            ">": operator.gt,
            ">=": operator.ge
        }
    
        # Create a boolean Series for initial filtering
        filtering_condition = pd.Series([True] * len(chunk), index=chunk.index)
        
        # Iterating over each filter condition in params['apply_filters']
        for column, condition in self.params.items():
            op = ops[condition['op']]
            value = condition['value']
            
            # Updating filtering_condition by applying each condition
            filtering_condition &= op(chunk[column], value)
    
        # Filter chunk based on initial conditions
        chunk_filtered = chunk[filtering_condition]
        chunk_filtered_out = chunk[~filtering_condition]
    
        # Check for additional columns and apply further filtering
        stringent_passed = pd.Series([True] * len(chunk_filtered), index=chunk_filtered.index)
        for column, condition in self.params.items():
            for suffix in ['_L', '_M', '_H']:
                additional_column = f'{column}{suffix}'
                if additional_column in chunk_filtered.columns:
                    op = ops[condition['op']]
                    value = condition['value']
                    stringent_passed &= op(chunk_filtered[additional_column], value)
    
        # Adding 'passed stringent' column
        chunk_filtered['passed_stringent'] = stringent_passed
        
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
 