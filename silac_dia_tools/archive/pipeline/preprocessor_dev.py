"""
Created on Mon Sep 18 11:43:50 2023

@author: rkerrid

Step 1: Module for filtering report.tsv output (DIA-NN version 1.8.1) with
SILAC settings as described in the README.md

Note: This script filters for contaminants by looking for the 'cont_' substring
in Protein.Groups so make sure your report.tsv is annotated in the same way or 
edit the remove_contaminants() funciton.

"""

import json
import os
import operator
from icecream import ic
import pandas as pd



class Preprocessor:
    def __init__(self, path, parameter_file, meta=None, config_dir=None):
        self.path = path
        self.chunk_size = 10000
        self.parameter_file = parameter_file
        self.meta = meta
        self.config_dir = config_dir or os.path.join(os.path.dirname(__file__), '..', 'configs')
        self.params = self._load_params()
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
            
        
        print('Filtering complete')
        return df, contams, filtered_set


    def remove_contaminants(self, df):
         # Adapt this method for Dask
         contams_mask = df['Protein.Group'].str.contains('cont_', case=False, na=False)
         contams_df = df[contams_mask]
         cleaned_df = df[~contams_mask]
         return cleaned_df, contams_df

    def apply_filters(self, df):
        # Initialize operator dict
        ops = {
            "==": operator.eq,
            "<": operator.lt,
            "<=": operator.le,
            ">": operator.gt,
            ">=": operator.ge
        }
    
        # Initialize the first condition
        first_column, first_condition = next(iter(self.params['apply_filters'].items()))
        op = ops[first_condition['op']]
        value = first_condition['value']
        filtering_condition = op(df[first_column], value)
    
        # Apply subsequent filter conditions
        for column, condition in list(self.params['apply_filters'].items())[1:]:
            op = ops[condition['op']]
            value = condition['value']
            current_condition = op(df[column], value)
            filtering_condition = filtering_condition & current_condition
    
        # Filter the DataFrame using the combined conditions
        df_filtered = df[filtering_condition]
        df_filtered_out = df[~filtering_condition]
        return df_filtered, df_filtered_out

    def drop_cols(self, df):
         # Adapt this method for Dask
         df['Genes'] = df['Genes'].fillna('')
         df['Protein.Group'] = df['Protein.Group'].str.cat(df['Genes'], sep='-')
         cols_to_keep = ['Run', 'Protein.Group', 'Stripped.Sequence', 'Precursor.Id', 
                         'Precursor.Charge', 'Lib.PG.Q.Value', 'Precursor.Quantity', 
                         'Precursor.Translated', 'Ms1.Translated']
         return df[cols_to_keep]

    def relable_run(self, df):
        # Adapt this method for Dask
        run_to_sample = dict(zip(self.meta_data['Run'].tolist(), self.meta_data['Sample'].tolist()))
        df['Run'] = df['Run'].map(run_to_sample)
        if df['Run'].isnull().any():
            raise ValueError("Some Run values in the report.tsv are not found in the metadata, please ensure metadata is correct.")
        return df


def _import_filtered(self, path):
    # Iterate through the file in chunks and apply preprocessing functions
    print('Beggining filtering')
    count = 1
    with open(f"{self.path}report_filtered.tsv", 'r', encoding='utf-8') as file:
        chunks = []
        

        for chunk in pd.read_table(file,sep="\t", chunksize=10000):

            chunks.append(chunk)
            
            # Update progress (optional)
            if self.update:
                print('chunk ', count,' processed')
            count+=1
        
        # Concatenate all chunks into a DataFrames
      
        df = pd.concat(chunks, ignore_index=True)
        return df

# import pandas as pd
# import json
# import os
# import operator
# import dask.dataframe as dd
# import datatable as dt
# from icecream import ic




# class Preprocessor:
#     def __init__(self, path, parameter_file, meta=None, config_dir=None):
#         self.path = path
#         self.parameter_file = parameter_file
#         self.meta = meta
#         self.config_dir = config_dir or os.path.join(os.path.dirname(__file__), '..', 'configs')
#         self.params = self._load_params()
#         self.chunk_size = 100000
#         self.relable_with_meta = self._confirm_metadata()
#         if self.relable_with_meta:
#             self.meta_data = pd.read_csv(f'{self.path}{self.meta}', sep=',')
#             print('Will relabel runs with metadata sample column')
#         self.update = True
        
#     def _load_params(self):
#         json_path = os.path.join(self.config_dir, self.parameter_file)
#         with open(json_path, 'r') as f:
#             params = json.load(f)
#             return params 
    
#     def _confirm_metadata(self):
#         if self.meta is None:
#             print("No metadata added, filtering will continue without relabeling")
#             return False
#         elif isinstance(self.meta, str):
#             print("Metadata added, looking for the following file:", self.meta)
#             meta_exists = self._check_directory()
#             if meta_exists:
#                 return True
#             else:
#                 print(f"Cannot find {self.meta} in {self.path}, check file name and location") 
#                 print("Filtering will continue without relabeling")
#         else:
#             print("File name is not a string, filering will continue without relabeling")
#             return False
    
#     def _check_directory(self):
#         file_list = os.listdir(self.path)
#         # Iterate through the list of filenames and check for a match
#         found = False
#         for filename in file_list:
#             if filename == self.meta:
#                 found = True
#                 print(f"CSV file '{self.meta}' found in {self.path}")
#                 return True
#         if not found:
#             print(f"CSV file '{self.meta}' not found in the directory.")
#             return False
    
#     def import_and_filter(self):
#         # Iterate through the file in chunks and apply preprocessing functions
#         print('Beggining filtering')
#         count = 1
#         with open(f"{self.path}report.tsv", 'r', encoding='utf-8') as file:
#             chunks = []
#             contams = []
#             filtered_set = []
    
#             for chunk in pd.read_table(file,sep="\t", chunksize=self.chunk_size):
#                 # if applicable, relable chunk Run colum with metadata
#                 if self.relable_with_meta:
#                     chunk = self.relable_run(chunk)
                
#                 # Apply filtering to each chunk
#                 chunk, contam = self.remove_contaminants(chunk)
#                 chunk, filtered_out = self.apply_filters(chunk)
                
#                 # Drop unnecesary columns
#                 chunk = self.drop_cols(chunk) 
                
#                 # Append chunks from respective filtering steps
#                 filtered_set.append(filtered_out)
#                 contams.append(contam)
#                 chunks.append(chunk)
                
#                 # Update progress (optional)
#                 if self.update:
#                     print('chunk ', count,' processed')
#                 count+=1
            
#             # Concatenate all chunks into a DataFrames
#             filtered_set = pd.concat(filtered_set, ignore_index=True)
#             contams = pd.concat(contams, ignore_index=True)
#             df = pd.concat(chunks, ignore_index=True)
            
        
#         # # Apply filtering to each chunk
#         # df, contam = self.remove_contaminants(df)
#         # df, filtered_out = self.apply_filters(df)
        
#         # # Drop unnecesary columns
#         # df = self.drop_cols(df) 
        
#         print('Filtering complete')
#         return df, contams, filtered_set

#     #Filtering
#     def remove_contaminants(self, chunk): # is self needed?
#         # Create a contaminants mask based on the cont_ string and make sure all values are boolean
#         contams_mask = chunk['Protein.Group'].str.contains('cont_', case=False, na=False)
#         if not all(isinstance(x, bool) for x in contams_mask):
#             print("contams_mask contains non-boolean values:", contams_mask[~contams_mask.isin([True, False])])

#         contams_df = chunk[contams_mask]  # Dataframe with only contaminants
#         cleaned_chunk = chunk[~contams_mask]  # Dataframe without contaminants
#         return cleaned_chunk, contams_df
        
#     def apply_filters(self, chunk):
#         # Initialize operator dict
#         ops = {
#             "==": operator.eq,
#             "<": operator.lt,
#             "<=": operator.le,
#             ">": operator.gt,
#             ">=": operator.ge
#         }

#          # Create a boolean Series with all True values and explicitly set its index
#         filtering_condition = pd.Series([True] * len(chunk), index=chunk.index)
        
#         # Iterating over each filter condition in params['apply_filters']
#         for column, condition in self.params['apply_filters'].items():
#             op = ops[condition['op']]
#             value = condition['value']
            
#             # Updating filtering_condition by applying each condition
#             filtering_condition &= op(chunk[column], value)

#         # Filter chunk and return both filtered and filtered out dfs
#         chunk_filtered = chunk[filtering_condition]
#         chunk_filtered_out = chunk[~filtering_condition]

#         return chunk_filtered, chunk_filtered_out

#     def  drop_cols(self, chunk): # is self needed?
#         chunk['Genes'] = chunk['Genes'].fillna('')
#         chunk['Protein.Group'] = chunk['Protein.Group'].str.cat(chunk['Genes'], sep='-')
#         cols_to_keep = [ 'Run',
#                           'Protein.Group',
#                           'Stripped.Sequence',
#                           'Precursor.Id', 
#                           'Precursor.Charge',
#                           'Lib.PG.Q.Value',
#                           'Precursor.Quantity',
#                           'Precursor.Translated',
#                           'Ms1.Translated'
#                           ]
#         chunk = chunk[cols_to_keep]
#         return chunk
    
#     def relable_run(self, chunk):
#         run_to_sample = dict(zip(self.meta_data['Run'], self.meta_data['Sample']))

#         # Apply the mapping to df2['Run'] and raise an error if a 'Run' value doesn't exist in df1
#         chunk['Run'] = chunk['Run'].map(run_to_sample)
#         if chunk['Run'].isna().any():
#             raise ValueError("Some Run values in the report.tsv are not found in the metadata, please ensure metadata is correct.")
            
#         return chunk
 