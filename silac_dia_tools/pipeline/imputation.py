# -*- coding: utf-8 -*-
"""
Created on Sun Nov 12 11:47:28 2023

@author: robbi
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 18:41:09 2023

@author: robbi
"""


import pandas as pd
import numpy as np
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
from icecream import ic

path = 'G:/My Drive/Data/data/no spikein data/'
# import metadata
def import_meta(path):
    metadata = pd.read_csv(f"{path}meta.csv")
    return metadata

def get_dataframes(path):
    total = pd.read_csv(f"{path}protein intensities/total_dlfq.csv", sep=',')
    light = pd.read_csv(f"{path}protein intensities/light_dlfq.csv", sep=',')
    nsp = pd.read_csv(f"{path}protein intensities/nsp_dlfq.csv", sep=',')
    return total, light, nsp
    

def replace_values(df):
        df.replace([np.inf, -np.inf], np.nan, inplace=True)
        df.iloc[:,1:] = df.iloc[:,1:].apply(lambda x: np.where(x < 0.001, np.nan, x))
        return df
        
def filter_for_valid_values(df, metadata):
    # Initialize a new column in df to track if the row should be kept
    df['keep_row'] = False
    
    for group in metadata['Treatment'].unique():
        sub_meta = metadata[metadata["Treatment"]==group]
        cols = sub_meta['Sample'].tolist()
        print(cols)
        
    
    # # Create a set of base column names by splitting the replicate numbers off and removing the 'Protein.Group'
    # base_columns = set(col.rsplit(' ', 1)[0] for col in df.columns if ' ' in col and col not in ['Protein.Group'])

    # # Iterate over each unique sample name (without replicate number)
    # for base in base_columns:
    #     col1 = f'{base} 1'
    #     col2 = f'{base} 2'
        
    #     # Check if both replicates have non-NaN values
    #     df['keep_row'] |= df[[col1, col2]].notna().all(axis=1)
    
    # # Select rows where 'keep_row' is True and drop the 'keep_row' column to clean up
    # df_filtered = df[df['keep_row']].drop(columns=['keep_row'])
    
    return df

def drop_non8h_cols(df):
    columns_to_keep = [col for col in df.columns if '8h' in col or col == 'Protein.Group']
    
    # Create a new DataFrame with only the filtered columns
    df_filtered = df[columns_to_keep]
    return df_filtered

# imputation
def create_distribution(data):
    mu, std = stats.norm.fit(data)
    mu = mu - 1.8 * std
    std = 0.25*std
   
    return mu, std
    
def impute(x, mu, std, imputed_only):
    
    if imputed_only:
        if np.isnan(x):
            return np.random.normal(mu, std)
        else: 
            return np.nan
        
    else:
        if np.isnan(x):
            return np.random.normal(mu, std)
        else:
            return x

def perform_imputation(df):
    cols = df.columns.values.tolist()[1:]
    imputed_values = pd.DataFrame(columns=df.columns.values.tolist())
    for condition in cols:
        data = df[condition].dropna()
        ic(np.all(np.isfinite(data)))
        ic(data)
        mu, std = create_distribution(data)
        
        imputed_values[condition] = df[condition].apply(lambda x: impute(x, mu, std, True))
        df[condition] = df[condition].apply(lambda x: impute(x, mu, std, False))
    
    return  df, imputed_values


def plot_histogram(df, imputed_values):
  # Plot histograms for each column in df1 and df2 on the same plot with different colors
    # Create a single figure with 2 rows and 3 columns
  
    # Plot histograms for each column in df and imputed_values
    for  col in df.columns.values.tolist()[1:]:
        

        plt.hist(df[col], bins=20, alpha=0.5, label='original data', color='blue')
        plt.hist(imputed_values[col], bins=20, alpha=0.5, label='imputed values', color='green')
        plt.title(col + ' Histogram')
        plt.xlabel('Log2 intensity')
        plt.ylabel('Frequency')
        # plt.legend()

  

        # Display the figure
        plt.show()

# import meta and dataframes
metadata = import_meta(path)
groups = metadata['Treatment'].unique()
total, light, nsp = get_dataframes(path)
# replace NaN and inf values
total = replace_values(total)
light = replace_values(light)
nsp = replace_values(nsp)
# filter for valid values

filter_for_valid_values(total, metadata)



# # Example usage with the provided DataFrame loading code
# nsp_df = pd.read_csv("G:/My Drive/Data/data/eIF4F optimization/protein intensities/nsp_href.csv", sep=',')
# total_df = pd.read_csv("G:/My Drive/Data/data/eIF4F optimization/protein intensities/total_href.csv", sep=',')
# light_df = pd.read_csv("G:/My Drive/Data/data/eIF4F optimization/protein intensities/light_href.csv", sep=',')

# # nsp_df = pd.read_csv("G:/My Drive/Data/data/no spikein data/protein intensities/nsp_dlfq.csv", sep=',')
# # total_df = pd.read_csv("G:/My Drive/Data/data/no spikein data/protein intensities/total_dlfq.csv", sep=',')
# # light_df = pd.read_csv("G:/My Drive/Data/data/no spikein data/protein intensities/light_dlfq.csv", sep=',')



# nsp_df = drop_non8h_cols(nsp_df)
# total_df = drop_non8h_cols(total_df)
# light_df = drop_non8h_cols(light_df)

# # replace = and -inf 
# total_df.replace([np.inf, -np.inf], np.nan, inplace=True)
# nsp_df.replace([np.inf, -np.inf], np.nan, inplace=True)
# light_df.replace([np.inf, -np.inf], np.nan, inplace=True)

# # Replace low values with NaN for filtering
# total_df.iloc[:,1:] = total_df.iloc[:,1:].apply(lambda x: np.where(x < 0.001, np.nan, x))
# nsp_df.iloc[:,1:] = nsp_df.iloc[:,1:].apply(lambda x: np.where(x < 0.001, np.nan, x))
# light_df.iloc[:,1:] = light_df.iloc[:,1:].apply(lambda x: np.where(x < 0.001, np.nan, x))

# # log2 before filtering and imputation
# total_df.iloc[:,1:] = np.log2(total_df.iloc[:,1:])
# nsp_df.iloc[:,1:] = np.log2(nsp_df.iloc[:,1:])
# light_df.iloc[:,1:] = np.log2(light_df.iloc[:,1:])

# # Filter the DataFrames
# total_df = filter_for_valid_values(total_df)
# nsp_df = filter_for_valid_values(nsp_df)
# light_df = filter_for_valid_values(light_df)

# # impute with gausian shift
# total_df, total_df_imputed = perform_imputation(total_df)
# nsp_df, nsp_df_imputed = perform_imputation(nsp_df)
# light_df, light_df_imputed = perform_imputation(light_df)

# plot_histogram(total_df, total_df_imputed)
# plot_histogram(nsp_df, nsp_df_imputed)
# plot_histogram(light_df, light_df_imputed)

# # base 2 exponentiation before saving
# total_df.iloc[:,1:] = 2**total_df.iloc[:,1:]
# nsp_df.iloc[:,1:] = 2**nsp_df.iloc[:,1:]
# light_df.iloc[:,1:] = 2**light_df.iloc[:,1:]

# nsp_df.to_csv("G:/My Drive/Data/data/eIF4F optimization/imputed intensities 8hr/nsp.csv", sep=',')
# light_df.to_csv("G:/My Drive/Data/data/eIF4F optimization/imputed intensities 8hr/light.csv", sep=',')
# total_df.to_csv("G:/My Drive/Data/data/eIF4F optimization/imputed intensities 8hr/total.csv", sep=',')

# # nsp_df.to_csv("G:/My Drive/Data/data/no spikein data/protein intensities/nsp.csv", sep=',')
# # light_df.to_csv("G:/My Drive/Data/data/no spikein data/protein intensities/light.csv", sep=',')
# # total_df.to_csv("G:/My Drive/Data/data/no spikein data/protein intensities/total.csv", sep=',')





