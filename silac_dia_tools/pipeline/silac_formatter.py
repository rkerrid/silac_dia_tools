
import pandas as pd
import numpy as np
from silac_dia_tools.pipeline.report import precursor_report


class SilacFormatter:
    def __init__(self, path):
        self.path = path
        self.filtered_report = None
        self.parsed_df = None
        self.combined_precursors = None
        self.stacked_intensities = None
        
    def format_silac_channels(self, filtered_report):
        self.parsed_df = self.parse_data_for_channel_info(filtered_report)
        self.combined_precursors = self.combine_modified_precursors()
        self.stacked_intensities = self.stack_intensities()
        
        # print('Generating report...')
        # precursor_report.silac_precursor_qc(self.df, self.path)
        # print('Saving precursors...')
        # new_path = f'{self.path}preprocessing/'
        # self.stacked_intensities.to_csv(new_path+'silac_precursors.tsv', sep='\t')
        print('Done!')
        return self.stacked_intensities


    def parse_data_for_channel_info(self, filtered_report):
        print('Parsing data for SILAC intensities')
        # Extracting SILAC labels from Precursor.ID and adding to separate column
        filtered_report['Label'] = filtered_report['Precursor.Id'].str.extract(r'\(SILAC-(K|R)-([HML])\)')[1]    
        filtered_report['Precursor'] = filtered_report['Stripped.Sequence'] + filtered_report['Precursor.Charge'].astype(str)
        
        # Splitting the data into two, one for 'Ms1.Translated' and another for 'Precursor.Translated'
        ms1_translated_df = filtered_report.copy()
        ms1_translated_df['intensity'] = filtered_report['Ms1.Translated']
        ms1_translated_df['quantity type'] = 'Ms1.Translated'
    
        precursor_translated_df = filtered_report.copy()
        precursor_translated_df['intensity'] = filtered_report['Precursor.Translated']
        precursor_translated_df['quantity type'] = 'Precursor.Translated'
        
        # Concatenate the two dataframes to create the parsed dataframe
        parsed_df = pd.concat([ms1_translated_df, precursor_translated_df], ignore_index=True)
        return parsed_df[['Run', 'Protein.Group', 'Precursor.Id', 'Precursor', 'Label', 'intensity', 'quantity type', 'Precursor.Quantity', 'Lib.PG.Q.Value']]

    def combine_modified_precursors(self):
        print('combining modified precursors')
        # Define aggregation functions for columns
        agg_functions = {
            'intensity': 'first',
            'Lib.PG.Q.Value': 'first',
        }
        # Aggregate data using groupby and agg function
        agg_df = self.parsed_df.groupby(['Run', 'Protein.Group',  'Precursor', 'quantity type', 'Label']).agg(agg_functions).reset_index()

        # Pivot the table to get 'H', 'M', 'L' intensities in separate columns
        pivoted_df = agg_df.pivot_table(index=['Run', 'Protein.Group', 'Precursor', 'quantity type'], columns='Label', values='intensity', fill_value=0).reset_index()

        # Rename columns based on the label
        pivoted_df.columns.name = None
        pivoted_df = pivoted_df.rename(columns={'H': 'H intensity', 'M': 'M intensity', 'L': 'L intensity'})

        # Merge the pivoted_df with the original df to get other columns back
        df = pd.merge(agg_df.drop(columns=['Label', 'intensity']), pivoted_df, on=['Run', 'Protein.Group',  'Precursor', 'quantity type'])

        df = df.drop_duplicates()
        return df

    def stack_intensities(self):
        df = self.combined_precursors
        columns_to_check = ['H intensity', 'M intensity', 'L intensity']
        for col in columns_to_check:
            if col not in df.columns:
                df[col] = 0
        df['Precursor.Quantity'] = df['H intensity'] + df['M intensity'] + df['L intensity']
        df['L to stack ratio'] = df['L intensity']/df['Precursor.Quantity']
        df['M to stack ratio'] = df['M intensity']/df['Precursor.Quantity']
        df['H to stack ratio'] = df['H intensity']/df['Precursor.Quantity']
        df = self.drop_nan(df)
        return df


    def drop_nan(self, df):
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

  








# def combine_modified_precursors(df):
#     print('combining modified precursors')
#     # Define aggregation functions for columns
#     agg_functions = {
#         'intensity': 'first',
#         'Lib.PG.Q.Value': 'first',
#     }
#     # Aggregate data using groupby and agg function
#     agg_df = df.groupby(['Run', 'Protein.Group',  'Precursor', 'quantity type', 'Label']).agg(agg_functions).reset_index()

#     # Pivot the table to get 'H', 'M', 'L' intensities in separate columns
#     pivoted_df = agg_df.pivot_table(index=['Run', 'Protein.Group', 'Precursor', 'quantity type'], columns='Label', values='intensity', fill_value=0).reset_index()

#     # Rename columns based on the label
#     pivoted_df.columns.name = None
#     pivoted_df = pivoted_df.rename(columns={'H': 'H intensity', 'M': 'M intensity', 'L': 'L intensity'})

#     # Merge the pivoted_df with the original df to get other columns back
#     df = pd.merge(agg_df.drop(columns=['Label', 'intensity']), pivoted_df, on=['Run', 'Protein.Group',  'Precursor', 'quantity type'])

#     df = df.drop_duplicates()
#     return df


# #Stacking intensities and calculating intensity to stack ratio
# def stack_intensities(df):
#     columns_to_check = ['H intensity', 'M intensity', 'L intensity']
#     for col in columns_to_check:
#         if col not in df.columns:
#             df[col] = 0
#     df['Precursor.Quantity'] = df['H intensity'] + df['M intensity'] + df['L intensity']
#     df['L to stack ratio'] = df['L intensity']/df['Precursor.Quantity']
#     df['M to stack ratio'] = df['M intensity']/df['Precursor.Quantity']
#     df['H to stack ratio'] = df['H intensity']/df['Precursor.Quantity']
#     return df

# def drop_nan(df):
#     cols_to_check = [
#     'Precursor.Quantity', 'H intensity', 'L intensity', 'M intensity', 
#     'L to stack ratio', 'M to stack ratio', 'H to stack ratio'
#     ]
#     # Create a mask where all values in the specified columns are either 0 or NaN
#     mask = df[cols_to_check].isin([0, np.nan]).all(axis=1)
    
#     # Use the mask to drop rows from the DataFrame
#     df = df[~mask]
    
#     #rename precursor for directLFQ
#     df = df.rename(columns={'Precursor':'Precursor.Id'})
#     return df[['Run', 'Protein.Group', 'Precursor.Id','Precursor.Quantity','Lib.PG.Q.Value','quantity type','H intensity','M intensity','L intensity','H to stack ratio', 'M to stack ratio', 'L to stack ratio']]






    
