import pandas as pd
import numpy as np
from tqdm import tqdm
from icecream import ic


class PrecursorRollup:
    def __init__(self, path):
        self.path = path
        self.protein_groups = None
        
    @staticmethod
    def select_precursor_translated(group):
        return group[group['quantity_type'] == 'Precursor.Translated']

    @staticmethod
    def select_ms1_translated(group):
        return group[group['quantity_type'] == 'Ms1.Translated']

    def calculate_protein_level_ratios(self, df):
  
        print("Calculating ratios from precursor information")
        protein_precursors = df.groupby(['Run', 'Protein.Group'])
        
        protein_data = []
        protein_count = 0
        protein_missed = 0
        
        for name, group in tqdm(protein_precursors, desc="Processing proteins"):
            precursor_group = self.select_precursor_translated(group)
            ms1_group = self.select_ms1_translated(group)
          
            if len(precursor_group) > 2 and len(ms1_group) > 2:
                # Compute the median of log2 ratios
                median_log2_ratio_l = np.median(np.log2(group['L_to_stack_ratio']))
                median_log2_ratio_m = np.median(np.log2(group['M_to_stack_ratio']))
                median_log2_ratio_h = np.median(np.log2(group['H_to_stack_ratio']))
                # Reverse log
                median_ratio_l = np.exp2(median_log2_ratio_l)
                median_ratio_m = np.exp2(median_log2_ratio_m)
                median_ratio_h = np.exp2(median_log2_ratio_h)
               
                total_intensity = np.sum(ms1_group['Precursor.Quantity'])
                # Create new row to add to new dataframe containing ratios and precursor translated quantity
                new_row = {
                    'Run': group['Run'].iloc[0],
                    'Protein.Group': group['Protein.Group'].iloc[0],
                   
                    'Total_intensity': total_intensity,
                    'L_to_stack_ratio': median_ratio_l,
                    'M_to_stack_ratio': median_ratio_m,
                    'H_to_stack_ratio': median_ratio_h,
                    'L_intensity': median_ratio_l * total_intensity,
                    'M_intensity': median_ratio_m * total_intensity,
                    'H_intensity': median_ratio_h * total_intensity
                    
                        
                }
                protein_data.append(new_row)
                protein_count += 1
               
            else:
                protein_missed += 1
        print('Total proteins counted ', protein_count)
        print('Total sets of precursors that didnt meet mimimum unique precursor requirements ', protein_missed, ' out of ', len(protein_precursors))
           
        protein_ratios = pd.DataFrame(protein_data)
        
        self.protein_groups = protein_ratios
        return self.protein_groups 
    
    def calculate_stringent_and_inclusive(self, df):
        # create deep coppies of df
       stringent = df.copy(deep=True)
       inclusive = df.copy(deep=True)
       # use 'passed stringent' to subset df based on filtering
       stringent = stringent[stringent['passed_stringent']]
       ic(stringent)
       inclusive = inclusive[inclusive['passed_stringent']==False]
       ic(inclusive)
       # calculate protein level ratios for each filtering algorythm
       stringent_protein_ratios = self.calculate_protein_level_ratios(stringent)
       inclusive_protein_ratios = self.calculate_protein_level_ratios(inclusive)
       
       return stringent_protein_ratios, inclusive_protein_ratios
   
    
   
    
   
    