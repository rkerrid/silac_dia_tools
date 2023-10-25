# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 16:51:01 2023

@author: rkerrid
"""
from silac_dia_tools.label_check import file_io as io
from silac_dia_tools.label_check import process_peptides as peptide
from silac_dia_tools.label_check import plotting
import pandas as pd


# work
path = 'C:/data/silac_dia_tools_files/data/spikein data/'
# home
path = 'G:/My Drive/Data/data/label check proof of concept/'

##Main
#path and file variables
raw_file = '20231010_AJW_HS_px07_HEK293_SILAC_heavycheck_3.raw'


#import MQ and thermo raw file data
msms = io.import_file(path + "msms.txt")
evidence = io.import_file(path + "evidence.txt")

#get peptide details
peptides = peptide.get_most_abundant(msms,evidence) 

# #get MS1 data  
ms1_df = io.get_raw_data(path, raw_file, 1000) # when accessing the same raw file, comment the import code accordingly

#combine data from MQ and thermo raw files with the option to save to a temp csv when adjusting plotting
for peptide_details in peptides:
    mz_values,intensity_values = peptide.extract_values(ms1_df, peptide_details)
    sequence = peptide_details['Sequence'].values[0]
    print(sequence)
    data = {"mz": mz_values,
            "intensity": intensity_values
            }
    df = pd.DataFrame(data)
    df.to_csv(path + "/temp data/" + sequence + ".csv", index=False)
    
for peptide_details in peptides:
    sequence = peptide_details["Sequence"].values[0]
    retention_time = peptide_details["Retention_time"].values[0]
    charge = peptide_details["Charge"].values[0]
    mz = peptide_details["m/z"].values[0]
    # adjuest mz by precursor apex offset
  
    AA_mass = peptide_details["AA_mass"].values[0]
    scan_number = peptide_details["scan_number"].values[0]
    expected_light_peak = mz-(AA_mass/charge)
    
    # Print extracted values
    print("Extracted values:", flush=True)
    print("mz:", mz, flush=True)
    print("Expected light peak:", expected_light_peak, flush=True)
    
    sequence = peptide_details['Sequence'].values[0]
    df = pd.read_csv(path + "/temp data/" + sequence + ".csv")
    mz = df.mz
    intensity = df.intensity
    #plot each peptide
    plotting.scan_plot(mz, intensity, peptide_details)
    
    
            
''' TO DO 
-test with ecoli (look for highest H)
    - implement function to call
    -clear out or color non relevant peaks
    -add more or less spectra to locat the peak of interest
-Refactor code
    -separate .py for dealing with raw files etc. (DONE)
    -clean up functions performing necessary actions
    -variable names, logics etc.
-Adjust plotting function
    -label peaks (DONE)
    -plot title
    - proline recycling expected peak
-add variables to adjust for import 
    -i.e. raw file, MQ file, which protein, etc.
-How to launch as tool?
'''
  
    
   
    