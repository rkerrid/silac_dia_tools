# -*- coding: utf-8 -*-
"""
Created on Sun Aug 13 14:17:47 2023

@author: rkerrid
"""




import file_io as io
import process_peptides as peptide
import pandas as pd

import plotting

##Main
#path and file variables
raw_file = 'Janice_20230919_RJK_ECsilac_12.raw'
path = "W:/RJK/label check script tests/ecoli new/" 

#import MQ and thermo raw file data
msms = io.import_file(path + "msms.txt")

msms = msms[msms["Raw_file"]== raw_file[:-4]]
evidence = io.import_file(path + "evidence.txt")


#get peptide details
pep = peptide.get_most_abundant_heavy(msms,evidence) 

#################
#get MS1 data  
ms1_df = io.get_raw_data(path, raw_file, 10000) # when accessing the same raw file, comment the import code accordingly

mz_values,intensity_values = peptide.extract_values(ms1_df, pep)
sequence = pep['Sequence'].values[0]
print(sequence)
data = {"mz": mz_values,
        "intensity": intensity_values
        }
df = pd.DataFrame(data)
df.to_csv(path + "/temp data/" + sequence + ".csv", index=False)
#############


sequence = pep["Sequence"].values[0]
retention_time = pep["Retention_time"].values[0]
charge = pep["Charge"].values[0]
mz = pep["m/z"].values[0]
# adjuest mz by precursor apex offset

AA_mass = pep["AA_mass"].values[0]
scan_number = pep["scan_number"].values[0]
expected_light_peak = mz-(AA_mass/charge)

# Print extracted values
print("Extracted values:", flush=True)
print("mz:", mz, flush=True)
print("Expected light peak:", expected_light_peak, flush=True)

sequence = pep['Sequence'].values[0]


df = pd.read_csv(path + "/temp data/" + sequence + ".csv")


mz = df.mz
intensity = df.intensity
#plot each peptide
plotting.scan_plot(mz, intensity, pep)
 
 

            
''' TO DO 
-Refactor code
    -separate .py for dealing with raw files etc. (DONE)
    -clean up functions performing necessary actions
    -variable names, logics etc.
-Adjust plotting function
    -label peaks
    -plot title
    - proline recycling expected peak
-add variables to adjust for import 
    -i.e. raw file, MQ file, which protein, etc.
-How to launch as tool?
'''
  
    
   