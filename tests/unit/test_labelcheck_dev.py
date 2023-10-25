# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 13:59:20 2023

@author: robbi
"""
from silac_dia_tools.label_check import file_io as io
from silac_dia_tools.label_check import process_peptides as peptide
from silac_dia_tools.label_check import plotting
from silac_dia_tools.label_check import process_peptide_details
import pandas as pd
from icecream import ic

# work
path = 'C:/data/silac_dia_tools_files/data/spikein data/'
# home
path = 'G:/My Drive/Data/data/label check proof of concept/'
path = 'G:/My Drive/Data/data/label check test/'

# # path of raw files
raw_files = io.list_raw_files(path)

#import MQ and thermo raw file data
msms = io.import_file(path + "msms.txt")
evidence = io.import_file(path + "evidence.txt")
ms1_dfs = io.get_raw_data(path, raw_files, 1000)
# ic(ms1_dfs)

# get peptide details
peptides = peptide.get_most_abundant_krp_peptide(msms,evidence) 

# process details and save temp csv
process_peptide_details.process_details(peptides, ms1_dfs, path)
# import and plot data
plotting.import_ms_data_for_plotting(path, peptides)