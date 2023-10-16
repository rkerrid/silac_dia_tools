# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 06:35:40 2023

@author: rkerrid

Example of using the  DIA SILAC pipeline without a H spike-in

for this example we will look at an metabolic pulse labelng experiment where the set of 
samples was labeled with a M pulse. For this set, a spike in wasn't added so instead we will
stack L and H or M intensities before processing with directLFQ

"""
# first set path to library and import modules (not an installable library yet)
import sys
sys.path.append('D:/Projects phd/General scripts for proteomics/SILAC DIA pipeline/')
#import modules
from silac_dia_tools.pipeline import filtering_diann_report as fdia
from silac_dia_tools.pipeline import format_silac_precursors as pdia
from silac_dia_tools.pipeline import extract_protein_level_ratios as rdia
from silac_dia_tools.pipeline import calculate_protein_intensities as idia

'''set the path of where your report.tsv is stored. New directories will be created at 
this location so creating a new folder for the entire pipeline is advised rather than 
the default diann output location'''

# assign the location of this folder to path variable
path = 'D:/data/SILAC DIA metabolic pulse labeling/with spike in/' 

# we can now import and filter this report.tsv
fdf = fdia.import_and_filter(path)

'''now to process the filtered report so that the corosponding H, M and L precursors are all in the same
row and ratios/unnormalized stacked intensities can be appended to the report  for further processing
 note: this step can take time especially for large files''' 
pdf = pdia.format_silac_channels(path)

'''now that we have a row wise report of precursors, we can calculate proteinlevel ratios and save
these ratios and protein groups to combine with protein level intensities later in the pipeline'''

proteins = rdia.calculate_protein_level_ratios(path)

''' now we have protein level ratios we can use either the raw SILAC intensities from each run as our 
M/H and L intensities or use directLFQ '''

# we can assume the spike in is always H but for dLFQ and unnorm functions we 
# need to specify whether the pulse was with H or M. The default is M but this 
# can be changed by adding the label='H' parameter to the following functions

idia.output_unnorm(path)
idia.output_dlfq(path)


