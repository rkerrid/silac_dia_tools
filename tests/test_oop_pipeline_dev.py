# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 11:52:35 2023

@author: rkerrid
"""

from icecream import ic
from silac_dia_tools.pipeline_dev.pipeline_dev import Pipeline
import pandas as pd 

'''attempting to format the silac channels first then filter afterwards. Filter columns to keep in this step are:
    Parameters used:
Precursor.Charge > 1
Mass.Evidence > 0.5
Global.PG.Q.Value < 0.01
Channel.Q.Value < 0.03
Translated.Q.Value < 0.03
Translated.Quality >= 0.05

additional columns may be required


'''
if __name__ == "__main__":
    path = 'G:/My Drive/Data/data/testing pipeline dev/'
    # pipeline = Pipeline('G:/My Drive/Data/data/eIF4F optimization/', 'filtering_parameters_strict.json', contains_reference = True, pulse_channel="M", meta="meta.csv")
    pipeline = Pipeline( f'{path}', 'filtering_parameters_strict.json', contains_reference = False, pulse_channel="H", meta="meta.csv")
   
   
