# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 11:52:35 2023

@author: rkerrid
"""

from icecream import ic
from silac_dia_tools.pipeline_dev.pipeline_dev import Pipeline as pileline_dev
from silac_dia_tools.pipeline.pipeline import Pipeline as pileline_old
import pandas as pd 

'''attempting to format the silac channels first then filter afterwards. Filter columns to keep in this step are:
    Parameters used:
        'Lib.PG.Q.Value'
Precursor.Charge > 1
Mass.Evidence > 0.5
Global.PG.Q.Value < 0.01
Channel.Q.Value < 0.03
Translated.Q.Value < 0.03
Translated.Quality >= 0.05

additional columns may be required


'''
if __name__ == "__main__":
    
    # path = 'G:/My Drive/Data/data/testing pipeline dev/bm old/'
    # pipeline_old = pileline_old( f'{path}', 'filtering_parameters_strict.json', contains_reference = False, pulse_channel="H", meta="meta.csv")
    # pipeline_old.run_href_pipeline()
    # pipeline_old.generate_reports()
    
    # path = 'G:/My Drive/Data/data/testing pipeline dev/bm/'
    # pipeline = pileline_dev( f'{path}', 'filtering_parameters_strict.json', contains_reference = False, pulse_channel="H", meta="meta.csv")
    # pipeline.preprocess_dev() # in href mode
    # pipeline.generate_reports()
    
    # path = 'G:/My Drive/Data/data/testing pipeline dev/eif4f old/'
    # pipeline_old = pileline_old( f'{path}', 'filtering_parameters_strict.json', contains_reference = True, pulse_channel="M", meta="meta.csv")
    # pipeline_old.run_href_pipeline()
    # pipeline_old.generate_reports()
    
    path = 'G:/My Drive/Data/data/testing pipeline dev/eif4f/'
    pipeline = pileline_dev( f'{path}', 'filtering_parameters_strict.json', contains_reference = True, pulse_channel="M", meta="meta.csv")
    print('begin writeing')
    pipeline.report.to_csv(f'{path}.report_formatted.tsv', sep=',')
    print('wrote to csv')
    pipeline.preprocess_dev(f'{path}/report_formatted.tsv') # in href mode
    pipeline.generate_reports()
   
   
