# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 11:52:35 2023

@author: rkerrid
"""

from icecream import ic
from silac_dia_tools.pipeline_new_version.pipeline_dev import Pipeline as pileline
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
    
    path = 'G:/My Drive/Data/data/testing pipeline dev/bm whole set/new/'
    print('initialize pipleine')
    pipeline = pileline( f'{path}', 'filtering_parameters_strict.json', meta='meta.csv')
    # pipeline.make_metadata()               
    pipeline.preprocess()
    
    # pipeline.run_href_pipeline()
    # pipeline.generate_reports()
    
    # path = 'G:/My Drive/Data/data/testing pipeline dev/bm whole set/'
    # pipeline_old = pileline_old( f'{path}', 'filtering_parameters_strict.json', contains_reference = True, pulse_channel="H", meta='meta.csv')
    # pipeline_old.run_href_pipeline()
    # pipeline_old.generate_reports()
    
    # path = 'G:/My Drive/Data/data/testing pipeline dev/bm old/'
    # pipeline_old = pileline_old( f'{path}', 'filtering_parameters_strict.json', contains_reference = False, pulse_channel="M", meta="meta.csv")
    # pipeline_old.run_href_pipeline()
    # pipeline_old.generate_reports()
    
    # path = 'G:/My Drive/Data/data/testing pipeline dev/bm/'
    # pipeline = pileline_dev( f'{path}', 'filtering_parameters_strict.json', contains_reference = False, pulse_channel="M", meta="meta.csv")
    # pipeline.preprocess_dev_dlfq() # in href mode
    # pipeline.generate_reports()
    
    # path = 'G:/My Drive/Data/data/testing pipeline dev/bm old/'
    # pipeline_old = pileline_old( f'{path}', 'filtering_parameters_strict.json', contains_reference = True, pulse_channel="M", meta="meta.csv")
    # pipeline_old.run_href_pipeline()
    # pipeline_old.generate_reports()
    
    # path = 'G:/My Drive/Data/data/testing pipeline dev/bm/'
    # pipeline = pileline_dev( f'{path}', 'filtering_parameters_strict.json', contains_reference = False, pulse_channel="M", meta="meta.csv")
    # pipeline.preprocess_dev_dlfq() # in href mode
    # pipeline.generate_reports()
    
    # path = 'G:/My Drive/Data/data/testing pipeline dev/eif4f old/'
    # pipeline_old = pileline_old( f'{path}', 'filtering_parameters_strict.json', contains_reference = True, pulse_channel="M", meta="meta.csv")
    # pipeline_old.run_href_pipeline()
    # pipeline_old.generate_reports()
    
    # path = 'G:/My Drive/Data/data/testing pipeline dev/eif4f/'
    # pipeline = pileline_dev( f'{path}', 'filtering_parameters_strict.json', contains_reference = True, pulse_channel="M", meta="meta_8h.csv")
    
    # # report = pipeline.preprocessor.import_no_filter(pipeline.filter_cols)
    # # print('saving subset tsv')
    # # report.to_csv(f'{path}subset_report.csv', sep='\t')
    # # # pipeline.report = pd.read_csv(f'{path}subset_report.csv')
    # pipeline.preprocess_dev() # in href mode
    # pipeline.generate_reports()
   
   
