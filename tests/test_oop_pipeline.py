# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 16:22:20 2023

@author: robbi
"""
from icecream import ic
from silac_dia_tools.pipeline.pipeline import Pipeline

if __name__ == "__main__":
    pipeline = Pipeline('G:/My Drive/Data/data/BM data/', 'filtering_parameters_strict.json', contains_reference = False, pulse_channel="H", meta="your_meta")
   
    
    pipeline.run_dlfq_pipeline()
    # pipeline.run_href_pipeline()
    pipeline.generate_reports()