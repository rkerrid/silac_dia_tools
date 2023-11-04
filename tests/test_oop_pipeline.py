# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 16:22:20 2023

@author: robbi
"""
from icecream import ic
from silac_dia_tools.pipeline.pipeline import Pipeline

if __name__ == "__main__":
    pipeline = Pipeline('G:/My Drive/Data/data/BM data/', 'filtering_parameters.json', contains_reference = False, pulse_channel="H", meta="your_meta")
   
    
    pipeline.preprocess()
    pipeline.format_channels()
    pipeline.roll_up_to_protein_level()
    pipeline.output_unnormalized()
    pipeline.output_dlfq()
    pipeline.output_href()