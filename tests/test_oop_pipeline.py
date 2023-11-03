# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 16:22:20 2023

@author: robbi
"""
from icecream import ic
from silac_dia_tools.pipeline.pipeliner import Pipeline

if __name__ == "__main__":
    pipeline = Pipeline('G:/My Drive/Data/data/BM data/', meta="your_meta")
   
    
    pipeline.preprocess()
    pipeline.format_channels()
    pipeline.roll_up_to_protein_level()
    pipeline.output_unnormormalized()
    pipeline.output_dlfq(pulse_channel='H')
    pipeline.output_href()