import os
import dask.dataframe as dd # add to setup
import tkinter as tk # add to setup
import pandas as pd
from pandastable import Table # add to setup
from icecream import ic 

from silac_dia_tools.pipeline.preprocessor_dev import Preprocessor
from silac_dia_tools.pipeline.silac_formatter import SilacFormatter
from silac_dia_tools.pipeline.precursor_rollup import PrecursorRollup
from silac_dia_tools.pipeline.calculate_intensities import IntensityCalculator
from silac_dia_tools.pipeline.report import filtering_report, precursor_report, protein_group_report, protein_intensities_report
from silac_dia_tools.pipeline.utils import manage_directories


class Pipeline:
    def __init__(self, path, parameter_file, contains_reference = True, pulse_channel="M", meta=None):
        # pipeline variables
        self.path = path
        self.pulse_channel = pulse_channel
        self.meta = meta
        self.parameter_file = parameter_file
        self.contains_reference = contains_reference
        
        # pipeline objects
        self.preprocessor = Preprocessor(self.path, self.parameter_file, self.meta)
        self.formatter = SilacFormatter(self.path)
        self.precursor_rollup = PrecursorRollup(self.path)
        self.intensity_calculator = IntensityCalculator(self.path, self.contains_reference, self.pulse_channel)
        
        # pipeline outputs
        self.params = self.preprocessor.params
        self.filtered_report = None
        self.filtered_out = None
        self.contaminants = None
        
        self.formatted_precursors = None
        self.protein_groups = None
        
        self.unnormalized_total_intensities = None
        self.unnormalized_nsp_intensities = None
        self.light_unormalized_intensities = None
        
        self.href_total_intensities = None
        self.href_nsp_intensities = None
        self.href_nsp_light = None
        
        self.dlfq_nsp_intensities = None
        self.dlfq_total_intensities = None
        self.dlfq_total_light = None
        
        
    def _preprocess(self):
        self.filtered_report, self.contaminants, self.filtered_out = self.preprocessor.import_and_filter()
        
    def _format_channels(self):
        self.formatted_precursors = self.formatter.format_silac_channels(self.filtered_report)
    
    def _roll_up_to_protein_level(self):
        self.protein_groups = self.precursor_rollup.calculate_protein_level_ratios(self.formatted_precursors)
        
        
    def _output_unnormalized(self):   
       self.unnormalized_total_intensities, self.unnormalized_nsp_intensities, self.light_unormalized_intensities = self.intensity_calculator.output_unnorm(self.protein_groups)
       
    def _output_href(self):
       self.href_total_intensities, self.href_nsp_intensities, self.href_nsp_light = self.intensity_calculator.output_href(self.protein_groups)
       
    def _output_dlfq(self):
       self.dlfq_total_intensities, self.dlfq_nsp_intensities, self.dlfq_total_light = self.intensity_calculator.output_dlfq(self.protein_groups, self.formatted_precursors)
    
    def make_metadata(self):
        root = tk.Tk()
        app = TestApp(self.path, root)
        root.protocol("WM_DELETE_WINDOW", app.on_closing)
        app.pack(fill="both", expand=True)  # Ensure the app fills the root window
        root.mainloop()

    def save_preprocessing(self):
        manage_directories.create_directory(self.path, 'preprocessing')
        
        print('Saving silac_precursors.tsv')
        self.formatted_precursors.to_csv(f'{self.path}preprocessing/silac_precursors.tsv',sep='\t')
        print('Saving protein_ratios.csv')
        self.filtered_report.to_csv(f'{self.path}preprocessing/protein_ratios.csv',sep='\t')
       
        
    def run_href_pipeline(self):
        self._preprocess()
        self._format_channels()
        self._roll_up_to_protein_level()
        self.save_preprocessing()
        self._output_unnormalized()
       
        self._output_href()
        
    def run_dlfq_pipeline(self):
        self._preprocess()
        self._format_channels()
        self._roll_up_to_protein_level()
        self.save_preprocessing()
        self._output_unnormalized()
        
        self._output_dlfq()
        
    def generate_reports(self): 
        print('Beginning filtering report')
        filtering_report.create_report(self.filtered_report, self.contaminants, self.filtered_out, self.path, self.params)
        print('Beginning precursor report')
        precursor_report.create_report(self.formatted_precursors, self.path, self.params)
        print('Beginning protein group report')
        protein_group_report.create_report(self.protein_groups, self.path, self.params)
        print('Beginning protein intensities report')
        file_list = [f for f in os.listdir(f'{self.path}protein intensities') if os.path.isfile(os.path.join(f'{self.path}protein intensities', f))]
        protein_intensities_report.create_report(self.path, self.params)
        
        

class TestApp(tk.Frame):
    def __init__(self, path, parent=None):
        self.path = path
        self.parent = parent
        tk.Frame.__init__(self, parent)
        self.main = self.master
        self.main.geometry('600x400+200+100')
        self.main.title('Data Entry Table')
        self.df = self.create_table_data(path)
        self.table = self.create_table()

    def create_table_data(self, path):
        # Predefined 'run' list
        print('Beginning dd import')
        dtype={'Channel.Evidence.Ms1': 'float64',
                'Channel.Evidence.Ms2': 'float64',
                'Channel.L': 'float64',
                'Channel.M': 'float64',
                'Channel.Q.Value': 'float64',
                'Mass.Evidence': 'float64',
                'Ms1.Area': 'float64',
                'Ms1.Profile.Corr': 'float64',
                'Ms1.Translated': 'float64',
                'Precursor.Normalised': 'float64',
                'Precursor.Quantity': 'float64',
                'Precursor.Translated': 'float64',
                'Quantity.Quality': 'float64'}
        
        df = dd.read_csv(f'{path}report.tsv', sep='\t',dtype=dtype)
        unique_runs = df['Run'].drop_duplicates().compute()
        print(f'unique runs: {unique_runs}')
        # Create DataFrame
        df = pd.DataFrame({'Run': unique_runs, 'Sample': ['' for _ in unique_runs], 'Treatment': ['' for _ in unique_runs]})
        return df

    def create_table(self):
        pt = Table(self, dataframe=self.df, showtoolbar=True, showstatusbar=True)
        pt.show()
        return pt

    def on_closing(self):
        # This function is called when the window is closed
        # You can add code here to handle the DataFrame
        print(self.table.model.df)  # Example: Print the DataFrame
        self.main.destroy()




