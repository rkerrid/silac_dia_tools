import os
from silac_dia_tools.pipeline.preprocessor import Preprocessor
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
        self.preprocessor = Preprocessor(self.path, self.parameter_file ,self.meta)
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
    
    def save_preprocessing(self):
        manage_directories.create_directory(self.path, 'preprocessing')
        print('Saving filtered_report.tsv')
        self.filtered_report.to_csv(f'{self.path}preprocessing/report_filtered.tsv',sep='\t')
        print('Saving silac_precursors.tsv')
        self.formatted_precursors.to_csv(f'{self.path}preprocessing/silac_precursors.tsv',sep='\t')
        print('Saving protein_ratios.csv')
        self.filtered_report.to_csv(f'{self.path}preprocessing/protein_ratios.csv',sep='\t')
       
        
    def run_href_pipeline(self):
        self._preprocess()
        self._format_channels()
        self._roll_up_to_protein_level()
        self._output_unnormalized()
       
        self._output_href()
        
    def run_dlfq_pipeline(self):
        self._preprocess()
        self._format_channels()
        self._roll_up_to_protein_level()
        self._output_unnormalized()
        
        self._output_dlfq()
        
    def generate_reports(self): 
        print('Beginning filtering report')
        filtering_report.create_report(self.filtered_report, self.contaminants, self.filtered_out, self.path, self.params)
        print('Beginning precursor report')
        precursor_report.create_report(self.formatted_precursors, self.path, self.params)
        print('Beginning protein group report')
        protein_group_report.create_report(self.protein_groups, self.path, self.params)
        file_list = [f for f in os.listdir(f'{self.path}protein intensities') if os.path.isfile(os.path.join(f'{self.path}protein intensities', f))]

        # all_intensities = [self.unnormalized_total_intensities,
        #                    self.unnormalized_nsp_intensities,
        #                    self.light_unormalized_intensities,
        #                    self.href_total_intensities,
        #                    self.href_nsp_intensities,
        #                    self.href_nsp_light,
        #                    self.dlfq_nsp_intensities,
        #                    self.dlfq_total_intensities,
        #                    self.dlfq_total_light] 
        
        protein_intensities_report.create_report(file_list, self.path, self.params)