from silac_dia_tools.pipeline.preprocessor import Preprocessor
from silac_dia_tools.pipeline.silac_formatter import SilacFormatter
from silac_dia_tools.pipeline.precursor_rollup import PrecursorRollup
from silac_dia_tools.pipeline.calculate_intensities import IntensityCalculator

class Pipeline:
    def __init__(self, path, meta=None):
        self.preprocessor = Preprocessor(path, meta)
        self.formatter = SilacFormatter(path)
        self.precursor_rollup = PrecursorRollup(path)
        self.intensity_calculator = IntensityCalculator(path)
        #protein roll up
        # calculate intensities
        self.filtered_report = None
        self.filtered_out = None
        self.contaminants = None
        
        self.formatted_report = None
        self.protein_groups = None
        
    def preprocess(self):
        self.filtered_report, self.contaminants, self.filtered_out = self.preprocessor.import_and_filter()
        
    def format_channels(self):
        self.formatted_report = self.formatter.format_silac_channels(self.filtered_report)
    
    def roll_up_to_protein_level(self):
        self.protein_groups = self.precursor_rollup.calculate_protein_level_ratios(self.formatted_report )
        
    def output_unnormalized(self):   
       self.intensity_calculator.output_unnorm(False)
       
    def output_href(self):
       self.intensity_calculator.output_href()
    
    def output_dlfq(self):
       self.intensity_calculator.output_href()