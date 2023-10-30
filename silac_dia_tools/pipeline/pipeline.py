from preprocessor import Preprocessor

class Pipeline:
    def __init__(self, path, meta=None):
        self.preprocessor = Preprocessor(path, meta)
    
    def execute(self):
        df, contams, filtered_set = self.preprocessor.import_and_filter()
        # Add other steps of your pipeline here