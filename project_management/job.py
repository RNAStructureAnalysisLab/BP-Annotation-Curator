# AUTHOR: Kristopher Church

from typing import Optional

class Job():
    
    def __init__(self, input_paths: Optional[list], output_paths: list, parameters: Optional[list]):
        self.input_paths = input_paths
        self.output_paths = output_paths
        self.parameters = parameters
        
    #--------------------------------------------------------------------------
    # Helper Functions
    #--------------------------------------------------------------------------
    
    