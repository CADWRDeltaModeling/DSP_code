# script by Lily Tomkovic to create boundary datasets for the meta-latinhypercube cases
# takes input yaml (format in progress to be somewhat forward compatible with eventual format of repo)
# produces case folders with any modified boundary inputs

import pandas as pd
from schimpy import schism_yaml
import os

import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))

# read in yaml
