# Module import
import sys
import numpy as np
import matplotlib.pyplot as plt

# Import paths file
pathPath = f"{os.environ['HOME']}/{'nHESEL-Processing-Scripts'}"
sys.path.append(pathPath)
import paths

# Append scripts folder
sys.path.append(paths.scriptPath)

# Import script
from extraction import extraction

# Now run the extraction script
extraction(dataName = "NAMEDATA")
