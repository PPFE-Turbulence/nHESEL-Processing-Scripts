# Module import
import sys
import numpy as np
import matplotlib.pyplot as plt

# Import paths file
pathPath = os.environ["HOME"]
sys.path.append(pathPath)
import paths

# Append scripts folder
sys.path.append(paths.scriptPath)

# Import script
from extraction import extraction

# Now run the extraction script
extraction(dataName = "NAMEDATA")
