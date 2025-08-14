### Libraries imports
import sys
import warnings

### Locals imports
# Import paths file
pathPath = os.environ["HOME"]
sys.path.append(pathPath)
import paths

# Append scripts folder
sys.path.append(paths.scriptPath)

# Import script(s)
from autoplots import autoplots
from multiplot import multiplot
from densityShows import densityShows
from extraPlots import qPlotsBasic, GammaPlotsBasic, LCFSPlots
from maxMeans import maxMeansPlot

### Code
# Try-except so that it doesn't break when running for multiple datasets
try:
    # Run functions
    autoplots(dataName = "NAMEDATA")
    multiplot(dataName = "NAMEDATA")
    densityShows(dataName = "NAMEDATA")
    qPlotsBasic(dataName = "NAMEDATA")
    GammaPlotsBasic(dataName = "NAMEDATA")
    LCFSPlots(dataName = "NAMEDATA")
    maxMeansPlot(dataName = "NAMEDATA", perpOrpara = "perp")
    maxMeansPlot(dataName = "NAMEDATA", perpOrpara = "para")
except:
    dataName = "NAMEDATA"
    warnings.warn(f"following dataName cannot be plotted: {dataName}")
