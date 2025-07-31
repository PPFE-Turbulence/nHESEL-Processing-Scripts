### Libraries imports
import sys

### Locals imports
# Import paths file
pathPath = "/home/s194112"
sys.path.append(pathPath)
import paths

# Import scripts path
sys.path.append(paths.scriptPath)

# Import relevant scripts
from processFolder import processFolder

### Code
# Then run the relevant processing script
processFolder(dataName = "NAMEDATA", plain = False, zMean = True, tMean = True, tzMean = True, zSum = False, tSum = False, tzSum = False, Init = True, scrapeoff = True)
