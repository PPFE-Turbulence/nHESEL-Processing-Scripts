# Contains *absolute* paths for all relevant folders, such as input files, raw data folder, source folders etc. for use in scripts
# Keeping all paths in one file allows for easy updating in all scripts

import os
pathPath = f"{os.environ['HOME']}/{'nHESEL-Processing-Scripts'}"

# Input files
inputPath = f"{pathPath}/HESEL_input"

# Raw data
dataPath = f"{pathPath}/HESEL_data"

# Post-processed data
processedPath = f"{pathPath}/HESEL_processed"

# Logs
logPath = f"{pathPath}/HESEL_logs"

# (HESEL) Source files
sourcePath = f"{os.environ['HOME']}"

# Scripts path
scriptPath = f"{pathPath}/scripts"

# Job Scripts path
jobScriptPath = f"{pathPath}/jobScripts"

# Template scripts path
templatePath = f"{pathPath}/templates"

# Temp files folder
tempPath = f"{pathPath}/tempFiles"

# Post-Processed data folder
postPath = f"{pathPath}/postProcessed"

# Post-Processed figurer data folder
postFigsPath = f"{pathPath}/postProcessed/Figurer"
