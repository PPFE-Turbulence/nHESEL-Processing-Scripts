
# Module import
import sys
import os
import numpy as np
from boutdata import collect
from configparser import ConfigParser
import time
from math import floor, ceil

# Import paths file
pathPath = f"{os.environ['HOME']}/{'nHESEL-Processing-Scripts'}"
sys.path.append(pathPath)
import paths

# Append scripts folder
sys.path.append(paths.scriptPath)

# Import submission script
from submissions import bsubmissions

class dumpProcessing:
    """
    Object used to process dump files from HESEL.

    The methods are:
     - readDumps(varNames, dumpPath = "Default", saveData = True, returnData = False, denormalise = True, overwrite = False, tind = False, xind = False, info = True, explicit = False)
     - get_gridvals()
     - get_times()
     - get_grids(scrapeoff = False)
     - get_lcfs()
     - get_lcfs_index(includeLCFS = False)
     - get_wall(scrapeoff = False)
     - get_option_value(varName)
     - denormalise_data(self, dataFrame, varName, explicit = False)

    Some legacy methods are present.

    """
    
    def __init__(self, dataName):
        self.dataName = dataName
        self.denormVals = {}
        return
    
    def readDumps(self, varNames, dumpPath = "Default", saveData = True, returnData = False, denormalise = True, overwrite = False, tind = False, xind = False, info = True, explicit = False):
        """
        Function to read BOUT++ dump files, merge them, and output simple numpy 3D array in format (t,x,y).
        Renormalisation to SI units is done using BOUT.inp file parameters from data folder.
    
        Parameters
        ----------
        varName : str or chr
            Name of variable to be loaded. Types are 'n' (density), 'phi' (electrostatic potential), Ti (ion temperature), and Te (electron temperature).
        dumpPath : str
            Path to folder with dump files from BOUT. Either relative (example: "./dumpFiles" for the folder 'dumpFiles' in the current working directory) or absolute (full path).
        saveData : bool
            Bool indicating whather data should be saved to processed data folder or not
        returnData : bool
            Bool indicating whether the data should be returned from the method or not. Do note that the data is pretty large (~4 GB per variable not uncommon)
        denormalise : bool
            Toggles renormalisation of variables to SI units. Should always be left on, but present for the sake of versatility.
        overwrite : bool
            Toggles whether the data already present in the targeted folder should be overwritten or not. False by default.
        tind : list containing two ints
            Lower and upper values of t to read data from
        xind : list conatining two ints
            Lower and upper values of x to read values from
        info : bool
            Bool indicating whether the info from loading data should be displayed or not. True by default.
    
        Raises
        ------
        AssertionError
            Check on varible types of varName and dumpPath.
    
        Returns
        -------
        data : list of BoutArray
            3D bout arrays containing datatype varName, in format (t,x,y).
    
        """
        # For default option of dump-path definition
        if dumpPath == "Default":
            dumpPath = f"{paths.dataPath}/{self.dataName}"

        if not isinstance(dumpPath, str):
                    raise AssertionError(f"dumpPath should be of type string. It is currently of type {type(dumpPath)}.")

        # Check tind
        if tind == False:
            tind = [0,-1] # All the data
            
        # Check xind
        if xind == False:
            xind = [0,-1]
        
        # Define empty list for variables to return
        if returnData == True:
            datas = []
        for varName in varNames:
            # Define save paths
            saveFolderPath = f"{paths.processedPath}/{self.dataName}"
            savePath = f"{saveFolderPath}/{varName}" # Used for np.save
            savePathName = f"{savePath}.npy" # Exact name for saved file

            # Check for variable type
            if not isinstance(varName, str): # Check whether varName is not of type str or a subclass of it
                raise AssertionError(f"varName '{varName}' should be of type string or character. It is currently of type {type(varName)}.")

            if explicit == True:
                print(f"Collecting data for variable: {varName}")
            # Read & collect data from dump files, for the given varName data type
            # 'n' for density for example
            # 'lil dirty fix for annoying output messages for Ti and Te
            if varName in ["Te", "Ti"]: # If this is breaking anything, uncomment from here...
                collectName = varName.lower()
            else:
                collectName = varName # ...to here, and change collectName to varName in the line below
            data = collect(collectName, path = dumpPath, tind = tind, xind = xind, xguards = False, yguards = False, info = info)
            data = data.squeeze() # One dimension is empty, as we only simulate in 2D for nHesel, so it is preferable we get rid of this dimension

            # Denormalise data
            if denormalise == True:
                data = self.denormalise_data(data, varName)

            # Add data to list of data
            if returnData == True:
                datas.append(data)
    
            # Check overwrite condition, and save to file
            if saveData == True and os.path.isfile(savePathName) == True and overwrite == False:
                print(f"Processed data for varName = '{varName}' already exists, and overwriting is set to {overwrite}. Skipping variable.\n")
            elif saveData == True:
                if os.path.isfile(savePathName) == True and overwrite == True:
                    print(f"Processed data for varName = '{varName}' already exists, and overwriting is set to {overwrite}. Overwriting data.\n")
                    
                # Save data as numpy array, either as overwriting or simply saving it
                # Print statement
                if explicit == True:
                    print(f"Saving data to: {savePathName}\n")

                # Create folder if it does not exist
                if not os.path.isdir(saveFolderPath): # First check if the processed data directory exists already
                    os.mkdir(saveFolderPath)
                
                # Save the data
                np.save(savePath, data)

        if explicit == True:
            print("Finished")
        
        # Return the data
        if returnData == True:
            return datas

    def get_gridvals(self, recalculate = False):
        """
        Method used to get the denormalised (2D) xz-grid spacing from the given data folder.

        Two dernormalised floats dx and dz are returned.
        """
        ## First check if the data is present
        fileNamedx = f"{paths.processedPath}/{self.dataName}/gridValsdx.npy"
        fileNamedz = f"{paths.processedPath}/{self.dataName}/gridValsdz.npy"
        if (recalculate == True) or ((os.path.isfile(fileNamedx) != True) or (os.path.isfile(fileNamedz) != True)):
            # If it isn't, then load from RAW data
            # First get the dx and dz variable
            dx = collect(varname="dx", xind=[0,0], xguards = False, yguards = False, info=False, path=f"{paths.dataPath}/{self.dataName}")
            dz = collect(varname="dz", xind=[0,0], xguards = False, yguards = False, info=False, path=f"{paths.dataPath}/{self.dataName}")
    
            # Then denormalize them (get them on SI units)
            dx = float(self.denormalise_data(dataFrame = dx, varName = "dx").item(0))
            dz = float(self.denormalise_data(dataFrame = dz, varName = "dz").item(0))

            # And remember to save the data for next usage
            np.save(fileNamedx, dx)
            np.save(fileNamedz, dz)

        else:
            # If it is, then load directly from the processed data
            dx = np.load(fileNamedx)
            dz = np.load(fileNamedz)

        return dx, dz
        
    def get_times(self, recalculate = False):
        """
        Returns an array with the denormalised times for all iterations along the t-axis in the given dataset.
        """
        ## First check if the data is present
        fileName = f"{paths.processedPath}/{self.dataName}/times.npy"
        if (recalculate == True) or (os.path.isfile(fileName) != True):
            # If it isn't, then load from RAW data
            # First get the t_array with all the times
            t_array = np.asarray(collect(varname="t_array", xguards = False, yguards = False, info=False, path=f"{paths.dataPath}/{self.dataName}"))
    
            # Then denormalise it
            t_array = self.denormalise_data(dataFrame = t_array, varName = "t_array")
            
            # And remember to save the data for next usage
            np.save(fileName, t_array)

        else:
            # If it is, then load directly from the processed data
            t_array = np.load(fileName)
        
        return t_array

    def get_grids(self, scrapeoff = False, recalculate = False):
        """
        Method used to get the xz-grid from the given data folder.

        scrape = True will return a grid corresponding to everything outside of (and including) the LCFS position

        Two 1D arrays xgrid and zgrid are returned.
        """
        ## First check if the data is present IN THE CORRECT FORMAT
        fileNamexgrid = f"{paths.processedPath}/{self.dataName}/xgridsScrapeoff{scrapeoff}.npy"
        fileNamezgrid = f"{paths.processedPath}/{self.dataName}/zgridsScrapeoff{scrapeoff}.npy"
        if (recalculate == True) or ((os.path.isfile(fileNamexgrid) != True) or (os.path.isfile(fileNamezgrid) != True)):
            # If it isn't, then load from RAW data
            # First get the size of the distance between points, in both directions
            dx, dz = self.get_gridvals()
    
            # Then define the lower bound of the grid (whether inside of outside the LCFS - LCFS is included due to flooring of int)
            if scrapeoff == True:
                # Get the LCFS index
                lowIndex = self.get_lcfs_index()
            else:
                lowIndex = 0
    
            # First get a single frame without xguards and zguards
            frame = collect(varname="n", path = f"{paths.dataPath}/{self.dataName}", xind=[lowIndex,-1], tind = [0,0], xguards = False, yguards = False, info = False)
    
            # Then extract the lengths from the dimensions of this frame
            nx = int(np.shape(frame)[1])
            nz = int(np.shape(frame)[3])
    
            # Then generate the gridvalues in the chosen direction (numerical coordinate of point from origo)
            xgrid = np.arange(0, dx*nx+dx, dx)
            zgrid = np.arange(0, dz*nz+dz, dz)

            # And save for next time
            # Remember to save for the scrapeoff type
            np.save(fileNamexgrid, xgrid)
            np.save(fileNamezgrid, zgrid)

        else:
            # If it is, then load it from the processed data
            # Remember to check for type (scrapeoff) when loading
            xgrid = np.load(fileNamexgrid)
            zgrid = np.load(fileNamezgrid)

        return xgrid, zgrid

    def get_lcfs(self, recalculate = False):
        """
        Method used to get the LCFS relative position in the system, as well as the actual position, relative to the inner part of the system.

        Returns the relative value x_lcfs and the actual value x_lcfsPos.
        """
        ## First check if the data is present IN THE CORRECT FORMAT
        fileNamexLCFS = f"{paths.processedPath}/{self.dataName}/x_lcfs.npy"
        fileNamexLCFSPOS = f"{paths.processedPath}/{self.dataName}/x_lcfsPOS.npy"
        if (recalculate == True) or ((os.path.isfile(fileNamexLCFS) != True) or (os.path.isfile(fileNamexLCFSPOS) != True)):
            # If it isn't, then load from RAW data
            # First get the relative position of the lcfs in the system
            x_lcfs = self.get_option_value("x_lcfs")
    
            # Then get the maximum position in the system, and use this to find what the relative position corresponds to
            xscale = np.max(self.get_grids()[0])
            x_lcfsPos = x_lcfs*xscale

            # Remember to save the data
            np.save(fileNamexLCFS, x_lcfs)
            np.save(fileNamexLCFSPOS, x_lcfsPos)

        else:
            # If it is, then simply load it
            x_lcfs = np.load(fileNamexLCFS)
            x_lcfsPos = np.load(fileNamexLCFSPOS)
        
        return x_lcfs, x_lcfsPos

    def get_lcfs_index(self, includeLCFS = False, recalculate = False):
        """
        Method used to return the lowest index corresponding to the LCFS (for inclusivity) in the simulation grid, along the x-direction.

        """
        ## First check if the variable had already been processed
        # Also remember to check for includeLCFS
        fileName = f"{paths.processedPath}/{self.dataName}/nLCFS_includeLCFS{includeLCFS}.npy"
        if (recalculate == True) or (os.path.isfile(fileName) != True):
            # If it hasn't, process it
            # Get the relative position of the LCFS
            x_lcfsRel = self.get_option_value("x_lcfs")
    
            # Now get a single frame without xguards and zguards
            frame = collect(varname="n", path = f"{paths.dataPath}/{self.dataName}", tind = [0,0], xguards = False, yguards = False, info = False)
    
            # Then extract the lengths from the dimensions of this frame
            nx = int(np.shape(frame)[1])
            nz = int(np.shape(frame)[3])
    
            # Then multiply the relative x position with the indices, and floor or ceil the value, depending on whether the LCFS should be included or not
            if includeLCFS == True:
                nLCFS = floor(nx*x_lcfsRel)
            elif includeLCFS == False:
                nLCFS = ceil(nx*x_lcfsRel)

            # Remember to save it for future usage
            # Also, remember to include the LCFS mode in saving
            np.save(fileName, nLCFS)

        else:
            # If it has already been processed, simply load it
            # Also remember to check for mode
            nLCFS = np.load(fileName)

        return nLCFS

    def get_wall(self, scrapeoff = False, recalculate = False):
        """
        Method used to get the wall relative position in the system, as well as the actual position, relative to the inner part of the system.

        Returns the relative/normalised value x_wall and the actual value x_wallPos.

        Scrapeoff argument determines whether the wall position should be given relative to the LCFS seperatrix or not.
        """
        ## First check if the variable had already been processed
        # Also remember to check for scrapeoff
        fileNamex = f"{paths.processedPath}/{self.dataName}/x_wall_scrapeoff{scrapeoff}.npy"
        fileNamexPos = f"{paths.processedPath}/{self.dataName}/x_wallPos_scrapeoff{scrapeoff}.npy"
        if (recalculate == True) or ((os.path.isfile(fileNamex) != True) or (os.path.isfile(fileNamexPos) != True)):
            # If it hasn't been processed, load from RAW data
            # First get the relative position of the 'wall'
            if scrapeoff == False:
                x_wall = self.get_option_value("x_wall")
            elif scrapeoff == True:
                x_wall = self.get_option_value("x_wall") - self.get_option_value("x_lcfs")
    
            # Then the x-grid max position
            xscale = np.max(self.get_grids()[0])
    
            # Finally use the max position as a scaling factor multiplied by the relative position, to get the actual position
            x_wallPos = x_wall*xscale

            # Remember to save for next time
            np.save(fileNamex, x_wall)
            np.save(fileNamexPos, x_wallPos)

        else:
            # If it has been processed before, simply load it
            # Remember to check for scrapeoff option
            x_wall = np.load(fileNamex)
            x_wallPos = np.load(fileNamexPos)
        
        return x_wall, x_wallPos

    def preprocess_data(self, recalculate = True):
        """
        Method calling the different methods for getting the wall position, lcfs index, grids, time, etc. Everything except the actual data extraction.

        This method can be called to ensure everything gets processed from the RAW data into the processed data folder, such that the RAW data is no longer nescessary.

        Note, this can also wipe the already present data as long as recalculate is set to true. Useful for when a new dataset is present.

        """
        print(f"Pre-processing dumpProcessing data for folder: {self.dataName}")
        self.get_gridvals(recalculate = recalculate)
        self.get_times(recalculate = recalculate)
        self.get_grids(scrapeoff = False, recalculate = recalculate)
        self.get_grids(scrapeoff = True, recalculate = recalculate)
        self.get_lcfs(recalculate = recalculate)
        self.get_lcfs_index(includeLCFS = False, recalculate = recalculate)
        self.get_lcfs_index(includeLCFS = True ,recalculate = recalculate)
        self.get_wall(scrapeoff = False, recalculate = recalculate)
        self.get_wall(scrapeoff = True, recalculate = recalculate)
        optionVals = {'n':'n0', 'ti':'Te0', 'te':'Te0', 'phi':'Te0', 'b':'B0', "t_array":"oci", "dx":"rhos", "dz":"rhos", "nx":"nx", "ny":"ny", "nz":"nz", "x_lcfs":"x_lcfs","x_wall":"x_wall", "rmajor":"Rmajor","rminor":"Rminor","q":"q"}
        for key in optionVals.keys():
            self.get_option_value(key, recalculate = recalculate)
        print("Finished Pre-Processing")

        return 0
        
        
    def get_option_value(self, varName, recalculate = False):
        """
        Function used to extract the corresponding denormalisation variable from the BOUT settings, read directly from the dump files, used *for normalization* purposes.

        Note that this does *not* return the actual value, but the value used to denormalise it, unless the variable name and dependency name is set to be the same in the normNames dict.

        TODO: Also simply loads the value from saved values, in case it has already been loaded once before.

        """
        # Case-sensitivity is handled here and by the boutdata library
        varNameLower = varName.lower()
        normNames = {'n':'n0', 'ti':'Te0', 'te':'Te0', 'phi':'Te0', 'b':'B0', "t_array":"oci", "dx":"rhos", "dz":"rhos", "nx":"nx", "ny":"ny", "nz":"nz", "x_lcfs":"x_lcfs","x_wall":"x_wall", "rmajor":"Rmajor","rminor":"Rminor","q":"q"} # Dictionary of variables to normalise with

        # Check for varName in above
        if varNameLower not in normNames.keys():
            raise Exception(f"'{varName}' is not an implemented variable in the dictionary. The current valid keys are: {normNames.keys()}\nTo use varName={varName} in this method, add the variable as a key to the 'sections' and 'normNames' dictionaries in this method, with the corresponding section and variable name for the BOUT.inp options file.")

        # Check if data already exists/has been extracted
        fileName = f"{paths.processedPath}/{self.dataName}/option_value_{normNames[varNameLower]}.npy"
        if (recalculate == True) or (os.path.isfile(fileName) != True):
            #print("Extracting value")
            # In case it hasn't, get data from files
            val = collect(varname=normNames[varNameLower], xguards = False, yguards = False, info=False, path=f"{paths.dataPath}/{self.dataName}")

            # And save the data for next time
            np.save(fileName, val)
        else:
            #print("Using known value")
            # In case it has, simply load it
            val = np.load(fileName)
        
        return val

    def get_option_value_OLD(self, varName):
        """
        Function used to extract a given variable option from an options file, and convert it to a float.

        Old version, no longer used, but still works.

        """        
        # Hardcoded variable names for finding quantities
        varNameLower = varName.lower() # remove case-sensitivity
        sections = {'n' : 'hesel', 'ti' : 'hesel', 'te' : 'hesel', 'phi' : 'hesel', 'b' : 'hesel'}
        normNames = {'n' : 'n0', 'ti' : 'Te0', 'te' : 'Te0', 'phi' : 'Te0', 'b' : 'Bt'}

        # Check for varName in above
        if varNameLower not in normNames.keys():
            raise Exception(f"'{varName}' is not an implemented variable in the dictionary. The current valid keys are: {normNames.keys()}\nTo use varName={varName} in this method, add the variable as a key to the 'sections' and 'normNames' dictionaries in this method, with the corresponding section and variable name for the BOUT.inp options file.")

        # Open input file
        file = open(f"{paths.dataPath}/{self.dataName}/BOUT.inp","r")

        # Write to temporary file (need section header for first variables)
        tempFile = open(f"{pathPath}/temp.inp","w+")
        tempFile.write("[diverse]\n")
        for line in file:
            tempFile.write(line)
        tempFile.close()

        # Read temporary file using configparser
        tempFile = open(f"{pathPath}/temp.inp")
        parser = ConfigParser()
        parser.read_file(tempFile)

        # Then delete the temporary file
        tempFile.close()
        os.remove(f"{pathPath}/temp.inp")

        # We are now ready to parse the config for variables
        # First get the names using keys from the predefined lists
        secName = sections[varNameLower]
        normName = normNames[varNameLower]

        # Then get the correct variable (as a string)
        normVar = parser[secName][normName]

        # Given it is a string, and we require a number, we must convert it. First, note that there may be a comment
        # as well as spaces and scientific notation included in the string above
        # Example: 1.7543e-27 # This is a comment
        # We therefore need to only get the "1.7543e-27" part of the string
        tempStr = ""
        commentStarted = False
        for c in normVar:
            # First check whether we've reached the comment
            if c == "#":
                commentStarted = True

            # Then write characters to string
            if commentStarted == False:
                if c != " ":
                    tempStr = tempStr + c

        # Finally convert to float and define normVar again
        normVar = float(tempStr)

        return normVar

    def denormalise_data(self, dataFrame, varName, explicit = False):
        """
        Used to renormalise dataframe (numpy-like array) for a given variable, using BOUT.inp file present in data folder.

        """
        # Define additional factors to add for proper SI-conversion (except for the temperature)
        # For potentials: normalisation unit is C/eV, so renormalisering is eV/C which means we must convert C to e using a factor 6,24150907446076e18
        factors = {}#{"phi":1/6.24150907446076e18}
        
        # Get the variable used to norm the data, if it hasn't been retrieved already
        if explicit == True:
            print("Getting denormalision value...")
        if varName not in self.denormVals:
            normVar = self.get_option_value(varName)
            self.denormVals[varName] = normVar # Add to dictionary of vals
        else:
            normVar = self.denormVals[varName]
        
        # Then, given the value used to normalise the data, de-normalise the data by multiplying the variable with the dataframe
        # Remember check for t_array, which is normalised not by dividing by a value, but multiplying, so we need to do the inverse operation
        if explicit == True:
            print("Denormalising data...")
        if varName == "t_array":
            dataFrame = dataFrame/normVar
        else:
            dataFrame = dataFrame*normVar

        # Add the conversion factor, if need be
        if varName in factors:
            dataFrame = dataFrame*factors[varName]
        
        return dataFrame

    def extract_all(self, overwrite = False):
        """
        Method used to extract the variables: n, Te, Ti, phi

        Old naivé method, use new methods from heselProcessing

        """
        # Define variables to extract
        extracts = ["n","Te","Ti","phi",'B']

        # Run data extraction
        self.readDumps(extracts,overwrite=overwrite)
        
        return

    def submit_extraction_job(self, overwrite = False):
        """
        Method used to submit a simple single-threaded task to the cluster, consisting of extracting the data from a single folder.

        """
        # First write a temporary script using this object, utiliting the processing template
        template = open(f"{paths.templatePath}/TEMPLATE_processing.py","r")
        tempFile = open(f"{pathPath}/tempProcessingFile.py","w+")
        for line in template:
            line = line.replace("NAMEDATA",f"{self.dataName}")
            line = line.replace("WRITEOVER",f"{overwrite}")
            tempFile.write(line)
        tempFile.close()

        # Change the file-permissions 
        os.chmod(f"{pathPath}/tempProcessingFile.py", 0o777)
        
        # Then define a bsub object using the written temporary file
        job = bsubmissions(jobName=f"Extract_{self.dataName}",scriptName=f"{pathPath}/tempProcessingFile.py")
        job.create_and_submit()

        # Then delete the temporary file
        time.sleep(2)
        os.remove("tempProcessingFile.py")

        return
