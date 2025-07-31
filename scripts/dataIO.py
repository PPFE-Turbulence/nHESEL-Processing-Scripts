# Initial imports
import sys
import os
import numpy as np

# Import paths file
pathPath = "/home/s194112"
sys.path.append(pathPath)
import paths

# Use dataProcessing script as dependency
# Append scripts folder
sys.path.append(paths.scriptPath)

# Import dataLoader script
from dumpProcessing import dumpProcessing

class dataIO:
    """
    Class used for loading, saving, storing, and ascessing data from BOUT outputs and processing.
    
    Init input is simply the name of the data folder (dataName), present in the folder with datafolders.
    
    Variables:
        dataName : str
            Name of the data folder to access
        processingObj : dumpProcessing type Object
            Object used to handle reading from dump files present in the data folder
        dataDict : dictionary
            Dictionary used to hold all data present in the dataIO object, either loaded, or added.
    
    Methods:
        __init__(dataName)
            Initialisation method. Defines the variables described above, which includes initialising a dumpProcessing object.
        load(varNames, minIndex, maxIndex)
            Used to load data for a given list of names for the given indices into memory (dataDict).
            Note this will overwrite whatever data was already present for the given variables.
        save(varNames)
            Saves the data of the given variables/keys into the processing folder, creating it if it does not already exist.
        returnKeys(varNames)
            Simply returns the data for the given list of keys for dataDict, in the order they are input in.
        addToKeys(varNames, varData)
            Adds the data of the given varNames and for the given varData to already present data, by *appending it along the first axis*.
            If no data is present, it simply adds it to the dictionary for the given keys.
                NOTE: APPENDS ALONG THE FIRST AXIS, MEANT TO BE USED FOR DATA ON FORMAT [t,---]
        overwriteKeys(varNames, varData)
            Overwrites whatever data is given for the given keys in varNames in the dictionary. Like adding data to an empty dictionary, except you remove the already present data first.
        clear_key(key)
            Method used to delete the data for a given key
        clear_data()
            Simple method used to delete all loaded data.
        
    """
    
    def __init__(self, dataName):
        """
        Init method for dataIO. Requires the name of a data folder
        present in the datafolder folder, defined from the paths file. Utilises the dumpProcessing object
        to load data and denormalise it.

        Parameters
        ----------
        dataName : str
            Name of data folder to load from

        Returns
        -------
        None.

        """
        self.dataName = dataName
        self.processingObj = dumpProcessing(dataName)
        self.dataDict = {}
        
        return
        
    def load(self, varNames, minIndex, maxIndex, xminIndex, xmaxIndex, explicit=False):
        """
        Method used to load a portion of a time-series of BOUT++ data,
        from the timestep minIndex to (and including) the timestep maxIndex.
        
        The variables to load (given in a list or general iterable), are saved using
        keys in a dictionary.

        Parameters
        ----------
        varNames : iterable with strings
            Iterable (list, array, etc.) containing keys (string) for loading data
        minIndex : int
            Lowest timestep index to load from dump data
        maxIndex : int
            Highest timestep index to load from dump data
        xminIndex : int
            Lowest x-value index to load from dump data
        xmaxIndex : int
            Highest x-value index to load from dump data

        Returns
        -------
        None.

        """
        # Add or overwrite data in dataDict iteratively
        for var in varNames:
            if explicit == True:
                print(f"Loading {var}, for x in [{xminIndex},{xmaxIndex}], t in [{minIndex},{maxIndex}]")
            data = self.processingObj.readDumps(varNames = [var], saveData = False, returnData = True, denormalise = True, tind = [minIndex, maxIndex], xind = [xminIndex, xmaxIndex], info=False)[0]
            self.dataDict[var] = data
        
        return
    
    def save(self, varNames, overwrite = False):
        """
        Method used to save the data for given keys in their current state to a numpy file, for later loading.
        Data can also be accessed directly from the self.dataDict instead of saving it.

        Parameters
        ----------
        varNames : iterable with strings
            Iterable (list, array, etc.) containing keys (string) for saving data
        overwrite : bool
            Whether files with same name should be overwritten

        Returns
        -------
        None.

        """
        # First check if the processed data folder exists, and create it if it does not exist
        processedPath = f"{paths.processedPath}/{self.dataName}"
        if not os.path.isdir(processedPath):
            os.mkdir(processedPath)
        
        # Save variables from dataDict using the given keys to numpy arrays
        for var in varNames:
            savePath = f"{paths.processedPath}/{self.dataName}/{var}.npy"
            if overwrite == False and os.path.isfile(savePath) == True:
                print(f"The following path already exists, and overwrite is set to false:\n{savePath}")
                pass
            else:
                np.save(file=savePath,arr=self.dataDict[var])
                
        
        return
    
    def returnKeys(self, varNames):
        """
        Method used to return the data for the given keys, stored in teh dataDict dictionary.

        Parameters
        ----------
        varNames : iterable with strings
            Iterable (list, array, etc.) containing keys (string) for which variables to return data from

        Returns
        -------
        varData : list of data
            List containing the data requested, in the order given in the input

        """
        
        varData = []
        for var in varNames:          
            data = self.dataDict[var]
            varData.append(data)
        
        return varData
    
    def addToKeys(self, varNames, varData):
        """
        Adds the given data in varData to the keys given in varNames, without overwriting it.
        The data is assumed to be on the format [t,---], and is therefore appended along the first axis.

        Parameters
        ----------
        varNames : iterable with strings
            Iterable (list, array, etc.) containing keys (string) indicating which variables to add data to
        varData : iterable containing bout data arrays
            Iterable (list, array, etc.) containing bout data arrays (numpy equivalent), which should be added to the data arrays

        Returns
        -------
        None.

        """
        IOtest = False
        
        # Go over each varname and its corresponding data, adding it to the data dictionary
        for i in range(len(varNames)):
            # Current variables
            var = varNames[i]
            newData = varData[i]
            
            # Check if the data even exist in the current dictionary
            if var not in self.dataDict: # If it doesn't, add it
                self.dataDict[var] = newData
            else: # If it does, combine it
                # Fetch current data
                oldData = self.dataDict[var]
                
                # Combine data
                if IOtest == True:
                    print(f"Shape of old data for key '{var}': {np.shape(oldData)}")
                    print(f"Shape of new data for key '{var}': {np.shape(newData)}")
                combinedData = np.append(oldData,newData,axis=0)
                
                # Overwrite old data for key
                self.dataDict[var] = combinedData
        
        return
    
    def overwriteKeys(self, varNames, varData):
        """
        Adds the given data in varData to the keys given in varNames.
        
        If the data is not present already, it is simply added to the dictionary already present.
        If the data for the given key *is* present, then it is overwritten.

        Parameters
        ----------
        varNames : iterable with strings
            Iterable (list, array, etc.) containing keys (string) indicating which variables to add/overwrite
        varData : iterable containing bout data arrays
            Iterable (list, array, etc.) containing bout data arrays (numpy equivalent) for each key in varNames

        Returns
        -------
        None.

        """
        for i in range(len(varNames)):
            # Define the current data
            var = varNames[i]
            data = varData[i]
            
            # Add it to the data dictionary, overwriting already present values
            self.dataDict[var] = data
        
        return
    
    def clear_key(self,key):
        """
        Simple method used to clear data for a given key.

        Parameters
        ----------
        key : str
            Key to delete data for.

        Returns
        -------
        None.

        """
        if key in self.dataDict:
            del self.dataDict[key]
        return
    
    def clear_data(self):
        """
        Method used to clear all loaded data from memory.

        Returns
        -------
        None.

        """
        del self.dataDict
        self.dataDict = {}
        return