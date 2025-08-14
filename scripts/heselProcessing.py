# Initial imports
import sys
import os
import numpy as np
from boutdata import collect

# Import paths file
pathPath = os.environ["HOME"]
sys.path.append(pathPath)
import paths

# Use dataProcessing script as dependency
# Append scripts folder
sys.path.append(paths.scriptPath)

# Import dataIO script
from dataIO import dataIO

class heselProcessing:
    """
    Class used for loading, manipulating, calculating, and handling data relevant to the HESEL model.
    Can be used for fetching data from BOUT files using the partialDataIO object, which is initialised at the start, or
    calculate individual variables, from present data, as well as a few meta-level functions.

    To minimize memory usage while loading variables, the mean and sum functions loads data in slices, and joins the data afterwards,
    greatly reducing the memory usage required for data extraction and calculations.

    Likewise care has been taken to minimize CPU usage by keeping loaded and calculated data in memory, right up until the point where they are no longer
    needed for further calculations, after which they are dropped from memory.

    The above memory-managment does not apply to the plainVariables method, as these is no simplifying loading and calculating the data directly,
    asides from loading one at a time. Do be aware of this method, and only use it when you really need it.

    Also note some documentation of the methods available, or their variables, may be outdated.

    Units given as output are:
    n [m^3], Te [eV], Ti [eV], psi [V], u [m/s], Gamma [1/m^2*s], q_e [J/m^2*s], q_i [J/m^2*s]
    
    
    Variables:
        dataName : str
            The name of the data folder to work with
        IO : dataIO object
            Internal name of dataIO object, used for loading, saving, and holding data
        dependencies : dictionary
            Dictionary used to define which variables are required to calculate other variables, such as the electric field 'E' requiring the electric potential 'phi'
        counters : dictionary
            Dictionary used to keep a count of how many times a given variable is going to be used during calculations.
            This is mostly meant to be used for a mix of memory management and processing optimisation. If a variable needs to be used
            in more than one calculation, then there is no reason to clear it from memory. It will therefore be kept in tempVars until no longer needed (its count reaches 0)
        tempVars : dictionary
            Dictionary used for storing calculated variables for future calculations, as to prevent having to calculate the same variable multiple times.
        lastMinIndexes : dictionary
            Dictionary storing which minimum indices were last used for a given variable. This is used to check whether the base-level
            data in dataIO's dataDict (data such as density, temperature, or electric potential - data without any dependencies)
            can simply be reused, or must be reloaded due to the indices having changed since last.
        lastMaxIndexes : dictionary
            Same as lastMinIndex, but for the max index for a given key instead. Both must match.
    Methods:
        __init__(dataName)
            Initialises all the above variables, and IO object.
        readyCounter(varNames)
            Resets the counters variable to all zeros, and then calls the countKeys method for each key input
        countKeys(key)
            Recursively checks how many dependent keys are required to get the key input, including the input dependent key itself (assuming it is a dependent key)
        dataFetch(varNames, minIndex = False, maxIndex = False, partialProcessing = True, reloadData = False)
            Method used to get any required data, whether that be base-data (does not require calculations, but loading using dataIO), or dependent data (requiring calculations).
            May call the calcVars method in the second case of needing to output dependent variables.
        get_dependent_variables(varNames, minIndex, maxIndex)
            Method used in dataFetch to get all dependent variables from a list of dependent variables, either by taking from tempVar storage, or calculating them anew.
        check_loaded_independent_variables(varNames, minIndex, maxIndex, reloadData)
            Method used in dataFetch to get which variables should be loaded into dataIO, due to them not being present for the current indices already. Only works on base-variables.
        addToDict(varNames, varData, overWrite = False)
            Method used to add the given data in varData to the dataIO dataDict dictionary, either by adding to already present data (assuming [t,---] format),
            or by overwriting already present data.
        E(phi)
            Method used to calculate the electric field for an electric potential matrix, assuming [t,---] format of input
        u(E,B)
            Method used to calculate the ExB drift velocity, assuming [t,---] format of inputs
        p(n,T)
            Method used to calculate the pressure, assuming [t,---] format of inputs
        Gammas(n,u)
            Method used to calculate the particle flux, assuming [t,---] format of inputs
        Gamma_perp(Gammas):
            Method used to extract the x-direction from Gammas. Simple interface calling the first index of Gammas
        Gamma_para(Gamma_perp)
            Method used to calculate the parallel flow of particles, using the perpendicular flow and conservation of mass
        q_perp(p,u)
            Method used to calculate the energy flux, assuming [t,---] format of inputs
        process_B(B, shape)
            Method used to expand the B-field from 1D [x] to 3D [t,x,z], assuming [x] format of input B-field, and [t,x,z] format of shape input
        calcVar(var, minIndex, maxIndex)
            Method used to first fetch, and then calculate and return a given *dependent* variable for a given input range.
        summedVariable(var, minIndex, maxIndex, axes)
            Method used to first calculate (using dataFetch) the given variable, and then summing it, before returning it
        generate_slices(indexRange, stepSize)
            Method used to generate and return a list of slices (a min and max index) for a given range, used to cut up a time-series into smaller parts
        averagedVariables(varNames, indexRange = False, sliceSize = 10, returnData = False, saveData = False, axes)
            Method used to get the given varName data on a summed format, using summedVariable, over the given index range, using a given stepsize
        summedVariables(varNames, indexRange = False, sliceSize = 10, returnData = False, saveData = False, axes)
            Method used to get the given varName data on a summed format, using summedVariable, over the given index range, using a given stepsize
        plainVariables(vvarNames, indexRange = False, returnData = False, saveData = False)
            Method used to get the given varName data, over the given index range, *without partial processing*
        zaveragedVariables
            Same as averagedVariables, but specifically called with axes 2 (2) to get mean over z [t,x,z,v]
        tzaveragedVariables
            Same as averagedVariables, but specifically called with axes 0 and 2 (0,2) to get mean over z and t [t,x,z,v]
        zsummedVariables
            Same as summedVariables, but specifically called with axes 2 (2) to get mean over z [t,x,z,v]
        tzsummedVariables
            Same as summedVariables, but specifically called with axes 0 and 2 (0,2) to get mean over z and t [t,x,z,v]
        
        
    
    """
    
    def __init__(self, dataName):
        self.dataName = dataName
        self.IO = dataIO(dataName) # Object used for loading, saving, and holding data
        
        # Variables used for storing information about the values which can be calculated
        self.dependencies = {"E" : ["phi"],"u" : ["E","B"],"pe" : ["n","Te"],"pi" : ["n","Ti"],"Gammas" : ["n","u"], "Gamma_perp" : ["Gammas"], "Gamma_para" : ["Gamma_perp"],"q_perpes" : ["pe","u"],"q_perpis" : ["pi","u"], "q_perpe" : ["q_perpes"], "q_perpi" : ["q_perpis"], "q_parae" : ["q_perpe"], "q_parai" : ["q_perpi"]} # Dictionary containing lists of which variables are dependencies for the given variables
        self.counters_old = {"E" : 0, "u" : 0, "pe" : 0,"pi" : 0, "Gammas" : 0, "Gamma_perp" : 0, "Gamma_para" : 0, "q_perpes" : 0,"q_perpis" : 0, "q_perpe" : 0,"q_perpi" : 0, "q_parae" : 0, "q_parai" : 0} # Dictionary containing which values have to be calculated to be accesible, and how many times they can be called before being removed from tempVars. This one is outdated, and only here for legacy purposes
        self.counters = {} # Dictionary containing which values have to be calculated to be accessible, and how many times they can be called before being removed from tempVars
        self.tempVars = {} # Dictionary to hold data temporarily, until their counters reaches 0
        
        # Define indices as empty dictionaries used for dataFetch object for first run
        self.lastMinIndexes = {}
        self.lastMaxIndexes = {}
        self.lastxMinIndexes = {}
        self.lastxMaxIndexes = {}
        
        return
    
    def readyCounter(self,varNames):
        """
        Resets and updates the counter dictionary self.counters to include the amount of times a given variable (with dependencies) is going to be calculated, given a list of varnames.

        Parameters
        ----------
        varNames : list of str
            List of strings, all present in self.dependencies. If any keys are input which are *not* part of self.dependencies, they will simply be ignored.

        Returns
        -------
        None.

        """
        # First, replace the current counters with an empty dictionary. Basically resetting the counter
        self.counters = {}
        
        # Call addCountersKeys to add all required keys to the key list, by traversing the dependencies tree recursively
        self.addCountersKeys(varNames)

        # Finally, count how many times we need to call each variable, whether that be calculating it, or simply getting it from being already calculated
        self.countKeys(varNames)

        return

    def addCountersKeys(self, keyNames):
        """
        Method which when called for a given list of keys checks the dependencies tree resursively, and adds the keys to the list of keys requiring to be calculated.

        Parameters
        ----------
        keyNames : list of str
            A list of keys, applicable to fetching or calculations

        Returns
        -------
        None.

        """

        # Go over each given key iteratively
        for key in keyNames:
            # Recursive case, if there are any dependencies, call addCounterKeys on them
            if key in self.dependencies.keys():
                self.addCountersKeys(self.dependencies[key])

            # Base-case, ensure the variable is part of the set of keys
            if key not in self.counters.keys():
                self.counters[key] = 0
                

    def countKeys(self, keyNames):
        """
        Method which when called for a given list of keys checks how many times they need to be calculated for the given variable order, in order to minimize recalculations and memory usage.

        Parameters
        ----------
        keyNames : list of str
            A list of keys, applicable to fetching or calculations

        Returns
        -------
        None.

        """
        # Check each given variable
        for key in keyNames:
            # First the recursive case: In case the count for the variable is 0, check if it has any dependent variables. If it does, call countKeys on them.
            if self.counters[key] == 0:
                if key in self.dependencies.keys():
                    self.countKeys(self.dependencies[key])

            # Now the base-case. Since this variable is needed, increase the count by 1
            self.counters[key] += 1
            
        return

    def readyCounter_old(self,varNames):
        """
        Resets and updates the counter dictionary self.counters to include the amount of times a given variable (with dependencies) is going to be calculated, given a list of varnames.

        Parameters
        ----------
        varNames : list of str
            List of strings, all present in self.dependencies. If any keys are input which are *not* part of self.dependencies, they will simply be ignored.

        Returns
        -------
        None.

        """
        # First, only run this for variables which are part of the self.dependencies dictionary
        varNames = [var for var in varNames if var in self.counters.keys()]
        
        # Now, reset the counters
        for key in self.counters_old.keys():
            self.counters_old[key] = 0
        
        # Then, for each variable, recursively check how many dependencies are present
        for var in varNames:
            self.countKeys(var)
    
    def countKeys_old(self, key):
        """
        Method which when called for a given key recursively checks self.dependencies for all instances of this key's
        dependencies, updating the self.counters variable as it does so. If the input 'key' does not have any dependencies,
        which themselves also have dependencies, then only the count for the 'key' variable goes up by one.
        
        Note, this method runs infinitely if a loop of variables is present.

        Parameters
        ----------
        key : str
            A key from the dictionary self.dependencies or self.counters.

        Returns
        -------
        None.

        """
        # First increase the count for this key by 1 [Base-case, nothing more happens than this]
        self.counters_old[key] += 1
        
        # Then get the key(s) that the key 'key' depends upon
        keys = self.dependencies[key]
        
        # Now go over each key that the key requires to be calculated
        for keyIter in keys:
            # [Recursive Case] If a key is present in the dependencies list call countKeys again for this key
            if keyIter in self.dependencies:
                self.countKeys(keyIter)
            # [Base-case] If the key isn't present, do nothing
            else:
                pass
        
        return
    
    def dataFetch(self, varNames, minIndex = None, maxIndex = None, xminIndex = None, xmaxIndex = None, partialProcessing = True, reloadData = False, testing = False):
        """
        Loads and returns the given variable names/keys in varNames.
        Note that loading a variable overwrites this data in the IO dataDict.

        Parameters
        ----------
        varNames : iterable with strings
            Iterable (list, array, etc.) containing keys (string) indicating which variables to load data from
        minIndex : int, optional (if partialProcessing == False)
            Lowest timestep to load from dump data. The default is False.
        maxIndex : int, optional (if partialProcessing == False)
            Highest timestep to load from dump data. The default is False.
        partialProcessing : bool, optional
            Bool indicating whether to load a segment of the data, or the whole dataset. True indicates loading a segment of the time-series, false all time. The default is True.
        reloadData : bool, optional
            Bool indicating whether to reload all data anew, or not. In case the statement is False, only data not already present will be loaded, or if the interval has changed from last time this function was run.

        Raises
        ------
        Exception
            A check is made whether minIndex and maxIndex has been explicitly defined given partialProcessing is True.

        Returns
        -------
        varData : list of denormalised BOUT arrays
            A list consisting of denormalised BOUT data, for the given keys in varNames.

        """
        if testing == True:
            print(f"\n// Fetching '{varNames}' //")
        ## First check partial processing and indices
        # Time min and max
        
        if partialProcessing == True and (minIndex == None or maxIndex == None):
            raise Exception(f"partialProcessing is set to {partialProcessing}, while minIndex and maxIndex are set to {minIndex} and {maxIndex}. Remember to define your indices, or turn off partial processing.")
        if partialProcessing == False:
            minIndex = 0
            maxIndex = -1
        
        # xmax and xmin
        if xminIndex == None or xmaxIndex == None:
            xminIndex = 0
            xmaxIndex = -1
        
        ## Now load all data with dependencies
        varDependencies = [var for var in varNames if var in self.dependencies]
        if testing == True:
            print(f"\nCurrent varNames (to be fetched): {varNames}")
            print(f"Current varDependencies: {varDependencies}")
            print(f"Current keys in IO memory (before load calls): {self.IO.dataDict.keys()}")
            print(f"Current keys in tempVar memory (before load calls): {self.tempVars.keys()}")
        calculatedVars = self.get_dependent_variables(varNames = varDependencies, minIndex = minIndex, maxIndex = maxIndex, xminIndex = xminIndex, xmaxIndex = xmaxIndex, testing = testing)
        
        ## Do a check for which independent data actually needs to be loaded anew, and which does not
        varNoDependencies = [var for var in varNames if var not in self.dependencies]
        loadNames = self.check_loaded_independent_variables(varNames = varNoDependencies, minIndex = minIndex, maxIndex = maxIndex, xminIndex = xminIndex, xmaxIndex = xmaxIndex, reloadData = reloadData, testing = testing)
        
        # Print statements for testing purposes
        if testing == True:
            print(f"Current varNoDependencies: {varNoDependencies}")
            print(f"Current loadNames (independent data which must be loaded to memory): {loadNames}")
            print(f"Current keys in IO memory (before independent-variables load call - may be called by get_dependent_variables): {self.IO.dataDict.keys()}")
        
        ## Now load and return data. Start by loading the data that requires it
        # Actually load the data for the given variables to the IO object
        self.IO.load(varNames = loadNames, minIndex = minIndex, maxIndex = maxIndex, xminIndex = xminIndex, xmaxIndex = xmaxIndex)
        
        # And read all the (independent/base) variables to a list, to be returned from this method
        varData = self.IO.returnKeys(varNoDependencies)

        if testing == True:
            print(f"Current keys in IO memory (after load): {self.IO.dataDict.keys()}")
            print(f"Current keys in tempVar memory (after load): {self.tempVars.keys()}")

        # After this, check the counters for all vars without dependencies. The counters must all be reduced by one, and if a counter reaches 0, remove the variable from self.IO (we don't need it anymore)
        for var in varNoDependencies:
            self.counters[var] -= 1
            if self.counters[var] <= 0:
                self.IO.clear_key(var)
        
        # Add the calculated variables to the varData list, while ensuring the resulting order corresponds to the input order of variables
        for i in range(len(varNames)):
            var = varNames[i]
            if var in self.dependencies: # So if it is one we didn't just load using IO
                calcedData = calculatedVars.pop(0)
                varData.insert(i,calcedData)

        if testing == True:
            print(f"Current keys in IO memory (after base-variable counter checks): {self.IO.dataDict.keys()}")
            print(f"Current keys in tempVar memory (after base-variable counter checks): {self.tempVars.keys()}")
        
        return varData
    
    def get_dependent_variables(self, varNames, minIndex, maxIndex, xminIndex, xmaxIndex, testing = False):
        """
        Method used to calculate and return dependent variables given in varNames. Note that the list of variables must only contain dependent variables.
        This method is utilised in the dataFetch method.

        Parameters
        ----------
        varNames : list of str
            List containing dependent variables keys to calculate and return.
        minIndex : int
            Lowest time-index to handle.
        maxIndex : int
            Highest time-index to handle.
        minIndex : int
            Lowest x-index to handle.
        maxIndex : int
            Highest x-index to handle.

        Raises
        ------
        Exception
            Exception raised if a base-variable is passed into this method.

        Returns
        -------
        calculatedVars : list of calculated data
            Data returned in order of keys input in varNames.

        """

        
        if testing == True:
            print(f"\n--- Current dependent variables list passed to get_dependent_variables: {varNames} ---")
        
        # First check if the variables are actually dependent variables, or base-variables. We only treat dependent variables in this method
        for var in varNames:
            if var not in self.dependencies:
                raise Exception(f"The variable {var} has been passed to the get_dependent_variables method, but this method only handles dependent-variables (variables with dependencies).")
        
        calculatedVars = []
        # Now check if the variable requires calculation
        for var in varNames:
            # First, check if the variable has *already* been calculated, and is therefore present in tempVars
            if var in self.tempVars:
                if testing == True:
                    print(f"Variable '{var}' is present in tempVars. No further calcution is required")
                
                # Add it to the list of calculated vars if so
                calculatedVars.append(self.tempVars[var])
                
                # Decrease the amount of uses by one
                self.counters[var] -= 1
                
                # Remove it if the amount of uses are now 0 or below
                if self.counters[var] <= 0:
                    if testing == True:
                        print(f"Variable '{var}' counter has reached 0 or below, and will be deleted")
                        print(f"tempVars before deletion of variables '{var}': {self.tempVars.keys()}")
                    
                    del self.tempVars[var]
                    
                    if testing == True:
                        print(f"tempVars after deletion of variables '{var}': {self.tempVars.keys()}")
            # If the variables has *not* been calculated, calculate it, reduce the counter by one, and add it if the counter is not at 0 for that variable (it must be used again)
            else:
                if testing == True:
                    print(f"Variable '{var}' is not present in tempVars, and requires calculation")
                # First get the variable
                calcedVar = self.calcVar(var = var, minIndex = minIndex, maxIndex = maxIndex, xminIndex = xminIndex, xmaxIndex = xmaxIndex, testing = testing)
                
                # Then add it to the list of calculated variables
                calculatedVars.append(calcedVar)
                
                # Reduce the counter by one
                self.counters[var] -= 1
                
                # If the counter is *not* at 0, add the calculated data to the list of temporary data
                if testing == True:
                    print(f"Current count for '{var}': {self.counters[var]}")
                if self.counters[var] > 0:
                    if testing == True:
                        print(f"Variable '{var}' counter has yet to reach 0, adding to tempVars")
                    self.tempVars[var] = calcedVar

        if testing == True:
            print(f"--- Returning Dependent Vals For '{varNames}' ---")
        
        return calculatedVars
    
    def check_loaded_independent_variables(self, varNames, minIndex, maxIndex, xminIndex, xmaxIndex, reloadData, testing = False):
        """
        Method used to check what data is already loaded and present in the dataIO object attached to this object, and usuably (same indices),
        as to prevent having to reload data already loaded. Unless reloadData is set to true, then everything gets reloaded.

        Parameters
        ----------
        varNames : list of str
            List containing dependent variables keys to calculate and return.
        minIndex : int
            Lowest time-index to handle.
        maxIndex : int
            Highest time-index to handle.
        reloadData : bool
            Bool determining whether all data should be reloaded, or only new data.

        Raises
        ------
        Exception
            Exception raised if a dependent variable is passed into this method.

        Returns
        -------
        loadNames : list of str
            List of keys to load in dataIO object.

        """
        
        ### First check if the variables are actually base-variables, or dependent variables. We only treat base-variables in this method
        for var in varNames:
            if var in self.dependencies:
                raise Exception(f"The variable {var} has been passed to the check_loaded_independent_variables method, but this method only handles base-variables (variables without dependencies).")

        if testing == True:
            print(f"\n<<< Given list of independent variables passed to check_loaded_independent_variables: {varNames} >>>")
        
        ### Do a check for which data actually needs to be loaded anew, and which does not
        ## If reloadData isn't true, check for all variables whether they need to be loaded or are already present
        if reloadData == False:
            loadNames = [] # List of variables which has yet to be loaded
            for var in varNames:
                if testing == True:
                    print(f"Variable '{var}' is in loaded dataDict: {var in self.IO.dataDict}")
                    if var in self.IO.dataDict:
                        print(f"t-indices for variable '{var}' are the same as before: {self.lastMinIndexes[var] == minIndex and self.lastMaxIndexes[var] == maxIndex}")
                        print(f"x-indices for variable '{var}' are the same as before: {self.lastxMinIndexes[var] == xminIndex and self.lastxMaxIndexes[var] == xmaxIndex}")
                # A bit long expression here, but essentially check if the variable exists in the dataIO (has it been loaded before), and then check if the last loaded indices were the same as now
                if var in self.IO.dataDict and (self.lastMinIndexes[var] == minIndex and self.lastMaxIndexes[var] == maxIndex) and (self.lastxMinIndexes[var] == xminIndex and self.lastxMaxIndexes[var] == xmaxIndex):
                    # If it is already present in dataDict and the indices are the same as last time, then we have already loaded this data before (it is 'base-data')
                    # Therefore, we do not need to add it to loadNames, as it is already loaded
                    if testing == True:
                        print(f"Passing {var}")
                    pass
                
                # If the variable is either not in the data dictionary, or one of the indices are different from last time, add it to variables to load
                else:
                    if testing == True:
                        print(f"Adding {var}")
                    # First append the variable name to list of new variables to calculate
                    loadNames.append(var)
                    
                    # Also remember to update the last loaded indices
                    self.lastMaxIndexes[var] = maxIndex
                    self.lastMinIndexes[var] = minIndex
                    self.lastxMaxIndexes[var] = xmaxIndex
                    self.lastxMinIndexes[var] = xminIndex
        
        ## In case reloadData is set to True, ensure all base-variables are loaded anew
        else:
            loadNames = varNames

        if testing == True:
            print(f"<<< Returning following names as loadNames: {loadNames} >>>\n")
        
        return loadNames
    
    def addToDict(self, varNames, varData, overWrite = False):
        """
        Adds or overwrites the given data for the input, in the corresponding order they are input in, for the IO object.

        Parameters
        ----------
        varNames : list of str
            Names/keys of variables which should be appended to, or overwritten.
        varData : list of data
            Data corresponding to the given names/keys of variables, to be appended or overwritten in the IO object.
        overWrite : bool, optional
            Whether to overwrite the data present, or append to it. Append if false, overwrite if true. The default is False.

        Returns
        -------
        None.

        """
        # Call IO to add the given data to the dataDict. Either by overwriting already present data, or appending it
        if overWrite == True:
            self.IO.overwriteKey(varNames = varNames, varData = varData)
        else:
            self.IO.addToKeys(varNames = varNames, varData = varData)
        
        return
    
    def E(self, phi):
        """
        Function used to calculate and return the E-field on the format [t,x,z,v], where v=[e_x,e_z].

        Output unit is in V/m

        Parameters
        ----------
        phi : Array
            Array containing electric potential on the format [t,x,z]

        Returns
        -------
        E : BOUT array
            E-field for the given indices, or partialprocessing is set to true. Otherwise E-field for whole time-series.
            Is calculated using the np.gradient function along x and z.
            E is on the format [t,x,z,v], where v=[e_x,e_z]

        """
        Etest = False
        # Edge-case fix for single-length slices, in which case the indices of phi will be [x,z] instead of [t,x,z], and axes to call needs to be downshifted by one
        if Etest == True:
            print(f"Size of phi: {np.shape(phi)}")
            
        if len(np.shape(phi)) == 3:
            if Etest == True:
                print("E does not have to be expanded")
            axes = (1,2)
        else:
            if Etest == True:
                print("E must be expanded")
            axes = (0,1)

        # Get the length along the first and second (x and z) axes
        dx, dz = self.IO.processingObj.get_gridvals()
        
        # Take the gradient along the axes (normally 1 and 2 - x and z), ignoring the t-axis. Result should be [v,t,x,z], or in single-slice case [v,x,z]
        E = np.gradient(phi,dx,dz,axis=axes) # Why is there no '-' in front? This would flip the direction of the field
        E = np.array(E)
        E = -E # Fixed here
        
        # Add single-length time-dimension in case it is not present [v,x,z] -> [v,t,x,z]
        if len(np.shape(phi)) != 3:
            if Etest == True:
                print(f"Size of E before expansion: {np.shape(E)}")
            E = np.expand_dims(E,axis=1)
            if Etest == True:
                print(f"Size of E after expansion: {np.shape(E)}")
            
        # Swap the axes to get [t,x,z,v] format, instead of [v,t,x,z]. This requires swapping 0 and 1 to get [t,v,x,z], 1 and 2 to get [t,x,v,z], and 2 and 3 to get [t,x,z,v]
        E = np.swapaxes(E, 0, 1)
        E = np.swapaxes(E, 1, 2)
        E = np.swapaxes(E, 2, 3)

        if Etest == True:
            print(f"Final shape of E: {np.shape(E)}")
        
        # Then return E
        return E
    
    def u(self, E, B):
        """
        Calculates the ExB-drift given a field E and B, on the format [t,x,z,v] for E and [t,x,z] for B.

        Output unit is in m/s.

        Parameters
        ----------
        E : Array
            E-field given on the format [t,x,z,v], where v=[e_x,e_y].
        B : Array
            B-field given on the format [t,x,z]. Note the B-field from nHesel may be 1D, and require an expansion to 3D as [x] -> [t,x,z].

        Returns
        -------
        u : Array
            Array of drift velocities for the given data on the format [t,x,z,v], where v=[e_x,e_z].

        """
        # First calculate the ExB product by using that B=<0,By,0> meaning ux = -E_z*B_y, uz = E_x*B_y (NOTE: INTERNAL Y AND Z, THESE ARE OPPOSITE EXTERNAL Y AND Z)
        # Remember E is on format [t,x,z,v], where v = [ex,ez]
        # First ux
        B2 = B**2 # Being done once for computational efficiency
        ux = E[:,:,:,1]*B/B2 # E_z*B_y (was -E_z*B_y)
        uz = -1*E[:,:,:,0]*B/B2 # -E_x*B_y (was E_x*B_y)
        
        # Combine them [v,t,x,z]
        u = np.array([ux,uz])
        
        # Swap axes, 0 and 1 [t,v,x,z], followed by 1 and 2 [t,x,v,z] and then 2 and 3 [t,x,z,v]
        u = np.swapaxes(u, 0, 1)
        u = np.swapaxes(u, 1, 2)
        u = np.swapaxes(u, 2, 3)
        
        # Then return u
        return u
    
    def p(self, n, T):
        """
        Calculates and returns the pressure in the system given a density 'n' and temperature 'T'.

        Output unit is in Pa.

        Parameters
        ----------
        n : Array
            Array containing particle densities given on the format [t,x,z].
        T : Array
            Array containing particle temperatures on the format [t,x,z].

        Returns
        -------
        p : Array
            Array containing particle pressure, given on the format [t,x,z]

        """
        # Calculate and return the pressure using the ideal gas law (so assuming the plasma is an ideal gas)
        p = n*T

        # Remember to add a conversion factor, as T is given in eV, and Pa = J/m^3, so Pa = 1.602176634e-19 eV/m^3
        p = p*1.602176634e-19

        # Aaaaand just in case we were given only a single-index slice to work with,
        # we only have 2 dimensions [x,z] where we would expect 3 [t,x,z], so we must add back the time-dimension
        if len(np.shape(p)) == 2:
            p = np.expand_dims(p,axis=0)
        
        return p
    
    def Gammas(self, n, u):
        """
        Calculates the particle flux along x and z given a density 'n' and velocity 'u'.

        Output unit is in 1/m^2*s

        Parameters
        ----------
        n : Array
            Array containing particle densities given on the format [t,x,z].
        u : Array
            Array containing particle velocities on the format [t,x,z,v], where v=[e_x,e_z].

        Returns
        -------
        Gammas : Array
            Array containing particle fluxes given on the format [t,x,z,v], where v=[e_x,e_z].

        """
        
        # First calculate the two Gamma_perp components
        Gamma_x = n * u[:,:,:,0] # [t,x,z]
        Gamma_z = n * u[:,:,:,1] # [t,x,z]
        
        # Then stack them together using a list [v,t,x,z]
        Gammas = np.array([Gamma_x, Gamma_z])
        
        # Then swap axes 0 and 1 [t,v,x,z], followed by 1 and 2 [t,x,v,z] and then 2 and 3 [t,x,z,v]
        Gammas = np.swapaxes(Gammas, 0, 1)
        Gammas = np.swapaxes(Gammas, 1, 2)
        Gammas = np.swapaxes(Gammas, 2, 3)
        
        # Finally return the calculated Gamma_perp values
        return Gammas

    def Gamma_perp(self, Gammas):
        """
        Interface to extract the x/perpendicular-perpendicular direction from the two Gammas calculated in the 'Gammas* method

        Output unit is in 1/m^2*s

        Parameters
        ----------
        Gammas : Array
            Array containing particle fluxes given on the format [t,x,z,v], where v=[e_x,e_z].

        Returns
        -------
        Gammas_perp : Array
            Array containing particle fluxes given on the format [t,x,z]

        """
        # Simply extract the x-direction (the perpendicular direction - indices for Gamma_perp are [t,x,z,v], v=[e_x,e_z] and we want [t,x,z] for the x-direction)
        Gamma_perp = Gammas[...,0]

        return Gamma_perp
        

    def Gamma_para(self, Gamma_perp):
        """
        Method to calculate the parallel (to field lines - y-direction aka out-of-plane) flow of particles, using the perpendicular flow and conservation of particles

        Output unit is in 1/m^2*s

        Parameters
        ----------
        Gamma_perp : Array
            Array containing particle fluxes given on the format [t,x,z]

        Returns
        -------
        Gammas_para : Array
            Array containing particle fluxes given on the format [t,x,z]

        """
        ## First the Gamma_perp derivative in x
        # Get the length along the first axis (x)
        dx = self.IO.processingObj.get_gridvals()[0]
        
        # Calculate the x-derivative (first axis - not zeroth - indices for Gamma_perp are currently [t,x,z])
        xDeriv_Gamma_perp = np.gradient(Gamma_perp,dx,axis=1)

        ## Now calculate the ballooning length
        # First get the major radius, and safety factor
        Rmajor = self.IO.processingObj.get_option_value("Rmajor")
        q = self.IO.processingObj.get_option_value("q")

        # Then find the ballooning length
        L_b = Rmajor*q

        # Finally multiply them together to get the parallel Gamma (particle flux along field lines - out of plane)
        Gamma_para = -L_b*xDeriv_Gamma_perp

        return Gamma_para
    
    def q_perps(self, p, u):
        """
        Calculates and returns the heat/energy fluxes in the system, given a density 'n' and pressure 'p'.

        Output unit is J/m^2*s

        Parameters
        ----------
        n : Array
            Array containing particle densities given on the format [t,x,z].
        p : Array
            Array containing particle pressures on the format [t,x,z].

        Returns
        -------
        q : Array
            Array containing perpendicular heat/energy flux, given on the format [t,x,z,v], v = [e_x, e_z]

        """        
        # First calculate the two Q (energy flux) components
        q_x = 3/2 * p * u[:,:,:,0]
        q_z = 3/2 * p * u[:,:,:,1]
        
        # Then stack them together using a list [v,t,x,z]
        q = np.array([q_x, q_z])
        
        # Then swap axes 0 and 1 [t,v,x,z], followed by 1 and 2 [t,x,v,z] and then 2 and 3 [t,x,z,v]
        q = np.swapaxes(q, 0, 1)
        q = np.swapaxes(q, 1, 2)
        q = np.swapaxes(q, 2, 3)
        
        # Finally return the values
        return q

    def q_perp(self, q_perps):
        """
        Extracts and returns the perpendicular (along x - outward) heat/energy flux in the system, given a density 'n' and pressure 'p'.

        Output unit is J/m^2*s

        Parameters
        ----------
        q_perps : Array
            Array containing perpendicular heat/energy fluxes, given on the format [t,x,z,v], v=[e_x,e_y]

        Returns
        -------
        q_perp : Array
            Array containing perpendicular heat/energy flux, given on the format [t,x,z]

        """        
        # Extract the first variable (e_x), and return it
        q_perp = q_perps[...,0]
        
        return q_perp

    def q_para(self, q_perp):
        """
        Calculates and returns the heat/energy flux in the system, given a density 'n' and pressure 'p'.

        Output unit is J/m^2*s

        Parameters
        ----------
        q_perp : Array
            Array containing perpendicular heat/energy flux, given on the format [t,x,z,v], v=[e_x,e_z]

        Returns
        -------
        q : Array
            Array containing parallel heat/energy flux, given on the format [t,x,z]

        """

        ## First extract the x
        
        ## First the q_perp derivative in x
        # Get the length along the first axis (x)
        dx = self.IO.processingObj.get_gridvals()[0]
        
        # Calculate the x-derivative (first axis - not zeroth - indices for q_perp are [t,x,z])
        xDeriv_q_perp = np.gradient(q_perp,dx,axis=1)

        ### THIS IS OLD THIS IS WRONG ###
        ## Now calculate the ballooning length
        # First get the major radius, and safety factor
        #Rmajor = self.IO.processingObj.get_option_value("Rmajor")
        #q = self.IO.processingObj.get_option_value("q")

        # Then find the ballooning length
        #L_b = Rmajor*q

        # Finally multiply them together to get the parallel q
        #q_para = -L_b*xDeriv_q_perp
        ### THIS IS OLD THIS IS WRONG ###

        # Then get the minor radius
        a = self.IO.processingObj.get_option_value("Rminor")

        # Finally multiply things together to get the parallel q
        q_para = -a*xDeriv_q_perp
        
        return q_para
    
    def process_B(self, B, shape):
        """
        Method to flip the direction of the B-field, and expand it from simply [x] to [t,x,z].
        This requires the input of the final shape B should have.

        Parameters
        ----------
        B : 1D Array
            B-field strength given for each value of 'x' in the system, so on the format [x]
        shape : tuple
            Tuple containing the final shape of B, on the form [t,x,z]

        Returns
        -------
        B_2D : 3D Array
            B-field expanded into the t- and z-direction, resulting in the format [t,x,z]
            Note that the B-field is assumed to be invarient in strength along t and z

        """
        ## First flip the direction of B from e_y to -e_y
        B = -1*B
        
        ## Then expand B to 3 dimensions (one swap could be avoided by doing the t-axis first, then the z-axis)
        # First repeat data along new axis (z - poloidal direction)
        B_2D = np.expand_dims(B,axis=1)
        B_2D = np.repeat(B_2D,shape[2],axis=1) # [x,z]
        
        # Then along second axis (t - time)
        B_2D = np.expand_dims(B_2D,axis=2)
        B_2D = np.repeat(B_2D,shape[0],axis=2) # [x,z,t]
        
        # We now have a B_2D on the format [x,z,t], but we need [t,x,z], so we swap axis 0 and 2 [t, z, x], then 1 and 2 [t, x, z]
        B_2D = np.swapaxes(B_2D, 0, 2)
        B_2D = np.swapaxes(B_2D, 1, 2) # Remember, this is along the y-direction
        
        return B_2D
    
    def calcVar(self, var, minIndex, maxIndex, xminIndex, xmaxIndex, testing = False):
        """
        Method used to collect required prerequisites for calculating a given variable, and then calculate them. Returns the given variable.

        """

        if testing == True:
            print(f"Calling calcVar for '{var}'")
        
        if var not in self.dependencies:
            applicableVars = [var for var in self.dependencies.keys()]
            raise Exception(f"'{var}' is not an applicable variable for calculation. Choose any from the following list instead:\n{applicableVars}")
        
        # Get required data to calculate the given variable, and calculate it
        if var == "E":
            # Load data
            data = self.dataFetch(["phi"], minIndex = minIndex, maxIndex = maxIndex, xminIndex = xminIndex, xmaxIndex = xmaxIndex, testing = testing)
            phi = data[0]
            
            # Get E
            dataCalc = self.E(phi)
        
        elif var == "u":
            # Load data
            data = self.dataFetch(["E", "B"], minIndex = minIndex, maxIndex = maxIndex, xminIndex = xminIndex, xmaxIndex = xmaxIndex, testing = testing)
            E, B = data[0], data[1]
            
            # Expand B
            B = self.process_B(B, shape=np.shape(E))
            
            # Get u
            dataCalc = self.u(E, B)
        
        elif var == "pe" or var == "pi":
            # Load data, for either case
            if var == "pe":
                data = self.dataFetch(["n","Te"], minIndex = minIndex, maxIndex = maxIndex, xminIndex = xminIndex, xmaxIndex = xmaxIndex, testing = testing)
            elif var == "pi":
                data = self.dataFetch(["n","Ti"], minIndex = minIndex, maxIndex = maxIndex, xminIndex = xminIndex, xmaxIndex = xmaxIndex, testing = testing)
            n, T = data[0], data[1]
            
            # Get p
            dataCalc = self.p(n, T)
        
        elif var == "Gammas":
            # Load data
            data = self.dataFetch(["n","u"], minIndex = minIndex, maxIndex = maxIndex, xminIndex = xminIndex, xmaxIndex = xmaxIndex, testing = testing)
            n, u = data[0], data[1]
            
            # Get Gammas
            dataCalc = self.Gammas(n, u)

        elif var == "Gamma_perp":
            # Load data
            data = self.dataFetch(["Gammas"], minIndex = minIndex, maxIndex = maxIndex, xminIndex = xminIndex, xmaxIndex = xmaxIndex, testing = testing)
            Gammas = data[0]

            # Get Gamma_perp
            dataCalc = self.Gamma_perp(Gammas)

        elif var == "Gamma_para":
            # Load data
            data = self.dataFetch(["Gamma_perp"], minIndex = minIndex, maxIndex = maxIndex, xminIndex = xminIndex, xmaxIndex = xmaxIndex, testing = testing)
            Gamma_perp = data[0]

            # Get Gamma_para
            dataCalc = self.Gamma_para(Gamma_perp)

        elif var == "q_perpes" or var == "q_perpis":
            # Load data, for either case
            if var == "q_perpes":
                data = self.dataFetch(["pe","u"], minIndex = minIndex, maxIndex = maxIndex, xminIndex = xminIndex, xmaxIndex = xmaxIndex, testing = testing)
            elif var == "q_perpis":
                data = self.dataFetch(["pi","u"], minIndex = minIndex, maxIndex = maxIndex, xminIndex = xminIndex, xmaxIndex = xmaxIndex, testing = testing)
            p, u = data[0], data[1]
            
            # Get q perpendicular
            dataCalc = self.q_perps(p, u)
        
        elif var == "q_perpe" or var == "q_perpi":
            # Load data, for either case
            if var == "q_perpe":
                data = self.dataFetch(["q_perpes"], minIndex = minIndex, maxIndex = maxIndex, xminIndex = xminIndex, xmaxIndex = xmaxIndex, testing = testing)
            elif var == "q_perpi":
                data = self.dataFetch(["q_perpis"], minIndex = minIndex, maxIndex = maxIndex, xminIndex = xminIndex, xmaxIndex = xmaxIndex, testing = testing)
            q_perps = data[0]
            
            # Get q perpendicular
            dataCalc = self.q_perp(q_perps)

        elif var == "q_parae" or var == "q_parai":
            # Load data, for either case
            # Remember the data is returned as a list, meaning we have to take the first element of it
            if var == "q_parae":
                data = self.dataFetch(["q_perpe"], minIndex = minIndex, maxIndex = maxIndex, xminIndex = xminIndex, xmaxIndex = xmaxIndex, testing = testing)
            elif var == "q_parai":
                data = self.dataFetch(["q_perpi"], minIndex = minIndex, maxIndex = maxIndex, xminIndex = xminIndex, xmaxIndex = xmaxIndex, testing = testing)
            q_perp = data[0]
            
            # Get q parallel
            dataCalc = self.q_para(q_perp)
        
        return dataCalc
    
    def averagedVariable(self, var, minIndex, maxIndex, xminIndex, xmaxIndex, axes = None, testing = False, independent = True, sumMode = False):
        """
        Method used for calculating or fething and returning a given variable, summed along the given axes.
        The base-assumption is that the data, whether calculated or loaded, is on the format [t,x,z,...],
        and the default axis to sum along is therefore z.

        Parameters
        ----------
        var : str
            Variable key to take sum of.
        minIndex : int
            Minimum time-index to get values for.
        maxIndex : int
            Maximum time-index to get values for.
        axes : tuple, optional
            Tuple containing the axes the sum along, with the assumption of the format [t,x,z,v]. The default is (2).

        Returns
        -------
        varAveraged : BOUT/numpy array
            Array of averaged data, along the given axes. Default results in data on format [t,x,...], potentially [t,x,v].

        """
        # Axes check
        if axes == None:
            raise Exception(f"'axes' is current set to {axes}. Please define the axes to mean over.")
        
        # Check if the key is present in self.counters, if independent is True. If not, add it.
        if independent == True and var not in self.counters:
            self.readyCounter(varNames = [var])

            if testing == True:
                print(f"'independent' is set to '{independent}', and the variable '{var}' was not present in self.counters. It has been added by calling self.readyCounter.")
        
        # Get the given variable from dataFetch (which may call calcVar)
        varData = self.dataFetch(varNames = [var], minIndex = minIndex, maxIndex = maxIndex, xminIndex = xminIndex, xmaxIndex = xmaxIndex, testing=testing)[0]

        # Now make sure all axes are present, even if they do not contain any data
        # First add the zeroth/time dimension back in, in case the indices were the same
        #print(f"Shape of data before minIndex == maxIndex check: {np.shape(varData)}")
        if minIndex == maxIndex:
            if testing == True:
                print(f"minIndex == maxIndex, adding zeroth dimension within averagedVariable to array: {var}")
            if var not in self.dependencies.keys(): # Temporary permanent solution. Some variables (the ones where they get calculated - dependet ones) are fine without the dimension addition, the ones that just get pulled directly from the data needs to get the dimension re-added
                varData = np.expand_dims(varData,0)

        #print(f"Shape of data before xminIndex == xmaxIndex check: {np.shape(varData)}")
        # Then add back in the x-dimension in case those indices were the same
        if xminIndex == xmaxIndex:
            if testing == True:
                print(f"xminIndex == xmaxIndex, adding first dimension within averagedVariable to array: {var}")
            varData = np.expand_dims(varData,1)

        # # Make sure the fourth axis is present, even if there is nothing to it, for consistency's sake
        # if len(np.shape(varData)) == 3:
        #     varData = np.expand_dims(varData,axis=3) # Count from 0
        #print(f"Shape of data before mode check: {np.shape(varData)}")
        
        # Then take the sum over the given axes, and return this (base is the z-direction, assuming [t,x,z,v]
        if sumMode == False:
            varAveraged = np.mean(varData,axis=axes)
        elif sumMode == True:
            varAveraged = np.sum(varData,axis=axes)

        #print(f"Shape of data after mode check: {np.shape(varData)}")

        # If 0 was part of the meaned axis, remember to add in this axis again though. We are expecting an output where t exists in further code
        # Only of it is not independent
        if type(axes) == int:
            tcheck = (0==axes)
        else:
            tcheck = (0 in axes)
        if tcheck and independent == False:
            varAveraged = np.expand_dims(varAveraged,0)
        
        return varAveraged
    
    def generate_slices(self, indexRange, stepSize):
        """
        Method used to generate a list of lists of min and max indices used for slicing up a series of data.
        
        Remember, endpoints are also loaded. The WHOLE range from min to and including max is loaded, when using BOUT.

        Parameters
        ----------
        indexRange : list of two int
            List containing two ints: a minimum index, and a maximum index, in this order.
        stepSize : int
            Size of the slices.

        Raises
        ------
        Exception
            An exception is raised if the indices are out of order, aka if the min index is higher than the max index

        Returns
        -------
        slices : list of slices
            List containing lists of min and max indices, used for slicing up an interval.

        """
        # First check if the indexes makes sense
        if indexRange[0] > indexRange[1]:
            raise Exception("The first index is current higher than the second index in indexRange.")
        
        # Count the amount of indexes to iterate over
        indexMin = indexRange[0]
        indexMax = indexRange[1]
        nIter = indexMax - indexMin + 1 # The +1 is required for the endpoint
        
        # Initial List
        sliceSize = stepSize
        leftover = nIter%sliceSize
        sliceAmount = int((nIter-leftover)/sliceSize)
        slices = []
        for i in range(sliceAmount):
            slices.append([indexMin+i*sliceSize,indexMin+(i+1)*sliceSize-1])
        
        # Adding the final slice
        if leftover != 0:
            slices.append([indexMin+sliceAmount*sliceSize,indexMin+sliceAmount*sliceSize+leftover-1])
        
        return slices
        
    def plainVariables(self, varNames, tindexRange = None, xindexRange = None, returnData = False, saveData = False, testing = False, info = False):
        """
        Method used to get or save the given varName data, over the given index range,
        using a given stepsize. Note this does not work well for larger variables,
        if you choose to want them returned, instead of saved.
        
        The reason for this is that to return variables, they must be kept in memory.
        If we simply wish to save them, we can just delete them after we are finished
        (not including base-variables, only dependent variables)

        Parameters
        ----------
        varNames : list of str
            List of variables to save or return, must correspond to data within the BOUT-data,
            which is defined in the dumpProcessing denormalise method local dictionary,
            or a dependent variable within this method.
        tindexRange : list of int, optional
            List or array of two ints, one being minimum t index, the other maximum t index. False means use the whole index. The default is False.
        xindexRange : list of int, optional
            List or array of two ints, one being minimum x index, the other maximum x index. False means use the whole index. The default is False.
        returnData : bool, optional
            Whether to return the data as a list, in the order given in varNames. The default is False.
        saveData : bool, optional
            Whether to save the data to the processed folder. The default is False.

        Returns
        -------
        dataToReturn : list of data
            List of data on format [t,x,z,...], corresponding to input varNames, in the same order.

        """
        
        # Check for sanity
        if returnData == False and saveData == False:
            print("Both returnData and saveData are false. Exiting.")
            return

        # Create plain var names
        plainNames = [f"{var}_plain" for var in varNames]
        
        # First, start the data counter (doing this means data is saved if needed for more than one calculation)
        self.readyCounter(varNames = varNames)
        
        # Then, handle the indices
        # Note that since we are not doing any splitting of the data, there is no point in doing any kind of partial processing
        # First time-indices
        if tindexRange == None:
            minIndex = 0
            maxIndex = -1
        else:
            minIndex = tindexRange[0]
            maxIndex = tindexRange[1]
            
        # Then spatial-indices
        if xindexRange == None:
            xminIndex = 0
            xmaxIndex = -1
        else:
            xminIndex = xindexRange[0]
            xmaxIndex = xindexRange[1]
        
        # If the data needs to be returned, add them to the dataIO object, and DON'T delete them
        # Otherwise, just load them and save them one at a time, deleting them as soon as they are saved
        # If they are not base-variables, calculate them, and add them to the list of things to return
        if testing == True:
            print(f"Initial counters: {self.counters}")
        for i in range(len(varNames)):
            var = varNames[i]
            plainVar = plainNames[i]
            if testing == True:
                print()
                print(f"==== Getting variable: {var} =====")
                print(f"Current counters: {self.counters}")
            if info == True:
                print(f"==== Getting variable: {var} =====")
            # Get the data, and add it to IO, overwriting to make sure we know which data we are working with
            # Unless it is base-data (then we don't need to write it twice, as it is already loaded)
            data = self.dataFetch(varNames = [var], minIndex = minIndex, maxIndex = maxIndex, xminIndex = xminIndex, xmaxIndex = xmaxIndex, testing = testing)[0]
            #if var in self.dependencies:
            #    if testing == True:
            #        print(f"Adding var '{plainVar}' to memory in plain loop")
            if testing == True:
                print(f"Adding var '{plainVar}' to memory in plain loop")
            self.IO.overwriteKeys([plainVar], [data])
            
            # Save the data, and don't keep it if we don't need to return it, and it is not a base-variable (these may be needed multiple times, and are not kept in this object, but the dataIO object)
            if saveData == True:
                self.IO.save([plainVar])
            if returnData == False and var in self.dependencies:
                self.IO.clear_key(plainVar)
        
        # Then, if we need to return any data, get it from the dataIO object, clear all IO data (base-data), and return it
        if returnData == True:
            dataToReturn = self.IO.returnKeys(varNames = plainNames)
            self.IO.clear_data()
            return dataToReturn
        # Otherwise, just clear all data, so it doesn't fill up memory space anymore (should no longer be needed with code update, though probably still fine to do.)
        else:
            self.IO.clear_data()
            return
    
    def averagedVariables(self, varNames, tindexRange = None, xindexRange = None, tsliceSize = 10, returnData = False, saveData = False, axes = None, testing = False, sumMode = False, customName = None, info = False):
        """
        Method used to get or save the given summed varName data, over the given index range, using a given stepsize, along the given axes.
        
        All data is removed from dataIO after getting the relevant data.

        Parameters
        ----------
        varNames : list of str
            List of variables to save or return, must correspond to data within the BOUT-data,
            which is defined in the dumpProcessing denormalise method local dictionary,
            or a dependent variable within this method.
        tindexRange : list of int, optional
            List or array of two ints, one being minimum t index, the other maximum t index. False means use the whole index. The default is False.
        tsliceSize : int, optional
            Integer determining how large t-slices to load at a time. The default is 10.
        returnData : bool, optional
            Whether to return the data as a list, in the order given in varNames. The default is False.
        saveData : bool, optional
            Whether to save the data to the processed folder. The default is False.
        axes : tuple
            Tuple containing which axes to take the mean over. Default is (2), meaning simply the y-axis.

        Returns
        -------
        sumData : list of data, optional
            List of averaged or summed data, for the given varNames (with _sum or _averaged added to all names). Is only returned if returnData == True.

        """
        # Axes check
        if axes == None:
            raise Exception(f"'axes' is current set to {axes}. Please define the axes to mean over.")
        
        # Check for sanity
        if returnData == False and saveData == False:
            print("Both returnData and saveData are false. Exiting.")
            return

        for var in varNames:
            varcheck = [v for v in varNames if v == var]
            if len(varcheck) > 1:
                print(f"The variable '{var}' is present multiple times in the input list. Only one instance is allowed. Exiting.")
                return
        
        # Create sum/averaged names
        if sumMode == False:
            sumNames = [f"{var}_averaged_{axes}" for var in varNames]
        elif sumMode == True:
            sumNames = [f"{var}_sum_{axes}" for var in varNames]

        # Adding the custom name if need be
        if customName != None:
            sumNames = [name + f"_{customName}" for name in sumNames]
        
        ## Following two steps are nescessary because we call dataFetch and store within dataIO
        # Then, clear all data in self.IO.dataDict for summed variables, in case they are there
        for sumvar in sumNames:
            self.IO.clear_key(sumvar)
        
        # Handle the index range
        # First time-indices
        if tindexRange == None:
            tsliceMinIndex = 0
            tsliceMaxIndex = -1
        else:
            tsliceMinIndex = tindexRange[0]
            tsliceMaxIndex = tindexRange[1]

        # Get the length of the time-series as the maximum index, in case -1 is set
        if tsliceMinIndex == -1:
            tsliceMinIndex = len(collect("t_array",info=False,path=f"{paths.dataPath}/{self.dataName}"))-1
        if tsliceMaxIndex == -1:
            tsliceMaxIndex = len(collect("t_array",info=False,path=f"{paths.dataPath}/{self.dataName}"))-1
            
        # Then spatial-indices
        if xindexRange == None:
            xminIndex = 0
            xmaxIndex = -1
        else:
            xminIndex = xindexRange[0]
            xmaxIndex = xindexRange[1]
        
        # Then get the slices for partial processing
        tslices = self.generate_slices(indexRange = [tsliceMinIndex, tsliceMaxIndex], stepSize = tsliceSize)

        # In case there's only one axis to work with, convert it to a list format (used for the 1-length edge-case)
        #if type(axes) == int:
        #    axes = [axes]
        
        # Given we are working with summed data, we can assume that the data is relatively small
        # Therefore, we just hold all the data at the same time, until we need to return it or save it (held in dataIO object)
        # Then we clear the data from the IO object at the end
        if testing == True:
            print(f"Slices generated: {tslices}")
        if info == True:
            if sumMode == False:
                print(f"Taking the mean along axes: {axes}")
            else:
                print(f"Taking the sum along axes: {axes}")
        
        # For each variable load its sum in slices, and *add* these slices to the dataIO object
        for indices in tslices:
            # First get the indices from the slice
            minIndex = indices[0]
            maxIndex = indices[1]
            # Data counters must be initialised each time the indices changes
            # (doing this means data is saved if needed for more than one calculation)
            self.readyCounter(varNames = varNames)
            
            if testing == True or info == True:
                print(f"\n==== Running Interval: {minIndex}-{maxIndex} =====")
                if testing == True:
                    print(f"Initilized counter: {self.counters}")
            # Get the summed data in slices
            for i in range(len(varNames)):
                var = varNames[i]
                sumvar = sumNames[i]

                if testing == True:
                    print(f"\nGetting data for varname '{var}'")
                    print(f"Current counters: {self.counters}")
                    print(f"Taking the mean along axes {axes}")

                if info == True:
                    print(f"Getting data for variable: {var}")
                
                # Then get the summed data using these indices
                data = self.averagedVariable(var = var, minIndex = minIndex, maxIndex = maxIndex, xmaxIndex = xmaxIndex, xminIndex = xminIndex, axes = axes, testing = testing, independent = False, sumMode = sumMode)
                
                # Add the summed data it to the dataDict in IO
                if testing == True:
                    if sumvar in self.IO.dataDict.keys():
                        print(f"Shape of variable to append to: {np.shape(self.IO.dataDict[sumvar])}")
                    print(f"Shape of data to append/add, without expansion: {np.shape(data)}")
                
                self.addToDict(varNames = [sumvar], varData = [data])

                if testing == True:
                    print(f"Shape of data after appending/adding: {np.shape(self.IO.dataDict[sumvar])}")

        # Now, in case we averaged over the t/0-axis, we need to do a fix on all the data, as the code above treats the data in small chunks (intervals)
        # This means that while we may get the mean or sum along t, we only get it for a small interval. We need to combine these intervals
        # For the mean
        if sumMode == False:
            if type(axes) == int:
                tcheck = (axes == 0)
            else:
                tcheck = (0 in axes)
            
            if tcheck:
                if testing == True:
                    print(f"The axis 0 has been called in averagedVariables. Fixing data.")
                
                for sumvar in sumNames:
                    # Load data, fix it, then overwrite data with the fixed version
                    dataToFix = self.IO.returnKeys(varNames = [sumvar])[0]
                    fixedData = self.tMeanFix(data = dataToFix, tindexRange = [tsliceMinIndex, tsliceMaxIndex], tsliceSize = tsliceSize)
                    self.IO.overwriteKeys(varNames = [sumvar], varData = [fixedData])

        # For the sum
        elif sumMode == True:
            if type(axes) == int:
                tcheck = (axes == 0)
            else:
                tcheck = (0 in axes)

            if tcheck:
                for sumvar in sumNames:
                    # Load data, fix it, then overwrite data with the fixed version
                    dataToFix = self.IO.returnKeys(varNames = [sumvar])[0]
                    fixedData = np.sum(dataToFix,axis=0)
                    self.IO.overwriteKeys(varNames = [sumvar], varData = [fixedData])
        
        # Then finally save and/or get the meaned or summed data to return it here
        if saveData == True:
            self.IO.save(sumNames)
        if returnData == True:
            # Get the data
            sumData = self.IO.returnKeys(sumNames)
            
            # Clear the meaned/summed data (and everything else) from IO
            self.IO.clear_data()
            
            return sumData
        
        # Clear the data from IO in any case
        self.IO.clear_data()
        
        return

    def summedVariables(self, varNames, tindexRange = None, xindexRange = None, tsliceSize = 10, returnData = False, saveData = False, axes = None, testing = False, customName = None, info = False):
        """
        Interface for the averagedVariables method, given summed as a name instead.

        """
        if returnData == True:
            sumData = self.averagedVariables(varNames = varNames, tindexRange = tindexRange, xindexRange = xindexRange, tsliceSize = tsliceSize, returnData = returnData, saveData = saveData, axes = axes, testing = testing, sumMode = True, customName = customName, info = info)
            return sumData
            
        else:
            self.averagedVariables(varNames = varNames, tindexRange = tindexRange, xindexRange = xindexRange, tsliceSize = tsliceSize, returnData = returnData, saveData = saveData, axes = axes, testing = testing, sumMode = True, customName = customName, info = info)

    def tMeanFix(self, data, tindexRange, tsliceSize):
        """
        This method is used to turn a slices t-mean into a single t-mean (mean along t-axis), by calculating the weighted average of the averages.

        Specifically, the code here first recalculates the sums for each t-interval, then it sums all the sums together,
        resulting in the sum along the t-axis (a linear operation). As we now have the sum along the t-axis, we can now simply divide by the total amount of time,
        and we get the average along t.
        """

        # Call averagedVariables, with the first axis being the meaned axis
        # This results in data which has been devided by a slice size for each interval. Simply multiply each value by the slice size, to get back the sum of each interval
        #    - This can be done by calling generate_slices or what it is called
        # We now have the actual values summed, for each slice
        # What we now wish to do is once more use the slice size and total amount of t-values to get the weighted mean
        #    - For each value, this is the slice size divided by the total amount of t-values
        # We now have the weighted average, which is the exact same as the usual average. Well done, now return that thing

        ## We first recover the sum for each interval along the first axis
        # First get the slices for the given interval
        slices = self.generate_slices(indexRange = tindexRange, stepSize = tsliceSize)

        # Then get the length of each slice. These should be the same, except for the last one. Remember the first and last value ARE included (endpoint inclusion)
        sliceSizes = [s[1] - s[0] + 1 for s in slices]

        # We now know how many values each meaned value was calculated from, and we now recover the original sums by multiplying with *N* (sliceSizes) for each interval
        data = data
        nDim = data.ndim
        for i in range(len(sliceSizes)):
            data[i,...] = sliceSizes[i]*data[i,...]

        ## We now have all the individual sums along the first (time) axis. We now calculate the weighted average along the first axis, using the slicesizes again
        # First get the sum of values for the weighted average
        nVals = tindexRange[1] - tindexRange[0] + 1

        # Now, since the summing is a linear operation, we simply sum along the first/t-axis
        data = np.sum(data,axis=0)

        # And then we divide the data by nVals (amount of time-intervals) to get the mean data
        data = data/nVals

        return data

    def zaveragedVariables(self, varNames, tindexRange = None, xindexRange = None, tsliceSize = 10, returnData = False, saveData = False, info = False, customName = None):
        """
        Method used to get or save the given summed varName data, over the given index range, using a given stepsize, along the second (z) axis.

        Parameters
        ----------
        varNames : list of str
            List of variables to save or return, must correspond to data within the BOUT-data,
            which is defined in the dumpProcessing denormalise method local dictionary,
            or a dependent variable within this method.
        tindexRange : list of int, optional
            List or array of two ints, one being minimum t index, the other maximum t index. False means use the whole index. The default is False.
        tsliceSize : int, optional
            Integer determining how large t-slices to load at a time. The default is 10.
        returnData : bool, optional
            Whether to return the data as a list, in the order given in varNames. The default is False.
        saveData : bool, optional
            Whether to save the data to the processed folder. The default is False.

        Returns
        -------
        averagedData : list of data, optional
            List of averaged data, for the given varNames (with _averaged added to all names). Is only returned if returnData == True.
        
        """
        
        # Check for sanity
        if returnData == False and saveData == False:
            print("Both returnData and saveData are false. Exiting.")
            return
        
        # Create sum names
        averagedNames = [f"{var}_averaged_z" for var in varNames]

        # Custom name addition
        if customName != None:
            averagedNames = [f"{name}_{customName}" for name in averagedNames]
        
        # Get all variables summed along z [t,x,z,---] format
        averaged_z = self.averagedVariables(varNames, tindexRange = tindexRange, xindexRange = xindexRange, tsliceSize = tsliceSize, returnData = True, saveData = False, axes = (2), info = info)
        
        # Then save or return the data
        if saveData == True:
            for i in range(len(varNames)):
                # Get the data
                averagedName = averagedNames[i]
                data = averaged_z[i]

                # Check if the directory exists, and create it if it doesn't
                if os.path.isdir(f"{paths.processedPath}/{self.dataName}") == False:
                    os.mkdir(f"{paths.processedPath}/{self.dataName}")

                    if info == True:
                        print(f"Processed directory does not exist.\nCreating directory for data.")
                
                # Save the data
                np.save(f"{paths.processedPath}/{self.dataName}/{averagedName}.npy",data)
        if returnData == True:
            return averaged_z
        return

    def taveragedVariables(self, varNames, tindexRange = None, xindexRange = None, tsliceSize = 10, returnData = False, saveData = False, info = False, customName = None):
        """
        Method used to get or save the given summed varName data, over the given index range, using a given stepsize, along the zeroth (t) axis.

        Parameters
        ----------
        varNames : list of str
            List of variables to save or return, must correspond to data within the BOUT-data,
            which is defined in the dumpProcessing denormalise method local dictionary,
            or a dependent variable within this method.
        tindexRange : list of int, optional
            List or array of two ints, one being minimum t index, the other maximum t index. False means use the whole index. The default is False.
        tsliceSize : int, optional
            Integer determining how large t-slices to load at a time. The default is 10.
        returnData : bool, optional
            Whether to return the data as a list, in the order given in varNames. The default is False.
        saveData : bool, optional
            Whether to save the data to the processed folder. The default is False.

        Returns
        -------
        averagedData : list of data, optional
            List of averaged data, for the given varNames (with _averaged added to all names). Is only returned if returnData == True.
        
        """
        
        # Check for sanity
        if returnData == False and saveData == False:
            print("Both returnData and saveData are false. Exiting.")
            return
        
        # Create sum names
        averagedNames = [f"{var}_averaged_t" for var in varNames]

        # Custom name addition
        if customName != None:
            averagedNames = [f"{name}_{customName}" for name in averagedNames]
        
        # Get all variables summed along t [t,x,z,---] format
        averaged_t = self.averagedVariables(varNames, tindexRange = tindexRange, xindexRange = xindexRange, tsliceSize = tsliceSize, returnData = True, saveData = False, axes = (0), info = info)
        
        # Then save or return the data
        if saveData == True:
            for i in range(len(varNames)):
                # Get the data
                averagedName = averagedNames[i]
                data = averaged_t[i]

                # Check if the directory exists, and create it if it doesn't
                if os.path.isdir(f"{paths.processedPath}/{self.dataName}") == False:
                    os.mkdir(f"{paths.processedPath}/{self.dataName}")

                    if info == True:
                        print(f"Processed directory does not exist.\nCreating directory for data.")
                
                # Save the data
                np.save(f"{paths.processedPath}/{self.dataName}/{averagedName}.npy",data)
        if returnData == True:
            return averaged_t
        return

    def tzaveragedVariables(self, varNames, tindexRange = None, xindexRange = None, tsliceSize = 10, returnData = False, saveData = False, info = False, customName = None):
        """
        Method used to get or save the given summed varName data, over the given index range, using a given stepsize, along the zeroth (t) and second (z) axes. Remember that the 1-length y-axis gets squezed away during dumpProcessing.

        Parameters
        ----------
        varNames : list of str
            List of variables to save or return, must correspond to data within the BOUT-data,
            which is defined in the dumpProcessing denormalise method local dictionary,
            or a dependent variable within this method.
        tindexRange : list of int, optional
            List or array of two ints, one being minimum t index, the other maximum t index. False means use the whole index. The default is False.
        tsliceSize : int, optional
            Integer determining how large t-slices to load at a time. The default is 10.
        returnData : bool, optional
            Whether to return the data as a list, in the order given in varNames. The default is False.
        saveData : bool, optional
            Whether to save the data to the processed folder. The default is False.

        Returns
        -------
        averagedData : list of data, optional
            List of averaged data, for the given varNames (with _averaged added to all names). Is only returned if returnData == True.
        
        """
        
        # Check for sanity
        if returnData == False and saveData == False:
            print("Both returnData and saveData are false. Exiting.")
            return
        
        # Create sum names
        averagedNames = [f"{var}_averaged_tz" for var in varNames]

        # Custom name addition
        if customName != None:
            averagedNames = [f"{name}_{customName}" for name in averagedNames]
        
        # Get all variables summed along tz [t,x,z,v] format
        averaged_tz = self.averagedVariables(varNames, tindexRange = tindexRange, xindexRange = xindexRange, tsliceSize = tsliceSize, returnData = True, saveData = False, axes = (0,2), info = info)
        
        # Then save or return the data
        if saveData == True:
            for i in range(len(varNames)):
                # Get the data
                averagedName = averagedNames[i]
                data = averaged_tz[i]
                
                # Check if the directory exists, and create it if it doesn't
                if os.path.isdir(f"{paths.processedPath}/{self.dataName}") == False:
                    os.mkdir(f"{paths.processedPath}/{self.dataName}")

                    if info == True:
                        print(f"Processed directory does not exist.\nCreating directory for data.")

                # Then save the data
                np.save(f"{paths.processedPath}/{self.dataName}/{averagedName}.npy",data)
        if returnData == True:
            return averaged_tz
        return
    
    def zsummedVariables(self, varNames, tindexRange = None, xindexRange = None, tsliceSize = 10, returnData = False, saveData = False, info = False, customName = None):
        """
        Method used to get or save the given summed varName data, over the given index range, using a given stepsize, along the second (z) axis.

        Parameters
        ----------
        varNames : list of str
            List of variables to save or return, must correspond to data within the BOUT-data,
            which is defined in the dumpProcessing denormalise method local dictionary,
            or a dependent variable within this method.
        tindexRange : list of int, optional
            List or array of two ints, one being minimum t index, the other maximum t index. False means use the whole index. The default is False.
        tsliceSize : int, optional
            Integer determining how large t-slices to load at a time. The default is 10.
        returnData : bool, optional
            Whether to return the data as a list, in the order given in varNames. The default is False.
        saveData : bool, optional
            Whether to save the data to the processed folder. The default is False.

        Returns
        -------
        sumData : list of data, optional
            List of summed data, for the given varNames (with _sum added to all names). Is only returned if returnData == True.
        
        """
        
        # Check for sanity
        if returnData == False and saveData == False:
            print("Both returnData and saveData are false. Exiting.")
            return
        
        # Create sum names
        sumNames = [f"{var}_sum_z" for var in varNames]

        # Custom name addition
        if customName != None:
            sumNames = [f"{name}_{customName}" for name in sumNames]
        
        # Get all variables summed along z [t,x,z,---] format
        summed_z = self.summedVariables(varNames, tindexRange = tindexRange, xindexRange = xindexRange, tsliceSize = tsliceSize, returnData = True, saveData = False, axes = (2), info = info)
        
        # Then save or return the data
        if saveData == True:
            for i in range(len(varNames)):
                # Get the data
                sumName = sumNames[i]
                data = summed_z[i]

                # Check if the directory exists, and create it if it doesn't
                if os.path.isdir(f"{paths.processedPath}/{self.dataName}") == False:
                    os.mkdir(f"{paths.processedPath}/{self.dataName}")

                    if info == True:
                        print(f"Processed directory does not exist.\nCreating directory for data.")

                # Save the data
                np.save(f"{paths.processedPath}/{self.dataName}/{sumName}.npy",data)
        if returnData == True:
            return summed_z
        return

    def tsummedVariables(self, varNames, tindexRange = None, xindexRange = None, tsliceSize = 10, returnData = False, saveData = False, info = False, customName = None):
        """
        Method used to get or save the given summed varName data, over the given index range, using a given stepsize, along the zeroth (t) axis.

        Parameters
        ----------
        varNames : list of str
            List of variables to save or return, must correspond to data within the BOUT-data,
            which is defined in the dumpProcessing denormalise method local dictionary,
            or a dependent variable within this method.
        tindexRange : list of int, optional
            List or array of two ints, one being minimum t index, the other maximum t index. False means use the whole index. The default is False.
        tsliceSize : int, optional
            Integer determining how large t-slices to load at a time. The default is 10.
        returnData : bool, optional
            Whether to return the data as a list, in the order given in varNames. The default is False.
        saveData : bool, optional
            Whether to save the data to the processed folder. The default is False.

        Returns
        -------
        sumData : list of data, optional
            List of summed data, for the given varNames (with _sum added to all names). Is only returned if returnData == True.
        
        """
        
        # Check for sanity
        if returnData == False and saveData == False:
            print("Both returnData and saveData are false. Exiting.")
            return
        
        # Create sum names
        sumNames = [f"{var}_sum_t" for var in varNames]

        # Custom name addition
        if customName != None:
            sumNames = [f"{name}_{customName}" for name in sumNames]
        
        # Get all variables summed along t [t,x,z,---] format
        summed_z = self.summedVariables(varNames, tindexRange = tindexRange, xindexRange = xindexRange, tsliceSize = tsliceSize, returnData = True, saveData = False, axes = (0), info = info)
        
        # Then save or return the data
        if saveData == True:
            for i in range(len(varNames)):
                # Get the data
                sumName = sumNames[i]
                data = summed_t[i]

                # Check if the directory exists, and create it if it doesn't
                if os.path.isdir(f"{paths.processedPath}/{self.dataName}") == False:
                    os.mkdir(f"{paths.processedPath}/{self.dataName}")

                    if info == True:
                        print(f"Processed directory does not exist.\nCreating directory for data.")

                # Save the data
                np.save(f"{paths.processedPath}/{self.dataName}/{sumName}.npy",data)
        if returnData == True:
            return summed_t
        return

    def tzsummedVariables(self, varNames, tindexRange = None, xindexRange = None, tsliceSize = 10, returnData = False, saveData = False, info = False, customName = None):
        """
        Method used to get or save the given summed varName data, over the given index range, using a given stepsize, along the zeroth (t) and second (z) axes. Remember that the 1-length y-axis gets sqeuuzed away during dumpProcessing.

        Parameters
        ----------
        varNames : list of str
            List of variables to save or return, must correspond to data within the BOUT-data,
            which is defined in the dumpProcessing denormalise method local dictionary,
            or a dependent variable within this method.
        tindexRange : list of int, optional
            List or array of two ints, one being minimum t index, the other maximum t index. False means use the whole index. The default is False.
        tsliceSize : int, optional
            Integer determining how large t-slices to load at a time. The default is 10.
        returnData : bool, optional
            Whether to return the data as a list, in the order given in varNames. The default is False.
        saveData : bool, optional
            Whether to save the data to the processed folder. The default is False.

        Returns
        -------
        sumData : list of data, optional
            List of summed data, for the given varNames (with _sum added to all names). Is only returned if returnData == True.
        
        """
        
        # Check for sanity
        if returnData == False and saveData == False:
            print("Both returnData and saveData are false. Exiting.")
            return
        
        # Create sum names
        sumNames = [f"{var}_sum_tz" for var in varNames]

        # Custom name addition
        if customName != None:
            sumNames = [f"{name}_{customName}" for name in sumNames]
        
        # Get all variables summed along z [t,x,z,v] format
        summed_tz = self.summedVariables(varNames, tindexRange = tindexRange, xindexRange = xindexRange, tsliceSize = tsliceSize, returnData = True, saveData = False, axes = (0,2), info = info)
        
        # Then save or return the data
        if saveData == True:
            for i in range(len(varNames)):
                # Get the data
                sumName = sumNames[i]
                data = summed_tz[i]

                # Check if the directory exists, and create it if it doesn't
                if os.path.isdir(f"{paths.processedPath}/{self.dataName}") == False:
                    os.mkdir(f"{paths.processedPath}/{self.dataName}")

                    if info == True:
                        print(f"Processed directory does not exist.\nCreating directory for data.")

                # Save the data
                np.save(f"{paths.processedPath}/{self.dataName}/{sumName}.npy",data)
        if returnData == True:
            return summed_tz
        return
