### Libraries imports
import sys
import numpy as np

### Locals imports
# Import paths file
pathPath = "/home/s194112/"
sys.path.append(pathPath)
import paths

# Import job submission script
sys.path.append(paths.scriptPath)

# Import scripts
from inputHandler import inputHandler
from submissions import simsubmissions

def dictionaryListInterpreter(dic, roundParameters):
    """
    Function used to convert an input dictionary such as {"INNERTEMP" : [100, 400, 100]} into the expanded format {"INNERTEMP" : [100, 200, 300, 400]}, by interpreting 3-value lists as initial value, final value, and step size. Note the end-point inclusivity, meaning the first and final point will always be included.

    """
    
    ## First construct the lists to sweep for each parameter
    parameterSweeps = {}
    for key in dic.keys():
        try:
            pLength = len(dic[key])
        except:
            if type(dic[key]) == int or type(dic[key]) == float:
                pLength = 1
        
        # For single-value parameter
        if pLength == 1:
            # Rounding
            if (roundParameters != None) and (key in roundParameters.keys()):
                parameterSweeps[key] = round(dic[key], roundParameters[key])
            # Not rounding
            else:
                parameterSweeps[key] = dic[key]

        # For multi-value parameter
        elif pLength == 3:
            # First extract the lower, upper, and step values
            lower = dic[key][0]
            upper = dic[key][1]
            step = dic[key][2]

            # Then construct the steps, ensuring endpoint-inclusiveness
            valList = np.arange(lower,upper,step)

            # Potential rounding
            if (roundParameters != None) and (key in roundParameters.keys()):
                for i in range(len(valList)):
                    valList[i] = round(valList[i], roundParameters[key])
            
            if valList[-1] != upper:
                valList = np.append(valList, upper)

            # Finally add it to the sweeps list
            parameterSweeps[key] = valList
            

        else:
            raise Exception(f"The key {key} has {pLength} elements, but either 1 or 3 were expected.")

    return parameterSweeps

def sweepCombinations(parameters, keys = None, index = None, branchListKeys = None, branchListVals = None):
    """
    Function used to generate all unique combinations of the given input parameters, in dictionary form, which are returned as a list of dictionaries.

    """
    
    # For each value in parameters[key], call sweepName with 'parameters', without 'key'
    # Unless 'key' is the only key present in 'parameters'
    
    # Pass a string down in recursive case, and add it to 'sweepNames' if base-case reached
    
    # Input check
    if len(parameters.keys()) < 1:
        raise Exception("There are less then 1 key in the input dictionary. Please check your inputs.")
    
    # Initialize keys
    if keys == None:
        keys = [k for k in parameters.keys()]

    # Initialize index
    if index == None:
        index = 0

    # Initialize lists
    if branchListKeys == None:
        branchListKeys = []
    if branchListVals == None:
        branchListVals = []
    
    # For loop
    nameList = []
    
    key = keys[index]
    for val in parameters[key]:

        # Base-case
        if index == len(keys)-1:
            # Construct dictionary with full branch path of keys and vals
            tempDict = {}
            for i in range(len(branchListKeys)):
                tempKey = branchListKeys[i]
                tempVal = branchListVals[i]
                tempDict[tempKey] = tempVal
            tempDict[key] = val

            # Then append dictionary to namelist
            nameList.append(tempDict)
        
        # Recursive case
        else:
            # First create two temporary lists to pass to recursive case
            # Keys list
            branchListKeysCopy = branchListKeys.copy()
            branchListKeysCopy.append(key)

            # Values list
            branchListValsCopy = branchListVals.copy()
            branchListValsCopy.append(val)

            # Then call function
            tempList = sweepCombinations(parameters = parameters, keys = keys, index = index+1, branchListKeys = branchListKeysCopy, branchListVals = branchListValsCopy)

            # Add returned namelist to namelist
            nameList.append(tempList)

    return nameList

def sweepName(parameters):
    """
    Function used to generate the base-name (without suffix) for use in the parameterSweep function.

    Input is a dictionary of parameters to generate name for.

    Note: values using '.' in their names will have the dot replaced by a comma (',') due to issues with file extensions.
    
    """
    
    # First define abridged key list
    keysAbridged = {"INNERTEMP" : "Tin", "NOUT" : "n", "INNERDENSITY" : "nin"}

    # Then required variables for loop
    sweepName = ""
    keys = [key for key in parameters.keys()]
    dicLen = len(keys)

    # Then construct the sweep name
    for i in range(dicLen):
        key = keys[i]
        
        if key in keysAbridged.keys():
            sweepName += f"{keysAbridged[key]}{parameters[key]}"
        else:
            sweepName += f"{key}{parameters[key]}"

        if i != dicLen-1:
            sweepName += "_"

    # Fix for file-extensions
    sweepName = sweepName.replace(".",",")

    return sweepName

def nameGenerator(dic, customSuffix = None, roundParameters = None):
    """
    Function used to mimic the names resulting from using the parameterSweep function.

    Rather useful in further processing or visualisation of multiple datafolder.

    The inputs are the same as in parameterSweep, but with only the dictionary and name.

    """
    # Interpret dictionary lists
    dic = dictionaryListInterpreter(dic, roundParameters)

    # Ensure no single-value (non-iterables) are present
    for key in dic.keys():
        # Try iterating
        try:
            for val in dic[key]:
                pass
        # If you can't, it is probably a singular value. Make it iterable.
        except:
            dic[key] = [dic[key]]
    
    # Then pass to sweepcombinations to generate list of combinations
    listCombies = sweepCombinations(dic)

    # Also flatten the list
    listCombies = np.array(listCombies).flatten()

    # Then apply possible rounding
    if roundParameters != None:
        for e in listCombies:
            for key in e:
                # If the key in each element is present in the roundParameters keys, round the value
                if key in roundParameters.keys():
                    e[key] = round(e[key], roundParameters[key])
                # Otherwise, do not round it
                else:
                    pass
    
    # Finally generate a list of sweep names
    listParamNames = listNames = [f"{sweepName(listDic)}" for listDic in listCombies]
    if customSuffix == None:
        if len(listCombies[0].keys()) == 1:
            listDataNames = [f"data_{sweepName(listDic)}_Singular" for listDic in listCombies]
        else:
            listDataNames = [f"data_{sweepName(listDic)}_Full" for listDic in listCombies]
    else:
        listDataNames = [f"data_{sweepName(listDic)}_{customSuffix}" for listDic in listCombies]
    
    # And for convenience, also generate the exact paths to the respective *processed* folders
    listProcessed = [f"{paths.processedPath}/{name}" for name in listDataNames]

    return listCombies, listParamNames, listDataNames, listProcessed


def parameterSweep(parameters, fullSweep = True, keepInp = True, customSuffix = None, includeSim = True, includeProcessing = True, explicitLists = False, roundParameters = None):
    """
    Function usedd to perform a sweep of the given parameters, with the given stepsizes.
    Can be used to perform a full sweep of the given parameter space (meaning all combinations),
    or simple single-line sweeps for each parameter, where the rest are set to default values.

    The list of parameters are current:
     - INNERTEMP [UNIT: eV] [DEFAULT: 100 eV]
        - Initial (and forced) temperature of electrons and ions at the inner part of the simulation. Used to calculate the initial pressures in the system.
     - INNERDENSITY [UNIT: e19 m^-3] [DEFAULT: 3.5]
        - Initial (and forced) density of the plasma at the inner part of the simulation.
     - NOUT [UNIT: 100 time units] [DEFAULT: 500]
        - Amount of time-outputs, with the total sim time being nout*timestep in the inp file. 'timestep' is 100
     - SEED [Unit: None] [DEFAULT: Current Unix Timestamp with ms precision]
        - Seed used to initialize simulations with (is a phase in the 'mixmode' function, see: https://bout-dev.readthedocs.io/en/latest/user_docs/variable_init.html)

    A more comprehensive description can be found in the inputHandler file (or rather, a more up-to-date list in case I've updated things again).

    The 'parameters' dictionary must contains either:
     - a singular [value] for the parameter
     - 3 values for the [lower bound, upper bound, stepsize]

    Note that this function is enpoint-inclusive, meaning the upper bound will *always* be included in the sweep.

    Inputs:
     - parameters [dictionary]
         - dictionary of parameter names, with each corresponding key containing either a singular value, or an upper bound, lower bound, and stepsize
     - fullSweep [bool, Default: True]
         - Whether the sweep should be done using all available parameter combinations, or simply for one parameter at a time, with the rest being default values
     - keepInp [boo, Default: True]
         - Whether the input files created during the sweep should be kept or not, after job submission
     - customSuffix [str, Default: None]
         - A possible custom suffix, instead of the _Singular suffix for single-parameter sweeps, or _Full for multi-parameter sweeps.
     - includeSim [bool, Default: True]
         - Whether to include the simulation in the sweep, or only do the postprocessing of sim data. Can be useful for re-processing with different parameters.
     - includeProcessing [bool, Default: True]
         - Whether to include the processing in the sweep. Can be useful for re-creating default figures in a lot of folder. Note setting includeSim = False is a good idea when doing this.
     - explicitLists [bool, Default: False]
         - Whether to turn off the parameter dictionary list comprehension. Setting this to true will disable the requirements of 1 or 3 elements in parameter lists, and they will be taken as-is, instead of being taken as lower, upper, step. Useful for direct inputs. Note this does NOT affect the fullsweep options or anything else. This is explicitly to allow for custom parameter value inputs.
     - roundParameters [dictionary, Default: None]
         - A list corresponding to the same values input in the dictionary. If nothing is input, no rounding is applied. Rounds to the given amount of digits for the given parameters in the dictionary. Note this is only implemented for the fullSwepp option.

    """
    ## First construct the lists to sweep for each parameter
    if explicitLists == False:
        ### WORKS, BUT CONSOLIDATED INTO ANOTHER FUNCTION ###
        # parameterSweeps = {}
        # for key in parameters.keys():
        #     try:
        #         pLength = len(parameters[key])
        #     except:
        #         if type(parameters[key]) == int or type(parameters[key]) == float: # Could probably use some better way of checking for single value type
        #             pLength = 1
            
        #     # For single-value parameter
        #     if pLength == 1:
        #         parameterSweeps[key] = parameters[key]
    
        #     # For multi-value parameter
        #     elif pLength == 3:
        #         # First extract the lower, upper, and step values
        #         lower = parameters[key][0]
        #         upper = parameters[key][1]
        #         step = parameters[key][2]
    
        #         # Then construct the steps, ensuring endpoint-inclusiveness
        #         valList = np.arange(lower,upper,step)
        #         if valList[-1] != upper:
        #             valList = np.append(valList, upper)
    
        #         # Finally add it to the sweeps list
        #         parameterSweeps[key] = valList
                
    
        #     else:
        #         raise Exception(f"The key {key} has {pLength} elements, but either 1 or 3 were expected.")
        ### WORKS, BUT CONSOLIDATED INTO ANOTHER FUNCTION ###

        parameterSweeps = dictionaryListInterpreter(dic = parameters, roundParameters = roundParameters)
    else:
        parameterSweeps = parameters

    ## Now do the actual sweeps
    # Start by constructing an inputHandler and simsubmission object
    inputObj = inputHandler()

    # Then do a quick fix ensuring all the values for all keys are iterables
    for key in parameterSweeps.keys():
        # Can data for key be iterated upon?
        try:
            for val in parameterSweeps[key]:
                # Do nothing
                pass
        # If not, assume it to be a singular value and make it iterable
        except:
            parameterSweeps[key] = [parameterSweeps[key]]
    
    # For the singular-parameter sweep
    if fullSweep == False:
        # Sweep over keys
        for key in parameterSweeps.keys():
            
            # Sweep over steps for key
            for val in parameterSweeps[key]:
                # First construct inputname and corresponding file
                if customSuffix == None:
                    inputName = sweepName({key : val}) + "_Singular"
                else:
                    inputName = sweepName({key : val}) + f"_{customSuffix}"
                inputObj.create_inp(name = inputName, parameters = {key : val}, overwrite = True)

                # Then submit a simulation job using this input file
                simObj = simsubmissions(inputFile = inputName, overwriteData = True)
                simObj.do_full_process(includeSim = includeSim, includeProcessing = includeProcessing)

                # Then get rid of input file if need be
                if keepInp != True:
                    inputObj.delete_inp(name = inputName)
            
                

    # For the full-parameter sweep
    else:
        # Call combinationsdict
        combinationsList = sweepCombinations(parameterSweeps)

        # Then flatten the list
        combinationsList = np.array(combinationsList).flatten()

        # Now apply the possible rounding
        if roundParameters != None:
            for key in roundParameters.keys():
                # Get the value for the rounding
                roundNumber = roundParameters[key]
                
                for combi in combinationsList:
                    # Round the applicable values
                    combi[key] = round(combi[key], roundNumber)

        # We now have a full 1D list, containing all possible parameter combinations
        # For each combination (read: element of list), get the sweep name and run a job submission
        for combination in combinationsList:
            # First construct inputname and corresponding file
            if customSuffix == None:
                inputName = sweepName(combination) + "_Full"
            else:
                inputName = sweepName(combination) + f"_{customSuffix}"
            inputObj.create_inp(name = inputName, parameters = combination, overwrite = True)

            # Then submit a simulation job using this input file
            simObj = simsubmissions(inputFile = inputName, overwriteData = True)
            simObj.do_full_process(includeSim = includeSim, includeProcessing = includeProcessing)

            # Then get rid of input file if need be
            if keepInp != True:
                inputObj.delete_inp(name = inputName)

            #print(inputName)
    
    print("Sweep submitted!")
