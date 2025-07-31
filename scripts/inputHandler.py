### Libraries imports
import sys
from pathlib import Path
import warnings
import os
from datetime import datetime

### Locals imports
# Import paths file
pathPath = "/home/s194112/"
sys.path.append(pathPath)
import paths

# Import job submission script
sys.path.append(paths.scriptPath)

class inputHandler:
    """
    Class used to create or delete an input file in the input folder, with the given name for the given parameters.
    parameters are given as a dictionary, with all non-specified parameters using default values, given below.
    For each string variable, give a corresponding integer or float value.
    
    Note that variable names are NOT case-sensitive.
    
    Template parameters:
     - INNERTEMP [UNIT: eV] [DEFAULT: 500 eV]
        - Initial (and forced) temperature of electrons and ions at the inner part of the simulation. Used to calculate the initial pressures in the system. Note, this value will overwrite whichever of INNERTEMP_ION and INNERTEMP_ELECTRON which have *not* been defined.
     - INNERTEMP_IONS [UNIT: eV] [DEFAULT: INNTERTEMP]
         - Used to specify the temperature of the ions. Used for calculating forced pressure profile. If not specified, uses INNERTEMP value.
     - INNERTEMP_ELECTRONS [UNIT: eV] [DEFAULT: INNERTEMP]
         - Used to specify the temperature of the electrons. Used for calculating forced pressure profile. If not specified, uses INNERTEMP value.
     - INNERDENSITY [UNIT: e19 m^-3] [DEFAULT: 3.5]
        - Initial (and forced) density of the plasma at the inner part of the simulation.
     - NOUT [UNIT: 100 time units] [DEFAULT: 600]
        - Amount of time-outputs, with the total sim time being nout*timestep in the inp file. 'timestep' is 100
     - SEED [Unit: None] [DEFAULT: Current Unix Timestamp with ms precision]
        - Seed used to initialize simulations with (is a phase in the 'mixmode' function, see: https://bout-dev.readthedocs.io/en/latest/user_docs/variable_init.html)
    
    """
    
    def __init__(self):
        # Import the template file as a string
        self.templateString = Path(f"{paths.templatePath}/TEMPLATE_input.inp").read_text()

        # Define the default parameters
        self.defaultParameters = {"INNERTEMP" : 500, "INNERDENSITY" : 3.5, 'NOUT' : 500, 'SEED' : None} # Remember to put everything in upper case

    def create_inp(self, name, parameters = None, overwrite = True):
        """
        Method used to create an input file in the input folder, with the given name for the given parameters.
        parameters are given as a dictionary, with all non-specified parameters using default values, given below.
        For each string variable, give a corresponding integer or float value.
    
        Note that variable names are NOT case-sensitive.

        See inputHandler main documentation for a list of available inputs.
    
        """

        # Start by updating the default SEED value
        self.defaultParameters["SEED"] = datetime.now().timestamp()
        
        # Now remove case-sensitivity of keys
        if parameters != None:
            for key in [k for k in parameters.keys()]:
                parameters[key.upper()] = parameters.pop(key)
        
        # Check for full default parameter usage
        if parameters == None:
            warnings.warn(f"The given parameters are defined as {parameters}, and the full default parameters will therefore be used.")
            parameters = (self.defaultParameters).copy()
            fullDefault = True

            # And handle the electron and ion temperatures as special cases
            parameters["INNERTEMP_IONS"] = parameters["INNERTEMP"]
            parameters["INNERTEMP_ELECTRONS"] = parameters["INNERTEMP"]
        else:
            fullDefault = False

        # Then partial default parameter usage
        if fullDefault == False:
            # If we are not utilizing full default parameters, replace each key not present in 'parameters' with the default value
            defaultKeys = []
            for key in self.defaultParameters.keys():
                if key not in parameters:
                    parameters[key] = self.defaultParameters[key]
                    defaultKeys.append(key)
            
            # Check for TINNER_ION and TINNER_ELECTRON due to special case of T_INNER overwrite                        
            # overwrite whichever key is *not* present in parameters
            if "INNERTEMP_IONS" not in parameters.keys():
                parameters["INNERTEMP_IONS"] = parameters["INNERTEMP"]
                defaultKeys.append("INNERTEMP_IONS")
                ionsDefined = False
            else:
                ionsDefined = True
            if "INNERTEMP_ELECTRONS" not in parameters.keys():
                parameters["INNERTEMP_ELECTRONS"] = parameters["INNERTEMP"]
                defaultKeys.append("INNERTEMP_ELECTRONS")
                electronsDefined = False
            else:
                electronsDefined = True

            # In case they are both defined, remove INNERTEMP from used default keys - don't do this before the above, as they won't have anythign to pull from otherwise
            if ionsDefined and electronsDefined:
                if "INNERTEMP" in defaultKeys:
                    defaultKeys.remove("INNERTEMP")

            # Print used default keys
            if len(defaultKeys) != 0:
                print(f"The following default keys have been utilised: {defaultKeys}")
                    

        # Then go over each parameter, replacing the value in the template string with the given value
        tempStr = self.templateString
        for key in parameters.keys():
            key = key.upper() # Remove case-sensitivity
            val = parameters[key]
            tempStr = tempStr.replace(f"<{key}>", str(val))

        # Then write to file
        if os.path.exists(f"{paths.inputPath}/{name}.inp") == True and overwrite == False:
            raise Exception(f"The file at {paths.inputPath}/{name}.inp already exists, and overwrite is set to {overwrite}.")
        else:
            file = open(f"{paths.inputPath}/{name}.inp", "w")
            file.write(tempStr)
            file.close()

        return 0
            
    def delete_inp(self, name):
        # First check if the file exists
        if os.path.exists(f"{paths.inputPath}/{name}.inp") == False:
            warnings.warn(f"The file at {paths.inputPath}/{name}.inp does not exist, and can therefore not be deleted.")
            return 1
        else:
            os.remove(f"{paths.inputPath}/{name}.inp")
            return 0