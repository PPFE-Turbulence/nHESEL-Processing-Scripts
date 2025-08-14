### Libraries imports
import sys
import numpy as np
import math as m
import matplotlib.pyplot as plt
import warnings
from scipy.optimize import fsolve
from scipy.integrate import simpson

### Locals imports
# Import paths file
pathPath = os.environ["HOME"]
sys.path.append(pathPath)
import paths

# Import job submission script
sys.path.append(paths.scriptPath)

# Import scripts
from downstreamMap import downstreamMap
from fxEstimates import fxEstimates
from maxMeans import maxMeansEnergyflux
from dumpProcessing import dumpProcessing
from FOLs import FOLs
from P_loss import P_loss
from axialDistances import axialDistances

### Script
def q(s, dataName, fxType = False, q0Type = "maxMeans", plasmaElongation = 1.29729, diffusionCorrection = 1.1, extendedOutput = False, overWriteS = None, overWritefx = None, overWriteLambda_q = None, overWriteq0 = None, q0Correction = 0.5, eOri = None, SoldMethod = False, Bold = False, info = False):
    """
    Function used to calculate the downstream transform of the heat profile from the midplane, unto the diverter.

    Input is a list of points on the diverter, relative to the downstream projected LCFS position, and the name of a processed datafolder.
    
    If extendedOutput is True, then the output is given as qs, S, B_pol, n_inner, a, I_p, kappa, lambda_q, fx.

    lambda_q, fx, and S, as well as q0 used in the calculations can all be overwritten by inputting overWrite{variablename}.

    The diffusioncorrection factor is what the power output is modified by, to account for diffusion which is not otherwise accounted for (not implemented).

    """

    ### Handle user input
    ## Map DFOL and PFOL to numbers
    # If str input
    if type(fxType) == str:
        fxType = fxType.lower()
        fxDict = {"pfol" : 0, "dfol" : 1, "const" : 2}
        fxType = fxDict[fxType]

    # Then check if valid numerical input
    if fxType not in [0,1, 2]:
        raise Exception(f"fxType is currently given as {fxType}, but the only allowed inputs are 'PFOL' and 'DFOL', or their corresponding unmbers 0 and 1.")

    ## Map q0Type to numbers
    # If str input
    if type(q0Type) == str:
        q0Type = q0Type.lower()
        q0Dict = {"max" : 0, "maxmeans" : 1, "lcfs" : 2}
        q0Type = q0Dict[q0Type]

    # Then check if valid numerical input
    if q0Type not in [0,1,2]:
        raise Exception(f"q0Type is currently given as {q0Type}, but the only allowed inputs are 'max' and 'maxMeans', or their corresponding unmbers 0 and 1.")
    
    ### Get all required values for projection q(s)
    ## Calculate the q0 value
    # First get both possible q0 values
    maxEnergy, maxMeanEnergy, toleranceState = maxMeansEnergyflux(dataName)

    # Check tolarance
    if toleranceState == True:
        warnings.warn("The tolerance on calculating the means energy flux was broken. This doesn't invalidate anything, but do be aware values from the transistion zone across the LCFS have been included in the mean, reducing its value!")
    
    ## Then get the fx value
    fx_PFOL, fx_DFOL = fxEstimates(dataName)
    fx_const = fxEstimates(dataName, constMode = True)

    if fxType == 0:
        fx = fx_PFOL
    elif fxType == 1:
        fx = fx_DFOL
    elif fxType == 2:
        fx = fx_const
    if overWritefx != None:
        fx = overWritefx

    ## Get the PFOL value
    PFOL, DFOL = FOLs(dataName, eOri = eOri)
    lambda_q = PFOL
    if overWriteLambda_q != None:
        lambda_q = overWriteLambda_q
    
    ## Finally calculate the diffusion factor S
    # First get the kappa from input
    kappa = plasmaElongation

    # Then load required quantities plasma current (I_p), minor radius (a), and edge electron density (the innermost part of the density n)
    # Plasma current
    I_p = 4.60e6 # Defined from excel file [Amps]

    # Minor radius
    dumpObj = dumpProcessing(dataName)
    a = dumpObj.get_option_value("Rminor") # [m]

    # Electron edge density
    n = np.load(f"{paths.processedPath}/{dataName}/n_averaged_z.npy") # Indices [t,x,z] -> [t,x]
    n_inner = n[0,0]/1e19 # Given in units of 1e19 m^-3

    # Vacuum permeability
    mu_0 = 1.25663706127e-6 # [N*Amp]
    
    # Calculate the poloidal magnetic field
    if Bold == True:
        B_pol = (mu_0*I_p)/(2*m.pi*a*m.sqrt((1+kappa**2)/2)) # Results in a bit too high value, but good to compare to
    else:
        B_pol = 0.51 # T

    # Finally combine all of this to get an approximation of the diffusion/advection factor S
    if overWriteS != None:
        S = overWriteS
    else:
        if SoldMethod == True:
            S = 0.09*(n_inner**1.02)*(B_pol**(-1.01)) # mm? # Old method, resulted in too low value
            S = S/1e3 # mm -> m # conversion of unit
        else:
            S = iterateS(dataName = dataName)

    ### Now apply the q(s) formula, to get the downstream q(s) profile
    # qs = np.array([(q0/2) * m.exp((S/(2*lambda_q))**2 - pos/(lambda_q*fx)) * m.erfc(S/(2*lambda_q) - pos/(S*fx)) for pos in s])
    qs = np.array([(1/2) * m.exp((S/(2*lambda_q))**2 - pos/(lambda_q*fx)) * m.erfc(S/(2*lambda_q) - pos/(S*fx)) for pos in s]) # Notice, q0 missing

    ### Goodie, now we can correct the above, by applying the q0 value, unless you want to overwrite it with something instead
    ## In case you want to overwrite it
    if overWriteq0 != None:
        q0 = overWriteq0
    ## In case you don't, we first get P_loss for this data, then we multiply each point with its radius,nd then integrate, after which we figure out what q0 should be at the diverter entrace to match P_loss
    else:
        # First get the R value for each s point
        Rs = axialDistances(s)

        # Then use this to get the qs*circumference values to get [W/m]
        q_times_circumferances = qs*2*m.pi*Rs

        # Now integrate these values to get [W]
        P_diverter = simpson(q_times_circumferances, s)

        # Then get the P_loss to compare to
        Pval = P_loss(dataName, diffusionCorrection = diffusionCorrection, eOri = eOri)

        # Finally use this to figure out what scale we need to multiply with
        q0 = Pval/P_diverter

        # And now for the qCorrection - basically the q0 right now is for *everything* from the midplane - this can be corrected for
        if q0Correction != None:
            q0 = q0*q0Correction

        if info == True:
            print(f"Ploss: {Pval/1e6} [MW]")
            print(f"Pdiverter: {P_diverter/1e6} [MW]")
            print(f"q0: {q0/1e6} [MW]")

    ## And then apply it
    qs = q0*qs
    
    ### In case addition information was requested
    if extendedOutput == True:
        qs = (qs, S, B_pol, n_inner, a, I_p, kappa, lambda_q, fx)

    # And info dump
    if info == True:
        print(f"B_pol: {B_pol}")
        print(f"n_inner: {n_inner}")
        print(f"q0: {q0}")
        print(f"S: {S}")
        print(f"fx: {fx}")
        print(f"2*lambda_q: {2*lambda_q}")
        print(f"(S/(2*lambda_q)): {(S/(2*lambda_q))}")
        print(f"(S/(2*lambda_q))**2: {(S/(2*lambda_q))**2}")
        print(f"s[20]/(lambda_q*fx): {s[20]/(lambda_q*fx)}")
        print()
    
    return qs

def intFOL_works(S, dataName, eOri = None, overWriteLambda_q = None):
    # Calculating the integral power fall-off length
    s = np.linspace(-5,5,num=500)
    qs = q(s, dataName = dataName, overWriteS = S, eOri = eOri, overWriteLambda_q = overWriteLambda_q)
    lambda_int = simpson(qs*s,s)/simpson(qs,s)
    return lambda_int

def intFOL_worksnot(S, dataName, eOri = None, overWriteLambda_q = None):
    # Calculating the integral power fall-off length as should be correct from paper (it is not!)
    s = np.linspace(-5,5,num=500)
    qs = q(s, dataName = dataName, overWriteS = S, eOri = eOri, overWriteLambda_q = overWriteLambda_q)
    lambda_int = simpson(qs, s)/np.max(qs)
    return lambda_int

def intFOL_worksmaybe(S, dataName):
    # Calculating the integral power fall-off length directly from numerical FWHM
    # First get the s and qs values
    s = np.linspace(-5,5,num=500)
    qs = q(s, dataName = dataName, overWriteS = S)

    #fuck
    return None

def iterateS(dataName, eOri = None, works = True, overWriteLambda_q = None):
    """
    Function to take an initial guess of S, and define two manners of calculating the shape of the integral, to get 

    """
    # First define an initial guess for S
    S_init = 0.05

    # Then calculate the PFOL length at midplane
    if overWriteLambda_q == None:
        PFOL = FOLs(dataName, eOri = eOri)[0]
    else:
        PFOL = overWriteLambda_q

    # Now define a function which calculates what the lambda_int value *should be* based upon the relation lambda_int = lambda_q + 1.64*S
    lambda_int = lambda S : PFOL + 1.64*S

    # And now define a function which is the difference between this value, and the calculated value from the integral
    if works == True:
        lambda_difference = lambda S : lambda_int(S) - intFOL_works(S, dataName = dataName, overWriteLambda_q = overWriteLambda_q)
    elif works == False:
        lambda_difference = lambda S : lambda_int(S) - intFOL_worksnot(S, dataName = dataName, overWriteLambda_q = overWriteLambda_q)

    # Now run an root-finding algorithm on this function
    S_found = fsolve(lambda_difference, S_init)[0]

    return S_found

def qsMax(dataName, sMin = -2, sMax = 2, sNum = 500, eOri = None):
    """
    Simple function to call a q(s) function, and return calculate its maximum value.
    
    """
    # First create the underlying s grid for q(s)
    s = np.linspace(sMin, sMax, num = sNum)

    # Then call q(s)
    qs = q(s, dataName, info = False, eOri = eOri)

    # And get its maximum value
    qsMax = np.max(qs)
    return qsMax
