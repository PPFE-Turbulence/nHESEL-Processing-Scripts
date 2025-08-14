### Libraries imports
import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import warnings

### Locals imports
# Import paths file
pathPath = f"{os.environ['HOME']}/paths"
sys.path.append(pathPath)
import paths

# Import scripts path
sys.path.append(paths.scriptPath)

# Import relevant scripts
from dumpProcessing import dumpProcessing

def maxMeansEnergyflux(dataName, maxWidth = 40, lcfsTolerance = 8, electronsOrIons = False, perpOrpara = "perp", eOri = None, seriesMode = False, tindexRange = None):
    """
    Returns the max value for the total energy flux for a given dataset, as well as the mean within a given index range of the max value, and a bool indicating whether the tolerance was broken.
    The total width of this range is maxWidth, and so the mean is taken within a range defined as [maxIndex-maxWidth/2, maxIndex+maxWidth/2].

    Whether the tolerance has been broken (maxVal too close to, or outside of, LCFS) is returned as a bool. False means all good, True means maxVal was too close.
    See maxMeans plot for visualisations. They can be found in all processed folders.

    electronsOrIons can be used to instead get the values for just the electrons or ions, by inputting "e" or "i" instead of False.

    The lcfsTolerance is a measure of how close to the LCFS we wish to exclude data. The default value of 8 means to exclude the data within 8
    index points from the LCFS. If the max index range is for some reason *above* this tolerance, the tolerance is ignored and the full range is used.

    This means a too large tolerance can result in utilizing data from the transistion zone, between the SOL and the layer before the seperatrix.

    """
    ## Interpret perpOrpara
    # Map it
    if type(perpOrpara) == str:
        perpOrpara = perpOrpara.lower()
        perpOrparaDict = {"perp" : 0, "para" : 1}
        perpOrpara = perpOrparaDict[perpOrpara]

    # Check it
    if perpOrpara not in [0, 1]:
        raise Exception(f"perpOrpara is currently set to {perpOrpara}, which is an invalid value. Choose either perp or para, or their corresponding value 0 and 1.")
    
    ## Get some general values
    # Define object
    dumpObj = dumpProcessing(dataName)
    
    # Call methods for values
    # Full range
    lcfsIndex = dumpObj.get_lcfs_index()
    x_lcfsRel, x_lcfsPos = dumpObj.get_lcfs() # LCFS
    x_wallRel, x_wallPos = dumpObj.get_wall() # Wall
    xgrid, zgrid = dumpObj.get_grids() # grids
    
    # Get LCFS-corrected full x-grid
    xgridLCFS = xgrid - x_lcfsPos
    
    # Scrapeoff grid & wall
    xgrid_scrapeoff, zgrid_scrapeoff = dumpObj.get_grids(scrapeoff=True)
    x_wallRel_scrapeoff, x_wallPos_scrapeoff = dumpObj.get_wall(scrapeoff=True)
    
    ## Load up the relevant data
    # Define what to laod
    if perpOrpara == 0:
        nameList = ["q_perpe", "q_perpi"]
    elif perpOrpara == 1:
        nameList = ["q_parae", "q_parai"]

    # For the tindexRange (in case none is given)
    if tindexRange == None:
        tindexRange = [100,-1]

    # Load, taking into account the given indices, and seriesMode
    if seriesMode == False: # Taking the t-average, so we now have [x]
        if tindexRange[1] == -1:
            data = [np.mean(np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_z.npy")[tindexRange[0]:], axis=0) for name in nameList]
            data_Init = [np.squeeze(np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_z_Init.npy"), axis=0) for name in nameList]
            data_scrapeoff = [np.mean(np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_z_scrapeoff.npy")[tindexRange[0]:], axis=0) for name in nameList]
            data_scrapeoff_Init = [np.squeeze(np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_z_scrapeoff_Init.npy"), axis=0) for name in nameList]
        else:
            data = [np.mean(np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_z.npy")[tindexRange[0]:tindexRange[1]], axis=0) for name in nameList]
            data_Init = [np.squeeze(np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_z_Init.npy"), axis=0) for name in nameList]
            data_scrapeoff = [np.mean(np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_z_scrapeoff.npy")[tindexRange[0]:tindexRange[1]], axis=0) for name in nameList]
            data_scrapeoff_Init = [np.squeeze(np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_z_scrapeoff_Init.npy"), axis=0) for name in nameList]
    else: # Not taking the t-average, so we now have [t, x]
        if tindexRange[1] == -1:
            data = [np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_z.npy")[tindexRange[0]:] for name in nameList]
            data_Init = [np.squeeze(np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_z_Init.npy"), axis=0) for name in nameList]
            data_scrapeoff = [np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_z_scrapeoff.npy")[tindexRange[0]:] for name in nameList]
            data_scrapeoff_Init = [np.squeeze(np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_z_scrapeoff_Init.npy"), axis=0) for name in nameList]
        else:
            data = [np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_z.npy")[tindexRange[0]:tindexRange[1]] for name in nameList]
            data_Init = [np.squeeze(np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_z_Init.npy"), axis=0) for name in nameList]
            data_scrapeoff = [np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_z_scrapeoff.npy")[tindexRange[0]:tindexRange[1]] for name in nameList]
            data_scrapeoff_Init = [np.squeeze(np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_z_scrapeoff_Init.npy"), axis=0) for name in nameList]
    
    # Extract the data
    q_perpe, q_perpi = data
    
    # Combine the data
    if eOri != None:
        if eOri == "e":
            q_perp = q_perpe
        elif eOri == "i":
            q_perp = q_perpi
    else:
        q_perp = q_perpe + q_perpi

    ### Commented out because implemented twice for some ood reason?
    # # electronsOrIons code
    # if electronsOrIons != False:
    #     if electronsOrIons == "e":
    #         q_perp = q_perpe
    #     elif electronsOrIons == "i":
    #         q_perp = q_perpi
    #     else:
    #         raise Exception(f"{electronsOrIons} is not a valid input for electronsOrIons. Choose either False, 'e', or 'i'.")
    ###
    
    # Get the max index and value, taking into account the series mode
    if seriesMode == False: # [x]
        maxIndex = np.argmax(q_perp)
        maxVal = np.max(q_perp)
    if seriesMode == True: # [t, x]
        maxIndex = np.argmax(q_perp, axis=1)
        maxVal = np.max(q_perp, axis=1)

    # Create upper and lower maxBounds (for mean)
    if seriesMode == False:
        # First define the lower bound
        maxBounds = [int(maxIndex-maxWidth/2)]
    
        # Then define the upper bound, such that it doesn't mess up and pick part of the gradient when close to the lcfs
        # Note the exception for when the max index is above the lcfs, or too close to it
        toleranceBroken = False
        if maxIndex < lcfsIndex-lcfsTolerance:
            if int(maxIndex+maxWidth/2) > lcfsIndex-lcfsTolerance:
                maxBounds.append(lcfsIndex-5)
            else:
                maxBounds.append(int(maxIndex+maxWidth/2))
        else:
            maxBounds.append(int(maxIndex+maxWidth/2))
            toleranceBroken = True
    else:
        # First define the lower bound for all t
        maxBounds = (maxIndex-maxWidth/2).astype(int)

        # Then the upper bound for all t
        potentialUpperMaxBounds = (maxIndex+maxWidth/2).astype(int)
    
        # Note the exception for when the max index is above the lcfs, or too close to it
        toleranceBroken = False
        uppperBounds = []
        for i in range(len(maxBounds)):
            potentialUpper = potentialUpperMaxBounds[i]
            if potentialUpper < lcfsIndex-lcfsTolerance:
                if potentialUpper > lcfsIndex-lcfsTolerance:
                    uppperBounds.append(lcfsIndex-5)
                else:
                    uppperBounds.append(potentialUpper)
            else:
                uppperBounds.append(potentialUpper)
                toleranceBroken = True
        uppperBounds = np.array(uppperBounds)

        # And append this along the second axis (after adding the second axis to existence for both arrays)
        maxBounds = np.expand_dims(maxBounds, axis=1)
        uppperBounds = np.expand_dims(uppperBounds, axis=1)
        maxBounds = np.append(maxBounds, uppperBounds, axis=1)

    
    if seriesMode == False:
        # Get the mean in a range around the max (+-maxWidth/2 indices - ignoring those close to LCFS)
        maxMean = np.mean(q_perp[maxBounds[0]:maxBounds[1]])
    else:
        # Get the mean in a range around the max (+-maxWidth/2 indices - ignoring those close to LCFS) for each timestep
        maxMean = []
        for i in range(np.shape(q_perp)[0]):
            maxMean.append(np.mean(q_perp[i,maxBounds[i,0]:maxBounds[i,1]]))
        maxMean = np.array(maxMean)

    # NaN value check
    if seriesMode == False:
        if (np.isnan(maxMean) == True) or (np.isnan(maxVal)):
            warnings.warn(f"maxMean or maxVal have NaN value, indicating an issue with loading the processed data. Please check your data folder at {paths.processedPath}/{dataName}. Returning 0 for maxMean and maxVal.")
            maxMean = 0
            maxVal = 0

    return maxVal, maxMean, toleranceBroken

def maxMeansPlot(dataName, maxWidth = 40, lcfsTolerance = 8, perpOrpara = "perp"):
    ## Interpret perpOrpara
    # Map it
    if type(perpOrpara) == str:
        perpOrpara = perpOrpara.lower()
        perpOrparaDict = {"perp" : 0, "para" : 1}
        mode = perpOrparaDict[perpOrpara]

    # Check it
    if mode not in [0, 1]:
        raise Exception(f"perpOrpara is currently set to {perpOrpara}, which is an invalid value. Choose either perp or para, or their corresponding value 0 and 1.")
    
    ## Get some general values
    # Define object
    dumpObj = dumpProcessing(dataName)
    
    # Call methods for values
    # Full range
    lcfsIndex = dumpObj.get_lcfs_index()
    x_lcfsRel, x_lcfsPos = dumpObj.get_lcfs() # LCFS
    x_wallRel, x_wallPos = dumpObj.get_wall() # Wall
    xgrid, zgrid = dumpObj.get_grids() # grids
    
    # Get LCFS-corrected full x-grid
    xgridLCFS = xgrid - x_lcfsPos
    
    # Scrapeoff grid & wall
    xgrid_scrapeoff, zgrid_scrapeoff = dumpObj.get_grids(scrapeoff=True)
    x_wallRel_scrapeoff, x_wallPos_scrapeoff = dumpObj.get_wall(scrapeoff=True)
    
    ## Load up the relevant data
    # Define what to lodd
    if mode == 0:
        nameList = ["q_perpe", "q_perpi"]
    elif mode == 1:
        nameList = ["q_parae", "q_parai"]
    tindexRange = [100,-1]

    # Load, taking into account the given indices
    if tindexRange[1] == -1:
        data = [np.mean(np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_z.npy")[tindexRange[0]:], axis=0) for name in nameList]
        data_Init = [np.squeeze(np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_z_Init.npy"), axis=0) for name in nameList]
        data_scrapeoff = [np.mean(np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_z_scrapeoff.npy")[tindexRange[0]:], axis=0) for name in nameList]
        data_scrapeoff_Init = [np.squeeze(np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_z_scrapeoff_Init.npy"), axis=0) for name in nameList]
    else:
        data = [np.mean(np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_z.npy")[tindexRange[0]:tindexRange[1]], axis=0) for name in nameList]
        data_Init = [np.squeeze(np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_z_Init.npy"), axis=0) for name in nameList]
        data_scrapeoff = [np.mean(np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_z_scrapeoff.npy")[tindexRange[0]:tindexRange[1]], axis=0) for name in nameList]
        data_scrapeoff_Init = [np.squeeze(np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_z_scrapeoff_Init.npy"), axis=0) for name in nameList]
    
    # Extract the data
    q_perpe, q_perpi = data
    
    # Change the unit of the data [W/m^2 -> kW/m^2]
    q_perpe /= 1e3
    q_perpi /= 1e3
    q_perp = q_perpe + q_perpi
    
    # Get the max index and value
    maxIndex = np.argmax(q_perp)
    maxVal = np.max(q_perp)

    # Create upper and lower maxBounds (for mean)
    # First define the lower bound
    maxBounds = [int(maxIndex-maxWidth/2)]

    # Then define the upper bound, such that it doesn't mess up and pick part of the gradient when close to the lcfs
    # Note the exception for when the max index is above the lcfs, or too close to it
    if maxIndex < lcfsIndex-lcfsTolerance:
        if int(maxIndex+maxWidth/2) > lcfsIndex-lcfsTolerance:
            maxBounds.append(lcfsIndex-5)
        else:
            maxBounds.append(int(maxIndex+maxWidth/2))
    else:
        maxBounds.append(int(maxIndex+maxWidth/2))
    
    # Get the mean in a range around the max (+-maxWidth/2 indices - ignoring those close to LCFS)
    maxMean = np.mean(q_perp[maxBounds[0]:maxBounds[1]])
    
    # Now plot the data
    plt.plot(q_perpe, label=nameList[0], color="red")
    plt.plot(q_perpi, label=nameList[1], color="blue")
    if mode == 0:
        plt.plot(q_perp, label="Total q_perp", color="black")
    if mode == 1:
        plt.plot(q_perp, label="Total q_para", color="black")
    xmin, xmax = plt.xlim()
    plt.axvline(lcfsIndex, color="r", linestyle="--", label="LCFS")
    plt.scatter(maxIndex, maxVal, color="orange", label="Max Value")
    plt.axvline(maxBounds[0], color="orange", linestyle = "-.", label = "Max-Index Mean Range")
    plt.axvline(maxBounds[1], color="orange", linestyle = "-.")
    plt.axhline(maxMean, (maxBounds[0]-xmin)/(xmax-xmin), (maxBounds[1]-xmin)/(xmax-xmin), color="green", linestyle='solid', label = "Mean Value in Max-Index Range")
    plt.xlabel("Index")
    plt.ylabel("Energy Flux [kW/m^2]")
    if mode == 0:
        plt.title(f"Max and max-mean perp energy flux values used for P_loss\n{dataName}")
    if mode == 1:
        plt.title(f"Max and max-mean parallel energy flux values\n{dataName}")
    plt.legend()
    plt.grid()
    if mode == 0:
        plt.savefig(f"{paths.processedPath}/{dataName}/maxMeanPlotPerp.png", dpi=300, bbox_inches='tight')
    if mode == 1:
        plt.savefig(f"{paths.processedPath}/{dataName}/maxMeanPlotPara.png", dpi=300, bbox_inches='tight')
    plt.close()