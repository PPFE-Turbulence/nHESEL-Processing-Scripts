# Module import
import sys
import numpy as np
import matplotlib.pyplot as plt
import warnings

# Import paths file
pathPath = "/home/s194112"
sys.path.append(pathPath)
import paths

# Append scripts folder
sys.path.append(paths.scriptPath)

# Function to process data
from heselProcessing import heselProcessing
from dumpProcessing import dumpProcessing

# Userwarning supression to avoid logspam
warnings.filterwarnings("ignore", category=UserWarning)


def generateLabels(data, ax, xORy, unit = None):
    """
    Function used to generate a list of tickmark labels for ax.set_tickmarks, from a list of data points. The unit gets appended at the end of the string.

    The labels to input into ax.set_tickmarks() is returned from this function.
    """
    
    # First get the current tickmarks
    if xORy == "x":
        labels = ax.get_xticklabels()
    elif xORy == "y":
        labels = ax.get_yticklabels()

    # Then get the text for the tickmarks
    labels = [item.get_text() for item in labels]

    # Remove the first and last label because they aren't used for some odd reason
    labels = labels[1:-1]

    # Now convert the labels to floats, and use this to figure out which 'percentage' of the interval we are at (normalise the interval)
    floats = np.array([float(string.replace('−', '-')) for string in labels]) # Wrong minus sign usage fix
    floats = floats - np.min(floats) # Ensure we have the right starting point
    perc = floats/np.max(floats)

    # Now extract the closest values to these percentages, by multiplying the percentages with the amount of values in the data array/list, and either flooring or ceiling them (rounding)
    scaled_percs = perc*(len(data)-1) # -1 due to counting from 0
    ints = [int(np.floor(val)) for val in scaled_percs]

    # Finally use the found integers to create the new labels, and return them
    for i in range(len(ints)):
        if unit != None:
            labels[i] = f"{data[ints[i]]:.1f} {unit}"
        else:
            labels[i] = f"{data[ints[i]]:.1f}"

    # Then insert a dummy-label
    labels.insert(0, -1)

    return labels

def densityShows(dataName):
    ## Get the data
    # Create required object
    processingObj = heselProcessing(dataName = dataName)
    
    # Get the required data
    n_averaged = np.load(f"{paths.processedPath}/{dataName}/n_averaged_z.npy")
    n_last = processingObj.plainVariables(varNames = ["n"], returnData = True, saveData = False, tindexRange = [-1,-1], info = False)[0]

    ## Get the grid and do some conversions
    # Create required object
    dumpObj = dumpProcessing(dataName = dataName)

    # Get grid and process
    lcfsrel, lcfsX = dumpObj.get_lcfs()
    xgrid, ygrid = dumpObj.get_grids()
    xgrid = np.flip(xgrid)
    xgrid -= lcfsX# Remove LCFS position from array
    xgrid *= 1e2 # Conversion [m] -> [cm]
    ygrid *= 1e2 # Conversion [m] -> [cm]
    tvals = dumpObj.get_times()
    tvals *= 1e3 # Conversion [s] -> [ms]
    
    ### THE Y-MEANED DATA PLOT ###
    # Plot all time in imshow, meaned along z-axis (y for them)
    fig, ax = plt.subplots()
    ax.imshow(np.rot90(n_averaged))
    #ax.invert_yaxis()
    plt.title(f"{dataName}[n] Averages")
    plt.xlabel("Time [ms]")
    plt.ylabel("Distance Relative to LCFS [cm]")
    
    # We need to draw the canvas, otherwise the labels won't be positioned and won't have values yet
    fig.canvas.draw()
    
    # Get the labels for the tickmarks, and replace them
    labelsx = generateLabels(tvals, ax, "x")#, "ms")
    labelsy = generateLabels(xgrid, ax, "y")#, "mm")
    ax.set_xticklabels(labelsx)
    ax.set_yticklabels(labelsy)
    
    # Make a horisontal line for the lcfs
    ylower, yupper = ax.get_ylim() # Get the upper and lower bound, may be inverted
    ydif = yupper - ylower # Get the total distance between upper and lower bound
    lcfsPos = ylower + ydif*lcfsrel# Determine the position of the lcfs line in the plot, as the scaled distance from ylower
    plt.axhline(lcfsPos, color="r",linestyle = "--", label="LCFS")
    plt.legend()
    
    # Then show the figure
    plt.savefig(f"{paths.processedPath}/{dataName}/DENSITIES_{dataName}",dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()
    
    
    
    ### THE LAST DENSITY DATA PLOT ###
    # Plot all time in imshow, meaned along z-axis (y for them)
    fig, ax = plt.subplots()
    ax.imshow(n_last)
    ax.invert_yaxis()
    plt.title(f"{dataName}[n] Final")
    plt.xlabel("y [cm]")
    plt.ylabel("Relative distance from LCFS [cm]")
    
    # We need to draw the canvas, otherwise the labels won't be positioned and won't have values yet
    fig.canvas.draw()
    
    # Get the labels for the tickmarks, and replace them
    labelsx = generateLabels(ygrid, ax, "x")#, "mm")
    labelsy = generateLabels(np.flip(xgrid), ax, "y")#, "mm")
    ax.set_xticklabels(labelsx)
    ax.set_yticklabels(labelsy)
    
    # Make a horisontal line for the lcfs
    ylower, yupper = ax.get_ylim() # Get the upper and lower bound, may be inverted
    ydif = yupper - ylower # Get the total distance between upper and lower bound
    lcfsPos = ylower + ydif*lcfsrel# Determine the position of the lcfs line in the plot, as the scaled distance from ylower
    plt.axhline(lcfsPos, color="r",linestyle = "--", label="LCFS")
    plt.legend()
    
    # Then show the figure
    plt.savefig(f"{paths.processedPath}/{dataName}/FINALDENSITY_{dataName}",dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()

    return
