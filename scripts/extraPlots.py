### Libraries imports
import sys
import numpy as np
import matplotlib.pyplot as plt

### Locals imports
# Import paths file
pathPath = os.environ["HOME"]
sys.path.append(pathPath)
import paths

# Import scripts folder
sys.path.append(paths.scriptPath)

# Scripts imports
from dumpProcessing import dumpProcessing
from FOLs import FOLs

def qPlotsBasic(dataName, includeFOL = True):
    # Define basic quantities
    minIndex = 100
    savePath = f"{paths.processedPath}/{dataName}"
    
    ## Load data & convert to proper unit quantities (from base SI)
    # First load the relevant quantities
    q_perpe = np.load(f"{paths.processedPath}/{dataName}/q_perpe_averaged_z_scrapeoff.npy")
    q_perpi = np.load(f"{paths.processedPath}/{dataName}/q_perpi_averaged_z_scrapeoff.npy")
    q_parae = np.load(f"{paths.processedPath}/{dataName}/q_parae_averaged_z_scrapeoff.npy")
    q_parai = np.load(f"{paths.processedPath}/{dataName}/q_parai_averaged_z_scrapeoff.npy")
    
    # Combine
    q_perp = q_perpe + q_perpi
    q_para = q_parae + q_parai
    q_total = q_perp + q_para

    # Reduce scale
    q_perp /= 1e6 # (W -> MW)
    q_para /= 1e6 # (W -> MW)
    q_total /= 1e6 # (W -> MW)

    # Pick out the data above 100 nouts
    q_perp = q_perp[minIndex:,:]
    q_para = q_para[minIndex:,:]
    q_total = q_total[minIndex:,:]
    
    # Calculate t-meaned values, and std values
    q_perp_mean = np.mean(q_perp, axis=0)
    q_perp_std = np.std(q_perp, axis=0)
    q_para_mean = np.mean(q_para, axis=0)
    q_para_std = np.std(q_para, axis=0)
    q_total_mean = np.mean(q_total, axis=0)
    q_total_std = np.std(q_total, axis=0)

    # And get the x-grid for the data
    dumpObj = dumpProcessing(dataName = dataName)
    grids = dumpObj.get_grids(scrapeoff = True)
    xgrid = grids[0][1:]

    # Convert to cm scale
    xgrid *= 100

    ## Also get the power fall-off length
    PFOL = FOLs(dataName = dataName)[0]
    #PFOL = FOLs(dataName = dataName, average = True)[0][0]

    # And remember to convert it from m to cm
    PFOL *= 100

    ## After loading and processing the data, it is time to plot
    # First q plot - perpendicular
    plt.figure()
    plt.plot(xgrid, q_perp_mean, color = "r", label = r"<q$_\perp$>$_{t,z}$")
    plt.plot(xgrid, q_perp_mean + 2*q_perp_std, color = "r", linestyle = "--", label = "95% Confidence Interval")
    plt.plot(xgrid, q_perp_mean - 2*q_perp_std, color = "r", linestyle = "--")
    plt.title(f"q$_\perp$\nNout > {minIndex}\nData: {dataName}")
    plt.ylabel(r"q$_\perp$ [MW/m$^2$]", fontsize=11)
    plt.xlabel("Position in SOL [cm]", fontsize=11)
    plt.grid()
    plt.legend()
    plt.savefig(f"{savePath}/q_perp_plot_{dataName}.png", dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()

    # Second q plot - parallel
    plt.figure()
    plt.plot(xgrid, q_para_mean, color = "b", label = r"<q$_{||}$>$_{t,z}$")
    plt.plot(xgrid, q_para_mean + 2*q_para_std, color = "b", linestyle = "--", label = "95% Confidence Interval")
    plt.plot(xgrid, q_para_mean - 2*q_para_std, color = "b", linestyle = "--")
    if includeFOL == True:# and True == False: # Disabling this temporarily # Re-enabling this
        plt.vlines(PFOL, ymin = np.min(q_para_mean - 2*q_para_std), ymax = np.max(q_para_mean + 2*q_para_std), color = "black", linestyle = "-.", label = "Power Fall-off Length")
        plt.title(f"q$_{||}$ & Power Fall-off Length\nNout > {minIndex}\nData: {dataName}")
    else:
        plt.title(f"q$_{||}$\nNout > {minIndex}\nData: {dataName}")
    plt.ylabel("q$_{||}$ [MW/m$^2$]", fontsize=11)
    plt.xlabel("Position in SOL [cm]", fontsize=11)
    plt.grid()
    plt.legend()
    plt.savefig(f"{savePath}/q_para_plot_{dataName}.png", dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()

    # Third q plot - parallel and total
    plt.figure()
    plt.plot(xgrid, q_para_mean, color = "b", label = r"<q$_{||}$>$_{t,z}$")
    plt.plot(xgrid, q_para_mean + 2*q_para_std, color = "b", linestyle = "--", label = "95% Confidence Interval")
    plt.plot(xgrid, q_para_mean - 2*q_para_std, color = "b", linestyle = "--")
    plt.plot(xgrid, q_total_mean, color = "purple", label = r"<q$_{total}$>$_{t,z}$")
    plt.plot(xgrid, q_total_mean + 2*q_total_std, color = "purple", linestyle = "--", label = "95% Confidence Interval")
    plt.plot(xgrid, q_total_mean - 2*q_total_std, color = "purple", linestyle = "--")
    if includeFOL == True:
        plt.vlines(PFOL, ymin = np.min(q_para_mean - 2*q_para_std), ymax = np.max(q_para_mean + 2*q_para_std), color = "black", linestyle = "-.", label = "PFOL")
        # plt.title(f"Parallel and Total q, & Power Fall-off Length\nNout > {minIndex}\nData: {dataName}")
    # else:
        # plt.title(f"Parallel and Total q\nNout > {minIndex}\nData: {dataName}")
    plt.ylabel("Heat Flux [MW/m$^2$]", fontsize=11)
    plt.xlabel("Position in SOL [cm]", fontsize=11)
    plt.grid()
    plt.legend()
    plt.savefig(f"{savePath}/q_total_plot_{dataName}.png", dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()

    return

def GammaPlotsBasic(dataName, includeFOL = True):
    # Define basic quantities
    minIndex = 100
    
    ## Load data & convert to proper unit quantities (from base SI)
    # First load the relevant quantities
    Gamma = np.load(f"{paths.processedPath}/{dataName}/Gamma_perp_averaged_z_scrapeoff.npy")[minIndex:,:] # From 600 and up

    # Reduce scale
    # Not needed

    # Calculate t-meaned values, and std values
    Gamma_mean = np.mean(Gamma, axis=0)
    Gamma_std = np.std(Gamma, axis=0)

    # And get the x-grid for the data
    dumpObj = dumpProcessing(dataName = dataName)
    grids = dumpObj.get_grids(scrapeoff = True)
    xgrid = grids[0][1:]

    # Convert to cm scale
    xgrid *= 100

    ## Also get the density fall-off length
    DFOL = FOLs(dataName = dataName)[1] # Average before takign FOL
    #DFOL = FOLs(dataName = dataName, average = True)[1][0] # Average after taking FOL

    # And remember to convert it from m to cm
    DFOL *= 100

    ## After loading and processing the data, it is time to plot
    # Gamma plot
    plt.figure()
    plt.plot(xgrid, Gamma_mean, color = "black", label = r"<$\Gamma_\perp$>$_{t,z}$")
    plt.plot(xgrid, Gamma_mean + 2*Gamma_std, color = "black", linestyle = "--", label = "95% Confidence Interval")
    plt.plot(xgrid, Gamma_mean - 2*Gamma_std, color = "black", linestyle = "--")
    if includeFOL == True:
        plt.vlines(DFOL, ymin = np.min(Gamma_mean - 2*Gamma_std), ymax = np.max(Gamma_mean + 2*Gamma_std), color = "blue", linestyle = "-.", label = "DFOL")
        plt.title(f"Perpendicular Density Flux & Density Fall-off Length\nNout > {minIndex}\nData: {dataName}")
    else:
        plt.title(f"Perpendicular Density Flux\nNout > {minIndex}\nData: {dataName}")
    plt.ylabel(r"Density Flux [1/m$^2$$\cdot$s]", fontsize=11)
    plt.xlabel("Position in SOL [cm]", fontsize=11)
    plt.grid()
    plt.legend()
    plt.savefig(f"{paths.processedPath}/{dataName}/Gamma_perp_plot_{dataName}.png", dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()

    return

def LCFSPlots(dataName, rollingAvg = True, rollingAvgSize = 10):
    # Define basic quantities
    minIndex = 100
    savePath = f"{paths.processedPath}/{dataName}"
    
    ## Loading & processing the data
    # First load and slice it to get LCFS
    q_perpe_LCFS = np.load(f"{paths.processedPath}/{dataName}/q_perpe_averaged_z_scrapeoff.npy")[:,0] # (indices [t,x])
    q_perpi_LCFS = np.load(f"{paths.processedPath}/{dataName}/q_perpi_averaged_z_scrapeoff.npy")[:,0]
    Gamma_perp_LCFS = np.load(f"{paths.processedPath}/{dataName}/Gamma_perp_averaged_z_scrapeoff.npy")[:,0]

    # Then combine
    q_perp_LCFS = q_perpe_LCFS + q_perpi_LCFS

    # And change units
    q_perp_LCFS /= 1e6 # W -> MW

    # Finally, get the time for the timesteps
    dumpObj = dumpProcessing(dataName = dataName)
    times = dumpObj.get_times()

    # And change units
    times *= 1e3 # s -> ms

    # Get the time the minIndex corresponds to
    minTime = ((minIndex-25)/np.shape(q_perpe_LCFS)[0])*times[-1]
    usualTime = (minIndex/np.shape(q_perpe_LCFS)[0])*times[-1] # minIndex/np.shape(q_perpe_LCFS)[0]) scales the index according to the amount of total indices, times[-1] then multiplies this by the total time
    
    ## Now plot
    # First the q_perp plot (what energy is crossing the seperatrix)
    plt.figure()
    plt.plot(times, q_perp_LCFS, label = "Average", color = "orange")
    plt.axvline(x = minTime, color = 'red', linestyle = "-.", label = 'Minimum Index Cutoff')
    plt.axvline(x = usualTime, color = 'b', linestyle = "--", label = 'Usual Index Cutoff')
    plt.xlabel("Time [ms]", fontsize=11)
    plt.ylabel("Perpendicular Heat Flux [MW/m$^2$]", fontsize=11)
    plt.title(f"Average q Across LCFS\nData: {dataName}")
    if rollingAvg == True:
        filtered = np.convolve(q_perp_LCFS, np.ones(rollingAvgSize)/rollingAvgSize, mode='same')
        # Edge effects 'fix'
        filtered[:int(rollingAvgSize/2)] = filtered[int(rollingAvgSize/2)]
        filtered[-int(rollingAvgSize/2):] = filtered[-int(rollingAvgSize/2)]
        plt.plot(times, filtered, label = f"N = {rollingAvgSize} Rolling Average of Average", color="black")
        plt.legend()
    plt.grid()
    plt.savefig(f"{savePath}/LCFS_q(t)_{dataName}.png", dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()

    # Then the Gamma_perp plot (what mass is crossing the seperatrix)
    plt.figure()
    plt.plot(times, Gamma_perp_LCFS, label = "Average", color = "purple")
    plt.axvline(x = minTime, color = 'red', linestyle = "-.", label = 'Minimum Index Cutoff')
    plt.axvline(x = usualTime, color = 'b', linestyle = "--", label = 'Usual Index Cutoff')
    plt.xlabel("Time [ms]", fontsize=11)
    plt.ylabel("Perpendicular Density Flux [1/m$^2$$\cdot$s]", fontsize=11)
    plt.title(f"Average Density Flux Across LCFS\nData: {dataName}")
    if rollingAvg == True:
        filtered = np.convolve(Gamma_perp_LCFS, np.ones(rollingAvgSize)/rollingAvgSize, mode='same')
        # Edge effects 'fix'
        filtered[:int(rollingAvgSize/2)] = filtered[int(rollingAvgSize/2)]
        filtered[-int(rollingAvgSize/2):] = filtered[-int(rollingAvgSize/2)]
        plt.plot(times, filtered, label = f"N = {rollingAvgSize} Rolling Average of Average", color = "black")
        plt.legend()
    plt.grid()
    plt.savefig(f"{savePath}/LCFS_Gamma(t)_{dataName}.png", dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()
    
    return
