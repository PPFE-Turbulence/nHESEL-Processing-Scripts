# Module import
import sys
import numpy as np
import matplotlib.pyplot as plt

# Import paths file
pathPath = os.environ["HOME"]
sys.path.append(pathPath)
import paths

# Append scripts folder
sys.path.append(paths.scriptPath)

# Import methods
from dumpProcessing import dumpProcessing

def autoplots(dataName, tindexRange = [100, -1]):
    ### First load the data from saved numpy files

    ### OLD LOAD USING EVERYTHING ###
    # # Load the data
    # nameList = ["n","Te","Ti","phi","u","Gamma_perp","Gamma_para","q_perpe","q_perpi", "pe", "pi","q_parae","q_parai"]
    # data = [np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_tz.npy") for name in nameList]
    # data_Init = [np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_tz_Init.npy") for name in nameList]
    # data_scrapeoff = [np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_tz_scrapeoff.npy") for name in nameList]
    # data_scrapeoff_Init = [np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_tz_scrapeoff_Init.npy") for name in nameList]
    ### OLD LOAD USING EVERYTHING ###

    ## Load the data
    # Define data to load
    nameList = ["n","Te","Ti","phi","u","Gamma_perp","Gamma_para","q_perpe","q_perpi", "pe", "pi","q_parae","q_parai"]

    # Then load all the data, taking the t-index range and meaning along this range (first axis mean)
    # In case the last index is -1
    if tindexRange[1] == -1:
        data = [np.mean(np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_z.npy")[tindexRange[0]:], axis=0) for name in nameList]
        data_Init = [np.squeeze(np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_z_Init.npy"), axis=0) for name in nameList]
        data_scrapeoff = [np.mean(np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_z_scrapeoff.npy")[tindexRange[0]:], axis=0) for name in nameList]
        data_scrapeoff_Init = [np.squeeze(np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_z_scrapeoff_Init.npy"), axis=0) for name in nameList]
    # In case it isn't
    else:
        data = [np.mean(np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_z.npy")[tindexRange[0]:tindexRange[1]], axis=0) for name in nameList]
        data_Init = [np.squeeze(np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_z_Init.npy"), axis=0) for name in nameList]
        data_scrapeoff = [np.mean(np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_z_scrapeoff.npy")[tindexRange[0]:tindexRange[1]], axis=0) for name in nameList]
        data_scrapeoff_Init = [np.squeeze(np.load(f"{paths.processedPath}/{dataName}/{name}_averaged_z_scrapeoff_Init.npy"), axis=0) for name in nameList]

    # Extract from data variables
    # First the full range
    n, Te, Ti, phi, u, Gamma_perp, Gamma_para, q_perpe, q_perpi, pe, pi, q_parae, q_parai = data
    n_Init, Te_Init, Ti_Init, phi_Init, u_Init, Gamma_Init, Gamma_para_Init, q_perpe_Init, q_perpi_Init, pe_Init, pi_Init, q_parae_Init, q_parai_Init = data_Init

    # Then the scrapeoff range
    n_scrapeoff, Te_scrapeoff, Ti_scrapeoff, phi_scrapeoff, u_scrapeoff, Gamma_perp_scrapeoff, Gamma_para_scrapeoff, q_perpe_scrapeoff, q_perpi_scrapeoff, pe_scrapeoff, pi_scrapeoff, q_parae_scrapeoff, q_parai_scrapeoff = data_scrapeoff
    n_scrapeoff_Init, Te_scrapeoff_Init, Ti_scrapeoff_Init, phi_scrapeoff_Init, u_scrapeoff_Init, Gamma_scrapeoff_Init, Gamma_para_scrapeoff_Init, q_perpe_scrapeoff_Init, q_perpi_scrapeoff_Init, pe_scrapeoff_Init, pi_scrapeoff_Init, q_parae_scrapeoff_Init, q_parai_scrapeoff_Init = data_scrapeoff_Init
    
    ### Now make time-meaned profiles for the following:
    ## n [m^3], Te [eV], Ti [eV], phi [V], u [m/s], Gamma_perp [1/m^2*s] q_e [J/m^2*s], q_i [J/m^2*s], pe [Pa], pi [Pa]
    ## First get the grid values, and the LCFS value
    # Define object
    dumpObj = dumpProcessing(dataName)
    
    # Call methods for values
    # Full range
    x_lcfsRel, x_lcfsPos = dumpObj.get_lcfs() # LCFS
    x_wallRel, x_wallPos = dumpObj.get_wall() # Wall
    xgrid, zgrid = dumpObj.get_grids() # grids

    # Get LCFS-corrected full x-grid
    xgridLCFS = xgrid - x_lcfsPos

    # Scrapeoff grid & wall
    xgrid_scrapeoff, zgrid_scrapeoff = dumpObj.get_grids(scrapeoff=True)
    x_wallRel_scrapeoff, x_wallPos_scrapeoff = dumpObj.get_wall(scrapeoff=True)

    ### Now plots
    ## [n] Density-plot
    fig, ax = plt.subplots()
    plt.plot(xgridLCFS[:-1], n, label="n", color="black")
    plt.plot(xgridLCFS[:-1], n_Init, label="Initial n", linestyle="--", color="black")
    ylower, yupper = ax.get_ylim()
    plt.vlines(x_lcfsPos-x_lcfsPos, label="LCFS", ymin = ylower, ymax = yupper, color="green", linestyle = "--")
    plt.vlines(x_wallPos-x_lcfsPos, label="Wall", ymin = ylower, ymax = yupper, color="olive")
    plt.xlabel("Position relative to LCFS [m]")
    plt.ylabel("Plasma Density [m^-3]")
    plt.title(f"Space- and Time-averaged Plasma Density\nData: {dataName}")
    plt.grid()
    plt.legend()
    plt.savefig(f"{paths.processedPath}/{dataName}/density_{dataName}")
    plt.close()
    
    ## [phi] Potential-plot
    fig, ax = plt.subplots()
    plt.plot(xgridLCFS[:-1], phi, label="phi", color="black")
    plt.plot(xgridLCFS[:-1], phi_Init, label="Initial phi", linestyle="--", color="black")
    ylower, yupper = ax.get_ylim()
    plt.vlines(x_lcfsPos-x_lcfsPos, label="LCFS", ymin = ylower, ymax = yupper, color="green", linestyle = "--")
    plt.vlines(x_wallPos-x_lcfsPos, label="Wall", ymin = ylower, ymax = yupper, color="olive")
    plt.xlabel("Position relative to LCFS [m]")
    plt.ylabel("Potential [V]")
    plt.title(f"Space- and Time-averaged Potential\nData: {dataName}")
    plt.grid()
    plt.legend()
    plt.savefig(f"{paths.processedPath}/{dataName}/potential_{dataName}")
    plt.close()
    
    ## [u] Velocity plots
    fig, ax = plt.subplots()
    plt.plot(xgridLCFS[:-1], u[:,1]/1000, label="u_pol", color="brown")
    plt.plot(xgridLCFS[:-1], u[:,0]/1000, label="u_rad", linestyle="-", color="black")
    ylower, yupper = ax.get_ylim()
    plt.vlines(x_lcfsPos-x_lcfsPos, label="LCFS", ymin = ylower, ymax = yupper, color="green", linestyle = "--")
    plt.vlines(x_wallPos-x_lcfsPos, label="Wall", ymin = ylower, ymax = yupper, color="olive")
    plt.xlabel("Position relative to LCFS [m]")
    plt.ylabel("Velocity [km/s]")
    plt.title(f"Space- and Time-averaged Plasma Particle Velocity\nData: {dataName}")
    plt.grid()
    plt.legend()
    plt.savefig(f"{paths.processedPath}/{dataName}/velocity_{dataName}")
    plt.close()
    
    ## [Te & Ti] Temperature-plot
    fig, ax = plt.subplots()
    plt.plot(xgridLCFS[:-1], Te, label="Te", color="red")
    plt.plot(xgridLCFS[:-1], Ti, label="Ti", color="blue")
    plt.plot(xgridLCFS[:-1], Te_Init, label="Initial Temperature", linestyle="--", color="black")
    ylower, yupper = ax.get_ylim()
    plt.vlines(x_lcfsPos-x_lcfsPos, label="LCFS", ymin = ylower, ymax = yupper, color="green", linestyle = "--")
    plt.vlines(x_wallPos-x_lcfsPos, label="Wall", ymin = ylower, ymax = yupper, color="olive")
    plt.xlabel("Position relative to LCFS [m]")
    plt.ylabel("Plasma Temperature [eV]")
    plt.title(f"Space- and Time-averaged Plasma Temperatures\nData: {dataName}")
    plt.grid()
    plt.legend()
    plt.savefig(f"{paths.processedPath}/{dataName}/temperatures_{dataName}")
    plt.close()
    
    ## [pe & pi] pressures-plot
    fig, ax = plt.subplots()
    plt.plot(xgridLCFS[:-1], pe/1000, label="pe", color="red")
    plt.plot(xgridLCFS[:-1], pi/1000, label="pi", color="blue")
    plt.plot(xgridLCFS[:-1], pe_Init/1000, label="Initial Pressure", linestyle="--", color="black")
    ylower, yupper = ax.get_ylim()
    plt.vlines(x_lcfsPos-x_lcfsPos, label="LCFS", ymin = ylower, ymax = yupper, color="green", linestyle = "--")
    plt.vlines(x_wallPos-x_lcfsPos, label="Wall", ymin = ylower, ymax = yupper, color="olive")
    plt.xlabel("Position relative to LCFS [m]")
    plt.ylabel("Plasma Pressure [kPa]")
    plt.title(f"Space- and Time-averaged Plasma Pressures\nData: {dataName}")
    plt.grid()
    plt.legend()
    plt.savefig(f"{paths.processedPath}/{dataName}/pressures_{dataName}")
    plt.close()
    
    ## [Gamma_perp & Gamma_para] Outward and parallel particle flux plot
    fig, ax = plt.subplots()
    plt.plot(xgrid_scrapeoff[:-1], Gamma_perp_scrapeoff, label="Perpendicular Gamma", color="black")
    plt.plot(xgrid_scrapeoff[:-1], Gamma_para_scrapeoff, label="Parallel Gamma", color="green")
    plt.plot(xgrid_scrapeoff[:-1], Gamma_para_scrapeoff+Gamma_perp_scrapeoff, label="Total Gamma", color="purple")
    ylower, yupper = ax.get_ylim()
    plt.vlines(0, label="LCFS", ymin = ylower, ymax = yupper, color="green", linestyle = "--")
    plt.vlines(x_wallPos_scrapeoff, label="Wall", ymin = ylower, ymax = yupper, color="olive")
    plt.xlabel("Position relative to LCFS [m]")
    plt.ylabel("Plasma Outward Particle Flux [1/m^2*s]")
    plt.title(f"Space- and Time-averaged Plasma Flux\nData: {dataName}")
    plt.grid()
    plt.legend()
    plt.savefig(f"{paths.processedPath}/{dataName}/particleFlux_{dataName}")
    plt.close()

    ## [Gamma_perp & Gamma_para] Outward and parallel particle flux plot for full simulation domain
    fig, ax = plt.subplots()
    plt.plot(xgridLCFS[:-1], Gamma_perp, label="Perpendicular Gamma", color="black")
    plt.plot(xgridLCFS[:-1], Gamma_para, label="Parallel Gamma", color="green")
    plt.plot(xgridLCFS[:-1], Gamma_para+Gamma_perp, label="Total Gamma", color="purple")
    ylower, yupper = ax.get_ylim()
    plt.vlines(x_lcfsPos-x_lcfsPos, label="LCFS", ymin = ylower, ymax = yupper, color="green", linestyle = "--")
    plt.vlines(x_wallPos-x_lcfsPos, label="Wall", ymin = ylower, ymax = yupper, color="olive")
    plt.xlabel("Position relative to LCFS [m]")
    plt.ylabel("Plasma Outward Particle Flux [1/m^2*s]")
    plt.title(f"Space- and Time-averaged Plasma Flux\nData: {dataName}")
    plt.grid()
    plt.legend()
    plt.savefig(f"{paths.processedPath}/{dataName}/particleFlux_FullDomain_{dataName}")
    plt.close()

    ## [Gamma_perp] Outward particle flux plot
    fig, ax = plt.subplots()
    plt.plot(xgrid_scrapeoff[:-1], Gamma_perp_scrapeoff, label="Perpendicular Gamma", color="black")
    ylower, yupper = ax.get_ylim()
    plt.vlines(0, label="LCFS", ymin = ylower, ymax = yupper, color="green", linestyle = "--")
    plt.vlines(x_wallPos_scrapeoff, label="Wall", ymin = ylower, ymax = yupper, color="olive")
    plt.xlabel("Position relative to LCFS [m]")
    plt.ylabel("Plasma Outward Particle Flux [1/m^2*s]")
    plt.title(f"Space- and Time-averaged Outward Plasma Flux\nData: {dataName}")
    plt.grid()
    plt.legend()
    plt.savefig(f"{paths.processedPath}/{dataName}/particleFlux_PerpOnly_{dataName}")
    plt.close()

    ## [Gamma_perp] Outward particle flux plot for full simulation domain
    fig, ax = plt.subplots()
    plt.plot(xgridLCFS[:-1], Gamma_perp, label="Perpendicular Gamma", color="black")
    ylower, yupper = ax.get_ylim()
    plt.vlines(x_lcfsPos-x_lcfsPos, label="LCFS", ymin = ylower, ymax = yupper, color="green", linestyle = "--")
    plt.vlines(x_wallPos-x_lcfsPos, label="Wall", ymin = ylower, ymax = yupper, color="olive")
    plt.xlabel("Position relative to LCFS [m]")
    plt.ylabel("Plasma Outward Particle Flux [1/m^2*s]")
    plt.title(f"Space- and Time-averaged Outward Plasma Flux\nData: {dataName}")
    plt.grid()
    plt.legend()
    plt.savefig(f"{paths.processedPath}/{dataName}/particleFlux_PerpOnly_FullDomain_{dataName}")
    plt.close()
    
    ## [qe & qi] Outward energy flux plots
    fig, ax = plt.subplots()
    plt.plot(xgrid_scrapeoff[:-1], q_perpe_scrapeoff/1000, label="q_perpe", color="red")
    plt.plot(xgrid_scrapeoff[:-1], q_perpi_scrapeoff/1000, label="q_perpi", color="blue")
    plt.plot(xgrid_scrapeoff[:-1], (q_perpi_scrapeoff+q_perpe_scrapeoff)/1000, label="q_perptotal", color="purple")
    ylower, yupper = ax.get_ylim()
    plt.vlines(0, label="LCFS", ymin = ylower, ymax = yupper, color="green", linestyle = "--")
    plt.vlines(x_wallPos_scrapeoff, label="Wall", ymin = ylower, ymax = yupper, color="olive")
    plt.xlabel("Position relative to LCFS [m]")
    plt.ylabel("Plasma Outward Energy Flux [kJ/m^2*s]")
    plt.title(f"Space- and Time-averaged Perpendicular Plasma Energy Flux\nData: {dataName}")
    plt.grid()
    plt.legend()
    plt.savefig(f"{paths.processedPath}/{dataName}/energyFlux_Perpendicular_{dataName}")
    plt.close()

    ## [qe & qi] Outward energy flux plots for full domain
    fig, ax = plt.subplots()
    plt.plot(xgridLCFS[:-1], q_perpe/1000, label="q_perpe", color="red")
    plt.plot(xgridLCFS[:-1], q_perpi/1000, label="q_perpi", color="blue")
    plt.plot(xgridLCFS[:-1], (q_perpi+q_perpe)/1000, label="q_perptotal", color="purple")
    ylower, yupper = ax.get_ylim()
    plt.vlines(x_lcfsPos-x_lcfsPos, label="LCFS", ymin = ylower, ymax = yupper, color="green", linestyle = "--")
    plt.vlines(x_wallPos-x_lcfsPos, label="Wall", ymin = ylower, ymax = yupper, color="olive")
    plt.xlabel("Position relative to LCFS [m]")
    plt.ylabel("Plasma Outward Energy Flux [kJ/m^2*s]")
    plt.title(f"Space- and Time-averaged Perpendicular Plasma Energy Flux\nData: {dataName}")
    plt.grid()
    plt.legend()
    plt.savefig(f"{paths.processedPath}/{dataName}/energyFlux_Perpendicular_FullDomain_{dataName}")
    plt.close()

    ## [q_||e & q_||i] Parallel (z-direction - out of plane) energy flow
    fig, ax = plt.subplots()
    plt.plot(xgrid_scrapeoff[:-1], q_parae_scrapeoff/1e6, label="q_parae", color="red")
    plt.plot(xgrid_scrapeoff[:-1], q_parai_scrapeoff/1e6, label="q_parai", color="blue")
    plt.plot(xgrid_scrapeoff[:-1], (q_parai_scrapeoff+q_parae_scrapeoff)/1e6, label="q_paratotal", color="purple")
    ylower, yupper = ax.get_ylim()
    plt.vlines(0, label="LCFS", ymin = ylower, ymax = yupper, color="green", linestyle = "--")
    plt.vlines(x_wallPos_scrapeoff, label="Wall", ymin = ylower, ymax = yupper, color="olive")
    plt.xlabel("Position relative to LCFS [m]")
    plt.ylabel("Plasma Outward Energy Flux [MJ/m^2*s]")
    plt.title(f"Space- and Time-averaged Parallel Plasma Energy Flux\nData: {dataName}")
    plt.grid()
    plt.legend()
    plt.savefig(f"{paths.processedPath}/{dataName}/energyFlux_Parallel_{dataName}")
    plt.close()

    ## [q_||e & q_||i] Parallel (z-direction - out of plane) energy flow for full range
    fig, ax = plt.subplots()
    plt.plot(xgridLCFS[:-1], q_parae/1e6, label="q_parae", color="red")
    plt.plot(xgridLCFS[:-1], q_parai/1e6, label="q_parai", color="blue")
    plt.plot(xgridLCFS[:-1], (q_parai+q_parae)/1e6, label="q_paratotal", color="purple")
    ylower, yupper = ax.get_ylim()
    plt.vlines(x_lcfsPos-x_lcfsPos, label="LCFS", ymin = ylower, ymax = yupper, color="green", linestyle = "--")
    plt.vlines(x_wallPos-x_lcfsPos, label="Wall", ymin = ylower, ymax = yupper, color="olive")
    plt.xlabel("Position relative to LCFS [m]")
    plt.ylabel("Plasma Outward Energy Flux [MJ/m^2*s]")
    plt.title(f"Space- and Time-averaged Parallel Plasma Energy Flux\nData: {dataName}")
    plt.grid()
    plt.legend()
    plt.savefig(f"{paths.processedPath}/{dataName}/energyFlux_Parallel_FullDomain_{dataName}")
    plt.close()

    return 0