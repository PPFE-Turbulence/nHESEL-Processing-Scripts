# Module import
import sys
import numpy as np
import matplotlib.pyplot as plt

# Import paths file
pathPath = "/home/s194112"
sys.path.append(pathPath)
import paths

# Append scripts folder
sys.path.append(paths.scriptPath)

# Import methods
from dumpProcessing import dumpProcessing


def multiplot(dataName, tindexRange = [100, -1]):
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
    ## n [m^3], Te [eV], Ti [eV], phi [V], u [m/s], Gamma_perp [1/m$^2$$\cdot$s] q_e [J/m$^2$$\cdot$s], q_i [J/m$^2$$\cdot$s], pe [Pa], pi [Pa]
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
    
    ## First create the plot with subplots
    fig, axes = plt.subplots(2, 3, figsize=(15, 10)) # Possibly adjust
    plt.subplots_adjust(hspace=0.5, wspace=0.5)
    ax1 = axes[0,0]
    ax2 = axes[0,1]
    ax3 = axes[0,2]
    ax4 = axes[1,0]
    ax5 = axes[1,1]
    ax6 = axes[1,2]
    
    ## [n] Density-plot
    ax1.plot(xgridLCFS[:-1]*100, n/1e19, label="n", color="black")
    ax1.plot(xgridLCFS[:-1]*100, n_Init/1e19, label="Initial n", linestyle="--", color="black")
    ylower, yupper = ax1.get_ylim()
    ax1.vlines((x_lcfsPos-x_lcfsPos)*100, label="LCFS", ymin = ylower, ymax = yupper, color="green", linestyle = "--")
    ax1.vlines((x_wallPos-x_lcfsPos)*100, label="Wall", ymin = ylower, ymax = yupper, color="olive")
    ax1.set_xlabel("Position relative to LCFS [cm]", fontsize=11)
    ax1.set_ylabel("Plasma Density [e19 m$^{-3}$]", fontsize=11)
    ax1.set_title("Density")
    ax1.grid()
    ax1.legend()
    
    ## [Te & Ti] Temperature-plot
    ax2.plot(xgridLCFS[:-1]*100, Te, label="T$_e$", color="red")
    ax2.plot(xgridLCFS[:-1]*100, Ti, label="T$_i$", color="blue")
    ax2.plot(xgridLCFS[:-1]*100, Te_Init, label="Initial Temperature", linestyle="--", color="black")
    ylower, yupper = ax2.get_ylim()
    ax2.vlines((x_lcfsPos-x_lcfsPos)*100, label="LCFS", ymin = ylower, ymax = yupper, color="green", linestyle = "--")
    ax2.vlines((x_wallPos-x_lcfsPos)*100, label="Wall", ymin = ylower, ymax = yupper, color="olive")
    ax2.set_xlabel("Position relative to LCFS [cm]", fontsize=11)
    ax2.set_ylabel("Plasma Temperature [eV]", fontsize=11)
    ax2.set_title("Temperatures")
    ax2.grid()
    ax2.legend()
    
    ## [pe & pi] pressures-plot
    ax3.plot(xgridLCFS[:-1]*100, pe/1000, label="p$_e$", color="red")
    ax3.plot(xgridLCFS[:-1]*100, pi/1000, label="p$_i$", color="blue")
    ax3.plot(xgridLCFS[:-1]*100, pe_Init/1000, label="Initial Pressure", linestyle="--", color="black")
    ylower, yupper = ax3.get_ylim()
    ax3.vlines((x_lcfsPos-x_lcfsPos)*100, label="LCFS", ymin = ylower, ymax = yupper, color="green", linestyle = "--")
    ax3.vlines((x_wallPos-x_lcfsPos)*100, label="Wall", ymin = ylower, ymax = yupper, color="olive")
    ax3.set_xlabel("Position relative to LCFS [cm]", fontsize=11)
    ax3.set_ylabel("Plasma Pressure [kPa]", fontsize=11)
    ax3.set_title("Pressures")
    ax3.grid()
    ax3.legend()
    
    ## [u] Velocity plots
    ax4.plot(xgridLCFS[:-1]*100, u[:,1]/1000, label="u$_{||}$", color="brown")
    #ax4.plot(xgridLCFS[:-1]*100, u[:,0]/1000, label="u_rad", linestyle="--", color="black")
    ylower, yupper = ax4.get_ylim()
    ax4.vlines((x_lcfsPos-x_lcfsPos)*100, label="LCFS", ymin = ylower, ymax = yupper, color="green", linestyle = "--")
    ax4.vlines((x_wallPos-x_lcfsPos)*100, label="Wall", ymin = ylower, ymax = yupper, color="olive")
    ax4.set_xlabel("Position relative to LCFS [cm]", fontsize=11)
    ax4.set_ylabel("Velocity [km/s]", fontsize=11)
    ax4.set_title("Velocity")
    ax4.grid()
    ax4.legend()
    
    ## [Gamma_perp & Gamma_para] Perpendicular & parallel particle flux plot
    ax5.plot((xgrid_scrapeoff[:-1])*100, Gamma_perp_scrapeoff, label=r"$\Gamma_\perp$", color="black")
    ax5.plot((xgrid_scrapeoff[:-1])*100, Gamma_para_scrapeoff, label=r"$\Gamma_{||}$", color="green")
    ax5.plot((xgrid_scrapeoff[:-1])*100, Gamma_para_scrapeoff+Gamma_perp_scrapeoff, label=r"Total $\Gamma$", color="purple")
    ylower, yupper = ax5.get_ylim()
    ax5.vlines(0, label="LCFS", ymin = ylower, ymax = yupper, color="green", linestyle = "--")
    ax5.vlines(x_wallPos_scrapeoff*100, label="Wall", ymin = ylower, ymax = yupper, color="olive")
    ax5.set_xlabel("Position relative to LCFS [cm]", fontsize=11)
    ax5.set_ylabel("Plasma Outward Particle Flux [1/m$^2$$\cdot$s]", fontsize=11)
    ax5.set_title("Mass Flux")
    ax5.grid()
    ax5.legend()
    
    ## [qe & qi] Outward q plots
    ax6.plot((xgrid_scrapeoff[:-1])*100, q_perpe_scrapeoff/1000, label=r"q$_{\perp,e}$", color="red")
    ax6.plot((xgrid_scrapeoff[:-1])*100, q_perpi_scrapeoff/1000, label=r"q$_{\perp,i}$", color="blue")
    ax6.plot((xgrid_scrapeoff[:-1])*100, (q_perpe_scrapeoff+q_perpi_scrapeoff)/1000, label=r"Total q$_{\perp}$", color="purple")
    ylower, yupper = ax6.get_ylim()
    ax6.vlines(0, label="LCFS", ymin = ylower, ymax = yupper, color="green", linestyle = "--")
    ax6.vlines(x_wallPos_scrapeoff*100, label="Wall", ymin = ylower, ymax = yupper, color="olive")
    ax6.set_xlabel("Position relative to LCFS [cm]", fontsize=11)
    ax6.set_ylabel("Plasma q$_\perp$ [kJ/m$^2$$\cdot$s]", fontsize=11)
    ax6.set_title("q")
    ax6.grid()
    ax6.legend()

    ## Data title
    plt.suptitle(f"Plots for data: {dataName}",size=24)
    
    ## Show the figure
    # Adjust layout so that subplots do not overlap
    plt.tight_layout()
    
    # Display the figure
    plt.savefig(f"{paths.processedPath}/{dataName}/MULTIPLOT_{dataName}",dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()

    return
