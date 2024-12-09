import xarray as xr
import numpy as np
import numpy.matlib

import time
import math

from scipy.interpolate import interp1d, griddata
from scipy.special import gammainc
from scipy.stats.distributions import chi2
from scipy.signal import welch

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.cm import ScalarMappable
from matplotlib.ticker import ScalarFormatter
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

import cmocean

g = 9.81 #  acceleration of gravity
creator_name = 'g.marechal'



def compute_bulk_l2s_SWIM(ds, id_start):
    
    """
    Purpose: prepare the directional wave spectrum observed with SWIM (L2S product) and compute the dominant wave direction/wavenumber and wave Bulk
    ---------
    Inputs: 
    ---------
    The Xarray Dataset (format: netcdf file)
    id_start: the index of the first wavenumber considered (to split partitions)

    Outputs:
    ---------
    theta_swim: The direction of the wave in North-East frame of coordinate
    sub_k: The wavenumber axis
    e_kth: the wavenumber-direction spectrum
    Lp: the dominant wavelength
    Dp: the dominant wavelength
    sigp: The directional spreading
    """

    

    sub_k = ds.k[id_start:].values # Remove short waves
    
    dk = np.gradient(sub_k) # the wavenumber step
    
    
    e_kth = (ds.wave_spectra[len(ds.time.values)//2:, id_start:]*(sub_k)**(-1)).values.T# the amplitude spec. for half a spectrum
    #e_kth = (ds.wave_spectra[:, id_start:]*(sub_k)**(-1)).values.T

    #theta_swim = ds.phi_geo[:].values * np.pi/180
    theta_swim = ds.phi_geo[len(ds.time.values)//2:].values * np.pi/180 # the direction projected onto N-E frame of coord.

    
    dth = abs(np.diff(theta_swim)[0]) # The directional step

    #print(f'Hs = {1/2*4*np.sqrt(np.nansum(np.nansum(e_kth.T*(sub_k)**(-2)*dk, axis = 1), axis= 0))*dth}')
    Hsnew = 4*np.sqrt(np.nansum(np.nansum(e_kth.T*dk, axis = 1), axis = 0)*dth)
    Etot = Hsnew**2/16

    #actual_ekth = e_kth.T*(sub_k)**(-2)
    
    e_kth = e_kth.T
    #norm_fac = e_kth.T*(sub_k)**(-2)
    ek = np.nansum(e_kth, axis = 0) * dth # The omnidirectional wave spectrum


    ##########
    # ---  Switch into frequency-direction spectrum
    ##########    
    g = 9.81
    cg = 1/2 * np.sqrt(g/sub_k) # The group velocity

    e_fth = (2*np.pi)/cg * e_kth

    ef = np.nansum(e_fth, axis = 0) * dth
    fren = (g * sub_k)**(1/2)/(2*np.pi)

    ##########
    # ---  Switch into frequency-direction spectrum
    ##########    
    g = 9.81
    cg = 1/2 * np.sqrt(g/sub_k) # The group velocity

    e_fth = (2*np.pi)/cg * e_kth

    ef = np.nansum(e_fth, axis = 0) * dth
    fren = (g * sub_k)**(1/2)/(2*np.pi)


    #########
    # --- Compute the 2D Peakedness param
    #########
    Q1 = np.zeros((np.size(fren)))
    dfreq = np.gradient(fren)

    Q2=0
    Etot =  np.nansum(np.nansum(e_fth * dfreq, axis = 1) * dth)

    #Etot0 =  np.nansum(np.nansum(e_fth_amp * dfreq, axis = 1) * dth)
    Q1 = np.zeros((len(fren)))
    for ind in range(np.size(fren)):
        Q1[ind]=np.sum(e_fth[:,ind]**2)*dth
        Q2=Q2+Q1[ind]*dfreq[ind]*g**2/(2*((np.pi*2)**4*fren[ind]**3 ))
    
        #   print('ft:',ind,fren[ind],wn[ind],'Cg:',Cg[ind],dfreq[ind]/(wn[ind]*dk[ind]),grav**2/(2*((np.pi*2)**4*fren[ind]**3 )) )
    Qkk=np.sqrt(Q2)/Etot

    print('Qkk =', Qkk, 'm')

    idp_th, idp_k = np.where(e_kth == np.amax(e_kth)) # The index of the dominant wavenumber/direction
    Lp = 2*np.pi/sub_k[idp_k][0] # The dominant wavelength
    Dp = theta_swim[idp_th] # The dominant direction
    print('Peak Wavelength =', Lp, 'm')
    print('Peak Direction =', Dp[0] * 180/np.pi, 'deg')

    ##########
    # ---  The directional spreading
    ##########

    C11 = np.nansum(e_fth, axis = 0) * dth
    #print(np.nansum((sub_k**2*e_fth).T*np.cos(theta_swim)))
    C22 = np.nansum((sub_k**2*e_fth).T*(np.cos(theta_swim))**2, axis = 1) * dth
    C33 = np.nansum((sub_k**2*e_fth).T*(np.sin(theta_swim))**2, axis = 1) * dth
    Q12 = np.nansum((sub_k*e_fth).T*np.cos(theta_swim), axis = 1) * dth
    Q13 = np.nansum((sub_k*e_fth).T*np.sin(theta_swim), axis = 1) * dth
       
    a = Q12/(np.sqrt((C22 + C33)*C11))
    b = Q13/(np.sqrt((C22 + C33)*C11))
    
    the_m = np.arctan2(b, a)
    #the_m = np.nansum(e_kth*theta_swim, axis = 1)/(np.nansum(e_kth, axis = 1)) # Works for Buoy
    sig = []                                 
    for k in range(len(sub_k)):                                  
        #sig.append((np.nansum(e_kth[k, :] * (theta_swim - the_m[k])**(2), axis = 0)/(np.nansum(e_kth[k,:], axis = 0)))**(1/2)) # works for Buoy
        sig.append(np.sqrt(2 * (1-np.sqrt(a[k]**2 + b[k]**2))))
        

    Qp = (np.sum(ef*dfreq))**2/(np.sum(ef**2*dfreq)) 

    #print(np.sqrt(a[k]**2 + b[k]**2))

    sig = np.array(sig)
    idp = np.where(ek==np.nanmax(ek))[0][0]

    sigp = sig[idp]*180/np.pi
    the_mp = the_m[idp] * 180/np.pi
    #plt.figure()
    #plt.plot(sub_k, sig * 180/np.pi)
    
    return theta_swim, sub_k, e_kth, Lp, Dp[0]*180/np.pi, the_mp, sigp, Qp, Qkk

def compute_bulk_l2s_SWIM_group(ds, id_start):
    
    """
    Purpose: prepare the directional wave spectrum observed with SWIM (L2S product)and compute the dominant wave direction/wavenumber and wave Bulk. This function is slightly modify for the wave group cases (short waves are removed)
    ---------
    Inputs: 
    ---------
    The Xarray Dataset (format: netcdf file)
    id_start: the index of the last wavenumber considered (to split partitions)

    Outputs:
    ---------
    theta_swim: The direction of the wave in North-East frame of coordinate
    sub_k: The wavenumber axis
    e_kth: the wavenumber-direction spectrum
    Lp: the dominant wavelength
    Dp: the dominant wavelength
    sigp: The directional spreading
    """
    

    sub_k = ds.k[:id_start].values # Remove short waves
    
    dk = np.gradient(sub_k) # the wavenumber step
    
    
    e_kth = 2*(ds.wave_spectra[len(ds.time.values)//2:, :id_start]*(sub_k)**(-1)).values.T# the amplitude spec. for half a spectrum
    theta_swim = ds.phi_geo[len(ds.time.values)//2:].values * np.pi/180 # the direction projected onto N-E frame of coord.
    #plt.figure()
    #plt.plot(sub_k, np.nansum(e_kth, axis = 1), '-o')
    
    dth = abs(np.diff(theta_swim)[0]) # The directional step

    #print(f'Hs = {1/2*4*np.sqrt(np.nansum(np.nansum(e_kth.T*(sub_k)**(-2)*dk, axis = 1), axis= 0))*dth}')
    Hsnew = 4*np.sqrt(np.nansum(np.nansum(e_kth.T*dk, axis = 1), axis = 0)*dth)
    Etot = Hsnew**2/16

    #actual_ekth = e_kth.T*(sub_k)**(-2)
    
    e_kth = e_kth.T
    #norm_fac = e_kth.T*(sub_k)**(-2)
    ek = np.nansum(e_kth, axis = 0) * dth # The omnidirectional wave spectrum
    ##########
    # ---  Switch into frequency-direction spectrum
    ##########    
    g = 9.81
    cg = 1/2 * np.sqrt(g/sub_k) # The group velocity

    e_fth = (2*np.pi)/cg * e_kth

    ef = np.nansum(e_fth, axis = 0) * dth
    fren = (g * sub_k)**(1/2)/(2*np.pi)

    ##########
    # ---  Switch into frequency-direction spectrum
    ##########    
    g = 9.81
    cg = 1/2 * np.sqrt(g/sub_k) # The group velocity

    e_fth = (2*np.pi)/cg * e_kth

    ef = np.nansum(e_fth, axis = 0) * dth
    fren = (g * sub_k)**(1/2)/(2*np.pi)



    #########
    # --- Compute the 2D Peakedness param
    #########
    Q1 = np.zeros((np.size(fren)))
    dfreq = np.gradient(fren)

    Q2=0
    Etot =  np.nansum(np.nansum(e_fth * dfreq, axis = 1) * dth)

    #Etot0 =  np.nansum(np.nansum(e_fth_amp * dfreq, axis = 1) * dth)
    Q1 = np.zeros((len(fren)))
    for ind in range(np.size(fren)):
        Q1[ind]=np.sum(e_fth[:,ind]**2)*dth
        Q2=Q2+Q1[ind]*dfreq[ind]*g**2/(2*((np.pi*2)**4*fren[ind]**3 ))
        #   print('ft:',ind,fren[ind],wn[ind],'Cg:',Cg[ind],dfreq[ind]/(wn[ind]*dk[ind]),grav**2/(2*((np.pi*2)**4*fren[ind]**3 )) )
    Qkk=np.sqrt(Q2)/Etot
    print('Qkk =', Qkk)


    idp_th, idp_k = np.where(e_kth == np.amax(e_kth)) # The index of the dominant wavenumber/direction
    Lp = 2*np.pi/sub_k[idp_k][0] # The dominant wavelength
    Dp = theta_swim[idp_th] # The dominant direction
    
    ##########
    # ---  The directional spreading
    ##########

    C11 = np.nansum(e_fth, axis = 0) * dth
    #print(np.nansum((sub_k**2*e_fth).T*np.cos(theta_swim)))
    C22 = np.nansum((sub_k**2*e_fth).T*(np.cos(theta_swim))**2, axis = 1) * dth
    C33 = np.nansum((sub_k**2*e_fth).T*(np.sin(theta_swim))**2, axis = 1) * dth
    Q12 = np.nansum((sub_k*e_fth).T*np.cos(theta_swim), axis = 1) * dth
    Q13 = np.nansum((sub_k*e_fth).T*np.sin(theta_swim), axis = 1) * dth
       
    a = Q12/(np.sqrt((C22 + C33)*C11))
    b = Q13/(np.sqrt((C22 + C33)*C11))
    
    the_m = np.arctan2(b, a)
    #the_m = np.nansum(e_kth*theta_swim, axis = 1)/(np.nansum(e_kth, axis = 1)) # Works for Buoy
    sig = []                                 
    for k in range(len(sub_k)):                                  
        #sig.append((np.nansum(e_kth[k, :] * (theta_swim - the_m[k])**(2), axis = 0)/(np.nansum(e_kth[k,:], axis = 0)))**(1/2)) # works for Buoy
        sig.append(np.sqrt(2 * (1-np.sqrt(a[k]**2 + b[k]**2))))
        

    Qp = (np.sum(ef*dfreq))**2/(np.sum(ef**2*dfreq)) 

    #print(np.sqrt(a[k]**2 + b[k]**2))

    sig = np.array(sig)
    idp = np.where(ek==np.nanmax(ek))[0][0]

    sigp = sig[idp]*180/np.pi
    the_mp = the_m[idp] * 180/np.pi
    #plt.figure()
    #plt.plot(sub_k, sig * 180/np.pi)
    
    return theta_swim, sub_k, e_kth, Lp, Dp[0]*180/np.pi, the_mp, sigp, Qp, Qkk



def plot_polar_spec_SWIM(ax_plot, cax_cbar, theta_SWIM, wavenumber_SWIM, spec_SWIM, vmin0, vmax0, side_spec, path_out, file_out):
    """
    Purpose: Plot the mean d wavenumber-direction wave spectrum observed by SWIM
    ---------
    
    Inputs: The Xarray Dataset (format: netcdf file)
    ---------
    ax_plot, cax_cbar: the axis for the polar plot, need to be defined beforehand
    cax_cbar: the position of the colorbar
    theta_SWIM: the direction of the wave in the north-east frame of coordinate
    wavenumber_SWIM: the considered wavenumbers
    spec_SWIM: the 2-D wavenumber-direction wave spectrum
    vmin0, vmax0: the vmin and vmax for pcolor function
    side_spec: The side of the wave spectrum you want to plot (left or right: depends on the convention you want)
    path_out: The directory where the figure is saved
    file_out: The name of the file

    Outputs: 
    ---------
    empty
    """
    ax = ax_plot
    
    if side_spec=='right':
        
        p1 = ax.contourf(theta_SWIM , wavenumber_SWIM, spec_SWIM,\
                     np.linspace(vmin0, vmax0, 30), extend = 'both', cmap = 'viridis')
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        ax.set_rlim([0, .1])
        #cax = fig.add_axes([.65, .9, .2, 0.02])
        # cbar = plt.colorbar(p1, cax = cax_cbar, orientation = 'horizontal', ticks = [vmin0, (vmax0+vmin0)/2 ,vmax0])
        # cbar.ax.set_title('E(k, $\\theta$) [$m^{2}$/rad]')

        angles_degrees = [1, 45, 90, 135, 179]

        angles_radians = np.deg2rad(angles_degrees)
        angles_radians = np.deg2rad(angles_degrees)

        # Define radius
        radius = 0.12
        for angle in angles_radians:
            ax.plot([angle, angle], [0, radius], color='w', linestyle = '--', linewidth = 1)

        #############,
        #---White circles
        #############
        radius_080m=2*np.pi/(80)
        radius_100m=2*np.pi/(100)
        radius_200m=2*np.pi/(200)
        radius_500m=2*np.pi/(500)
        #radius_300m=(2*np.pi)/300.


        #circle1 = plt.Circle((0.0, 0.0), radius_080m, color='w',transform=ax.transData._b,linestyle='--',fill=False)
        #ax.add_patch(circle1)
        circle2 = plt.Circle((0.0, 0.0), radius_100m, color='w',transform=ax.transData._b,linestyle='--',fill=False)
        ax.add_patch(circle2)
        circle3 = plt.Circle((0.0, 0.0), radius_200m, color='w',transform=ax.transData._b,linestyle='--',fill=False)
        ax.add_patch(circle3)
        circle4 = plt.Circle((0.0, 0.0), radius_500m, color='w',transform=ax.transData._b,linestyle='--',fill=False)
        ax.add_patch(circle4)
        #ax.text(140* np.pi/180, radius_080m-0.004,'80 m',color='w', fontsize = 9)
        ax.text(140 * np.pi/180, radius_100m-.004,'100 m',color='w', fontsize = 9)
        ax.text(140 * np.pi/180, radius_200m-.004,'200 m',color='w', fontsize = 9)
        #ax.text(140 * np.pi/180, radius_500m-.004,'500 m',color='w', fontsize = 9)
        angle_labels = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W','NW']
        ax.set_thetagrids(angles = range(0, 360, 45),\
                          labels = angle_labels,fontsize=12)
        
        ax.set_xlim([0, np.pi])
        props = dict(boxstyle='round', facecolor='w', alpha=1)

        #ax.text(1.25, .13, f"Lp = {Lp:.0f} m \nDp = {(Dp):.0f} deg\nDirM =  {the_mp:.0f} deg\nDirSpr = {sigp:.0f} deg\nQkk = {Qkk:.0f} m",\
        #        bbox=props)
        #ax.text(330 * np.pi/180, .2, f"Lp = {Lp:.0f} m \nDp = {(Dp):.0f} deg\nDirM =  {the_mp:.0f} deg\nDirSpr = {sigp:.0f} deg\nQkk = {Qkk:.0f} m",\
        #        bbox=props) # for typhoon

        
    if side_spec=='left':

        p1 = ax.contourf(theta_SWIM + np.pi , wavenumber_SWIM, spec_SWIM,\
                         np.linspace(vmin0, vmax0, 30), extend = 'both', cmap = 'viridis')
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        ax.set_ylim([0, .1])
        #cax = fig.add_axes([.68, .9, .2, 0.02])
        # cbar = plt.colorbar(p1, cax = cax_cbar, orientation = 'horizontal', ticks = [vmin0, (vmax0+vmin0)/2 ,vmax0])
        # cbar.ax.set_title('E(k, $\\theta$) [$m^{2}$/rad]')
        angles_degrees = [181, 225, 270, 305, 0]

        angles_radians = np.deg2rad(angles_degrees)
        angles_radians = np.deg2rad(angles_degrees)

        # Define radius
        radius = 0.12
        for angle in angles_radians:
            ax.plot([angle, angle], [0, radius], color='w', linestyle = '--', linewidth = 1)

        #############,
        #---White circles
        #############
        radius_080m=2*np.pi/(80)
        radius_100m=2*np.pi/(100)
        radius_200m=2*np.pi/(200)
        radius_500m=2*np.pi/(500)
        #radius_300m=(2*np.pi)/300.


        #circle1 = plt.Circle((0.0, 0.0), radius_080m, color='w',transform=ax.transData._b,linestyle='--',fill=False)
        #ax.add_patch(circle1)
        circle2 = plt.Circle((0.0, 0.0), radius_100m, color='w',transform=ax.transData._b,linestyle='--',fill=False)
        ax.add_patch(circle2)
        circle3 = plt.Circle((0.0, 0.0), radius_200m, color='w',transform=ax.transData._b,linestyle='--',fill=False)
        ax.add_patch(circle3)
        circle4 = plt.Circle((0.0, 0.0), radius_500m, color='w',transform=ax.transData._b,linestyle='--',fill=False)
        ax.add_patch(circle4)
        #ax.text(20* np.pi/180, radius_080m-0.004,'80 m',color='k', fontsize = 9)
        #ax.text(20 * np.pi/180, radius_100m-.004,'100 m',color='k', fontsize = 9)
        #ax.text(20 * np.pi/180, radius_200m-.004,'200 m',color='k', fontsize = 9)
        #ax.text(20 * np.pi/180, radius_500m-.004,'500 m',color='k', fontsize = 9)
        #angle_labels = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W','NW']
        #ax.text(230* np.pi/180, radius_080m-0.004,'80 m',color='k', fontsize = 9)
        ax.text(230 * np.pi/180, radius_100m-.004,'100 m',color='w', fontsize = 9)
        ax.text(230 * np.pi/180, radius_200m-.004,'200 m',color='w', fontsize = 9)
        #ax.text(230 * np.pi/180, radius_500m-.004,'500 m',color='k', fontsize = 9)
        angle_labels = ['N', 'NW', 'W', 'SW', 'S']
        ax.set_thetagrids(angles =np.arange(359, 181, -44),
                          labels = angle_labels, fontsize=12)
        ax.set_xlim([np.pi, 2*np.pi])


        ax.grid(False)
        ax.set_frame_on(False)
        ax.set_yticks([])
        props = dict(boxstyle='round', facecolor='w', alpha=1)

    # Decomment it to display the Bulk param.
    
        #ax.text(280 * np.pi/180, .32, f"Lp = {Lp:.0f} m \nDp = {(Dp + 180):.0f} deg\nDirM =  {the_mp + 180:.0f} deg\nDirSpr = {sigp:.0f} deg\nQkk = {Qkk:.0f} m",\
        #       bbox=props)

        #ax.text(220 * np.pi/180, .28, f"Lp = {Lp:.0f} m \nDp = {(Dp + 180):.0f} deg\nDirM =  {the_mp + 180:.0f} deg\nDirSpr = {sigp:.0f} deg\nQkk = {Qkk:.0f} m",\
        #       bbox=props) # for ACC
        
        #ax.text(330 * np.pi/180, .28, f"Lp = {Lp:.0f} m \nDp = {(Dp + 180):.0f} deg\nDirM =  {the_mp + 180:.0f} deg\nDirSpr = {sigp:.0f} deg\nQkk = {Qkk:.0f} m",\
        #        bbox=props) # for typhoon
        #ax.text(300 * np.pi/180, .26, f"Lp = {Lp:.0f} m \nDp = {(Dp + 180):.0f} deg\nDirM =  {the_mp + 180:.0f} deg\nDirSpr = {sigp:.0f} deg\nQkk = {Qkk:.0f} m",\
        #        bbox=props) # for group
   #cbar = plt.colorbar(p1, cax = cax , orientation =  'horizontal', ticks = [np.amin(spec_SWIM), (np.amin(spec_SWIM) + np.amax(spec_SWIM))/2, np.amax(spec_SWIM)])
   #cbar.ax.xaxis.set_major_formatter(FuncFormatter(format_ticks))
   #cbar.ax.set_title('$\mathcal{S}$(k, $\\theta$)')
    ax.grid(False)
    ax.set_frame_on(False)
    ax.set_yticks([])
    
    plt.savefig(path_out + '/' + file_out, dpi = 300, bbox_inches = 'tight', transparent = False)
    
    
def format_ticks(x, pos):
    return f'{x:.2f}'



def remove_linear_trend(data):
    
    """
    Purpose: the SSH obs from swot has a trend within its 50 km swath. This function remove it.
    ---------
    
    Input: The 2d SSH anomaly (SSH_karin - SSH_karin.mean())
    ---------
    Output: The detrended data
    ---------
    """
    # Get the shape of the data
    rows, cols = data.shape
    
    # Initialize an array to store the detrended data
    detrended_data = np.zeros_like(data)
    
    # Loop through each column
    for i in range(cols):
        # Fit a linear model (y = mx + c) to each column
        x = np.arange(rows)  # Independent variable (row index)
        y = data[:, i]        # Dependent variable (column data)
        
        # Compute the slope (m) and intercept (c) using linear regression
        m, c = np.polyfit(x, y, deg=1)
        
        # Generate the linear trend for the column
        trend = m * x + c
        
        # Subtract the linear trend from the column data
        detrended_data[:, i] = y - trend
    
    return detrended_data


    
