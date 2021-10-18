# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 19:56:15 2020

@author: 陳煒淮
"""
from __future__ import print_function
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import math
from scipy.integrate import cumtrapz
from astropy.io import fits
from astropy.table import Table, Column, MaskedColumn
from astropy.io import ascii
from fit_kinematic_pa import fit_kinematic_pa
from cap_plot_velfield import plot_velfield
from cap_symmetrize_velfield import symmetrize_velfield
import sys
from datahandling import filehandling,match,comparePA,newPAdiff,rewrite
from fit_kpa_kstest import fit_kpa_ks
from scipy import stats
from scipy import interpolate,stats
from scipy.spatial import distance
from scipy.interpolate import interp1d, interp2d
from pprint import pprint
from plotbin.display_bins_generators import display_bins_generators
from plotbin.display_pixels import display_pixels
from plotbin.display_bins import display_bins
from vorbin.voronoi_2d_binning import voronoi_2d_binning,_sn_func
from datahandling import findPAdiff
from measurejPA_class import jPA
from kinemetry.kinemetry import kinemetry
from kinemetry.run_kinemetry_examples import plot_kinemetry_profiles_velocity,plot_kinemetry_maps

def kPA(plateifu, data, source,re_criterion_list=[1,0.5,0.3],plot=True,snthreshold=3,dap='SPX',measure='Capellari',para='MaNGA',binning=True,mangadir=None):
    assert source in ['STAR', 'HA', 'O3']
    assert dap in [None,'SPX','HYB10','VOR10'] # GAU & MASTARHC2 or SPX or HYB10 or VOR10
    assert measure in ['Capellari', 'YM', 'KS', 'Kinemetry','None'] # to get kPA
    assert para in ['NSA', 'MaNGA'] # Re, delta, BA source
    assert binning in [True, False] # apply Voronoi binning

    if measure=='KS':col=5+1
    if measure!='KS':col=3+1
    if plot:plt.figure(figsize=(3.5*col,10))

    arcsecpp = {
        'NVSS':15.,
        'FIRST':1.8,
        'FIRST CATALOG':1.8,
        'GIVENPA':1.8,
        'VLASS':1.
    }

    if measure!='None':my_decision,jpa = data[(plateifu,'my_decision','str')], data[(plateifu,'jPA','float')]
    if plot:
        # Open radio image and plot
        if my_decision in ['first catalog']:radio_image=fits.open(f'/Volumes/SDrive/yenting_pa_alignment/radio_and_optical_images/first/ifu={plateifu}_ra={ra:.3f}_dec={dec:.3f}_first.fits')['Primary'].data
        else:radio_image=fits.open(f'/Volumes/SDrive/yenting_pa_alignment/radio_and_optical_images/{my_decision}/ifu={plateifu}_ra={ra:.3f}_dec={dec:.3f}_{my_decision}.fits')['Primary'].data
        size = (radio_image.shape[0]-1)/2*arcsecpp[my_decision]
        side = np.arange(-size, size+arcsecpp[my_decision], arcsecpp[my_decision])
        assert abs(side[0])-abs(side[-1])<arcsecpp[my_decision]/2
        assert len(side)==radio_image.shape[0]
        xnew, ynew = map(np.ravel, np.meshgrid(side, side))  
        plt.subplot(3,col,1)
        radio_image-=np.min(radio_image)-1e-5
        plt.tricontour(-xnew, ynew,-2.5*np.log10(radio_image.ravel()/np.max(radio_image.ravel())), levels=np.arange(0,5,0.5),linewidths=1)
        plt.plot(-side,side*np.tan(np.radians(jpa+90)), 'k--', linewidth=1) # PA
        
        plt.subplot(3,col,col)
        optical = mpimg.imread('/Volumes/SDrive/yenting_pa_alignment/radio_and_optical_images/optical2/name=%s_optical.jpg'%(plateifu))
        plt.imshow(optical)
        plt.title('SDSS optical image')
        plt.axis('off')
        # if my_decision=='first catalog':
        #     jPA(plateifu, data=data, source='first', perc = data[(plateifu,'perc','float')], deg = 45 ,crop=0,redshift=data[(plateifu,'z_1','float')],givenpa=jpa,top=data[(plateifu,'top','int')])
        # else:jPA(plateifu, data=data, source=my_decision, perc = data[(plateifu,'perc','float')], deg = 45 ,crop=0,redshift=data[(plateifu,'z_1','float')],givenpa=None,top=data[(plateifu,'top','int')])
        # display_pixels(xnew, ynew, radio_image.ravel(),pixelsize=np.min(distance.pdist(np.column_stack([xnew, ynew]))),cmap='plasma')
    
    # lower bound for number of pixels
    crit_pixel = {
        '12':50,
        '91':50,
        '61':50,
        '37':30,
        '19':20
    }
    #Ha-6564
    halpha_channel={
        'MASTARHC2':23,
        'GAU':18
    }
    #Oiii-5008
    oiii_channel={
        'MASTARHC2':16,
        'GAU':13
    }
    
    if dap=='SPX' or dap=='HYB10' or dap=='VOR10':
        ha_ch = halpha_channel['MASTARHC2']
        oiii_ch = oiii_channel['MASTARHC2']
    elif dap==None:
        ha_ch = halpha_channel[data[(plateifu,'manga_data')]]
        oiii_ch = oiii_channel[data[(plateifu,'manga_data')]]

    # print('starting',plateifu)
    output=[plateifu]

    # Open MaNGA Maps file. Extract x,y, flux, and STAR velocity.
    if mangadir==None:galaxy=fits.open(f'/Volumes/SDrive/yenting_pa_alignment/MaNGA/{dap}-MILESHC-MASTARSSP/manga-{plateifu}-MAPS-{dap}-MILESHC-MASTARSSP.fits.gz')
    else:galaxy=fits.open(f'{mangadir}/manga-{plateifu}-MAPS-{dap}-MILESHC-MASTARSSP.fits.gz')

    if source=='STAR':
        F_map = galaxy['SPX_MFLUX'].data
        F_ivar = galaxy['SPX_MFLUX_IVAR'].data
        F_mask = galaxy['SPX_SNR'].data
        V_map = galaxy['STELLAR_VEL'].data
        V_ivar=galaxy['STELLAR_VEL_IVAR'].data
        mask = galaxy['STELLAR_VEL_MASK'].data
    if source=='HA':
        F_map = galaxy['EMLINE_GFLUX'].data[ha_ch]
        F_ivar = galaxy['EMLINE_GFLUX_IVAR'].data[ha_ch]
        F_mask = galaxy['EMLINE_GFLUX_MASK'].data[ha_ch]
        V_map = galaxy['EMLINE_GVEL'].data[ha_ch]
        V_ivar = galaxy['EMLINE_GVEL_IVAR'].data[ha_ch]
        mask = galaxy['EMLINE_GVEL_MASK'].data[ha_ch]
    if source=='O3':
        F_map = galaxy['EMLINE_GFLUX'].data[oiii_ch]
        F_ivar = galaxy['EMLINE_GFLUX_IVAR'].data[oiii_ch]
        F_mask = galaxy['EMLINE_GFLUX_MASK'].data[oiii_ch]
        V_map = galaxy['EMLINE_GVEL'].data[oiii_ch]
        V_ivar = galaxy['EMLINE_GVEL_IVAR'].data[oiii_ch]
        mask = galaxy['EMLINE_GVEL_MASK'].data[oiii_ch]

    # rather than direct indicing as done in YM's, a true middle is given by MaNGA
    X = galaxy['SPX_SKYCOO'].data[0] #left=east=+RA
    Y = galaxy['SPX_SKYCOO'].data[1]
    if para=='MaNGA':
        r = galaxy['SPX_ELLCOO'].data[1]
    p50 = data[(plateifu,'NSA_ELPETRO_TH50_R','float')]
    if para=='NSA':
        delta = np.radians(data[(plateifu,'NSA_ELPETRO_PHI','float')])
        ba = data[(plateifu,'NSA_ELPETRO_BA','float')]
        if np.isnan(p50):
            delta,ba,p50=0,1,V_map.shape[0]/4
            print(f'{plateifu} using default rather than nsa values')
        xpos = X*np.cos(delta)-Y*np.sin(delta)
        ypos = X*np.sin(delta)+Y*np.cos(delta)
        r = np.sqrt((xpos/ba)**2+ypos**2)/p50
    v0 = interp2d(X[int(V_map.shape[0]/2)],Y.T[int(V_map.shape[1]/2)],V_map)(0,0)
    Ferr = 1/np.sqrt(F_ivar)
    snr = F_map*np.sqrt(F_ivar)
    
    V_map -= v0 
    # https://www.sdss.org/dr15/algorithms/bitmasks/#MANGA_DAPPIXMASK
    mask_any = np.heaviside(mask,0)
    mask = np.array([[int(bin(et).replace('0b','')) for et in row] for row in mask])
    mask=[(mask//10**i%2).astype(int) for i in range(31)]

    applybinning = False
    if binning:
        idx = ((V_ivar>0.)&(mask[30]==0.)).nonzero()
        at,bt,ct,dt,et=map(np.ravel,[X[idx],Y[idx],F_map[idx],Ferr[idx],V_map[idx]])
        # at,bt,ct,dt,et=map(np.ravel,[X[idx],Y[idx],abs(V_map[idx]),1/np.sqrt(V_ivar[idx]),V_map[idx]])
        assert _sn_func(np.flatnonzero(dt > 0), ct, dt) > snthreshold,f'{_sn_func(np.flatnonzero(dt > 0), ct, dt)} should be larger than {snthreshold}'
        if np.nanmin(ct/dt)>=snthreshold:
            xBin, yBin, sn, V_mapBin = at,bt,ct/dt,et
        if np.nanmin(ct/dt)<snthreshold:
            binNum, xBin, yBin, xBar, yBar, sn, nPixels, scale = voronoi_2d_binning(at,bt,ct,dt,snthreshold,plot=False)
            xyIn = np.column_stack([at,bt])
            xyOut = np.column_stack([xBin,yBin])
            V_mapBin = interpolate.griddata(xyIn, et, xyOut, method='linear')
            applybinning=True

    # plotting
    if plot:
        plt.subplot(3,col,1)
        if not applybinning:
            idx = (mask[30]==0.).nonzero()
            at,bt,ct,et=map(np.ravel,[X[idx],Y[idx],V_map[idx],r[idx]])
            # plt.clabel(plt.tricontour(at,bt,et,levels=re_criterion_list[::-1],colors='black'))
            plt.tricontour(at,bt,et,levels=re_criterion_list[::-1],colors='black')
            display_pixels(at,bt,ct, vmin=np.nanmin(ct),vmax=np.nanmax(ct),colorbar=True)
            plt.title(f'velocity map, unbinned (m30)\njPA={jpa}$^\circ$({my_decision})\nNSA TH50_R={p50} arcsec',size=10)
        if applybinning:
            idx = ((V_ivar>0.)&(mask[30]==0.)).nonzero()
            at,bt,ct=map(np.ravel,[X[idx],Y[idx],r[idx]])
            # plt.clabel(plt.tricontour(at,bt,ct,levels=re_criterion_list[::-1],colors='black'))
            plt.tricontour(at,bt,ct,levels=re_criterion_list[::-1],colors='black')
            display_bins_generators(xBin, yBin, V_mapBin, at, bt, vmin=np.nanmin(V_mapBin),vmax=np.nanmax(V_mapBin),colorbar=True)
            plt.title(f'velocity map, binned (m30,+e)\njPA={jpa}$^\circ$({my_decision})\nNSA TH50_R={p50} arcsec',size=10)
        plt.xlim([-2*np.max(abs(at)),2*np.max(abs(at))])
        plt.ylim([-2*np.max(abs(bt)),2*np.max(abs(bt))])
        plt.gca().invert_xaxis()

        plt.subplot(3,col,col+1)
        if not applybinning:
            idx = ((V_ivar>0.)&(mask[30]==0.)).nonzero()
            at,bt,ct=map(np.ravel,[X[idx],Y[idx],snr[idx]])
            plt.tricontour(at,bt,ct,levels=[snthreshold],colors='red')
            display_pixels(at,bt,ct, vmin=np.nanmin(ct),vmax=np.nanmax(ct),colorbar=True,cmap='viridis')
            if source=='STAR':plt.title(f'mean g band SNR, unbinned\n(m30,+e)\nred line has SNR of {snthreshold}',size=10)
            else:plt.title(f'{source} SNR, unbinned\n(m30,+e)\nred line has SNR of {snthreshold}',size=10)
        if applybinning:
            idx = ((V_ivar>0.)&(mask[30]==0.)).nonzero()
            at,bt,ct=map(np.ravel,[X[idx],Y[idx],snr[idx]])
            plt.tricontour(xBin, yBin, sn,levels=[snthreshold],colors='red')
            display_bins_generators(xBin, yBin, sn, at, bt, vmin=np.nanmin(sn),vmax=np.nanmax(sn),colorbar=True)
            if source=='STAR':plt.title(f'mean g band SNR, binned\n(m30,+e)\nred line has SNR of {snthreshold}',size=10)
            else:plt.title(f'{source} SNR, binned\n(m30,+e)\nred line has SNR of {snthreshold}',size=10)
        plt.gca().invert_xaxis()

        plt.subplot(3,col,2*col+1)
        at,bt,ct=map(np.ravel,[X,Y,mask_any])
        display_pixels(at,bt,ct, vmin=0,vmax=1,colorbar=False,cmap='viridis')
        plt.gca().invert_xaxis()
        plt.title('all possible mask')


    for b,re_criterion in enumerate(re_criterion_list):
        if b>2:plot=False
        # binning or not differs in noise related
        # ks or not differs in re related
        if binning : 
            index = ((r<re_criterion)&(mask_any!=1)).nonzero()
            index_kin = ((mask_any!=1)).nonzero()
            # index_ks = ((mask_any!=1)).nonzero()
        if not binning : 
            index = ((snr>=snthreshold)&(V_ivar>0.)&(r<re_criterion)&(mask_any==0)).nonzero()
            index_kin = ((snr>=snthreshold)&(V_ivar>0.)&(mask_any==0)).nonzero()
            # index_ks = ((snr>snthreshold)&(V_ivar>0.)&(mask_any!=1)).nonzero()
                
        countgoodpix = len(index[0])
        applybinning = False
        if binning:
            at,bt,ct,dt,et=map(np.ravel,[X[index],Y[index],F_map[index],Ferr[index],V_map[index]])
            # at,bt,ct,dt,et=map(np.ravel,[X[index],Y[index],abs(V_map[index]),1/np.sqrt(V_ivar[index]),V_map[index]])
            assert _sn_func(np.flatnonzero(dt > 0), ct, dt) > snthreshold,f'{_sn_func(np.flatnonzero(dt > 0), ct, dt)} should be larger than {snthreshold}'
            if np.min(ct/dt)>snthreshold:
                xBin, yBin, sn, V_mapBin = at,bt,ct/dt,et
            if np.min(ct/dt)<snthreshold:
                binNum, xBin, yBin, xBar, yBar, sn, nPixels, scale = voronoi_2d_binning(at,bt,ct,dt,snthreshold,plot=False)
                idx = np.column_stack([xBin,yBin])
                xyIn = np.column_stack([at,bt])
                xyOut = np.column_stack([xBin,yBin])
                V_mapBin = interpolate.griddata(xyIn, et, xyOut, method='linear')
                applybinning=True
        
        # get kinematic PA
        if measure=='KS':
            if plot:
                ax2=plt.subplot(3,col,col*b+col-1)
                ax1=plt.subplot(3,col,col*b+col-2)
            if not plot:ax1,ax2=None,None
            if not binning:[PA,dstat,pval,pxl,bins,counthalf1,counthalf2]= fit_kpa_ks(X[index],Y[index],V_map[index],float(re_criterion),ax1=ax1,ax2=ax2,binning=binning,debug=False)
            if binning:[PA,dstat,pval,pxl,bins,counthalf1,counthalf2]= fit_kpa_ks(xBin,yBin,V_mapBin,float(re_criterion),ax1=ax1,ax2=ax2,binning=binning)
            output.extend([PA,dstat,pval,pxl,countgoodpix])
        if measure=='Kinemetry':
            pix = np.nanmin([re_criterion*p50/0.5,data[(plateifu,'nomial FOV radius (arcsec)','float')]*np.sqrt(3)])
            if not binning:k = kinemetry(xbin=X[index_kin], ybin=Y[index_kin], moment=V_map[index_kin], error=np.sqrt(1/V_ivar[index_kin]), scale=0.5, radius=np.full(3,pix), cover=0.0, plot=True)
            if binning:k = kinemetry(xbin=xBin, ybin=yBin, moment=V_mapBin, scale=0.5, radius=np.full(3,pix), cover=0.0, plot=False)
            PA,angErr,vSyst = -k.pa[0]%180,k.er_pa[0],k.vsys
            k0 = k.cf[:,0]
            k1 = np.sqrt(k.cf[:,1]**2 + k.cf[:,2]**2)
            k5 = np.sqrt(k.cf[:,5]**2 + k.cf[:,6]**2)
            k51 = k5/k1
            erk1 = (np.sqrt( (k.cf[:,1]*k.er_cf[:,1])**2 + (k.cf[:,2]*k.er_cf[:,2])**2 ))/k1
            erk5 = (np.sqrt( (k.cf[:,5]*k.er_cf[:,5])**2 + (k.cf[:,6]*k.er_cf[:,6])**2 ))/k5
            erk51 = ( np.sqrt( ((k5/k1) * erk1)**2 + erk5**2  ) )/k1 
            # plot_kinemetry_profiles_velocity(k)
            # plot_kinemetry_maps(X[index], Y[index], V_map[index], k)
            output.extend([k.rad[0],PA,angErr,k51[0],erk51[0]])
        if measure=='Capellari':
            if not binning:PA,angErr,vSyst = fit_kinematic_pa(X[index],Y[index],V_map[index], debug=False, nsteps=361,quiet=False, plot=False)
            if binning:PA,angErr,vSyst = fit_kinematic_pa(xBin,yBin,V_mapBin, debug=False, nsteps=361,quiet=False, plot=False)
            PA=-PA%180.        
            output.extend([PA,angErr,vSyst])
        if measure=='YM':
            print('PA_'+str(source)+'_'+str(re_criterion)+'RE')
            PA = float(data[(str(plateifu),'PA_'+str(source)+'_'+str(re_criterion)+'RE')])
            angErr = float(data[(str(plateifu),'PA_'+str(source)+'_ERR_'+str(re_criterion)+'RE')])
            vSyst = float(data[(str(plateifu),'PA_'+str(source)+'_VELOFF_'+str(re_criterion)+'RE')])
            output.extend([PA,angErr,vSyst])
        if measure in ['Kinemetry','Capellari','YM']:
            half1_list = [e for i,e in enumerate(V_map[index]) if X[index][i]*np.sin(np.radians(PA))+Y[index][i]*np.cos(np.radians(PA))>0]
            half2_list = [e for i,e in enumerate(V_map[index]) if X[index][i]*np.sin(np.radians(PA))+Y[index][i]*np.cos(np.radians(PA))<0]
            if len(half1_list)!=0 and len(half2_list)!=0 : [dstat,pval] = stats.ks_2samp(half1_list, half2_list) # The null hypothesis is that the two distributions are identical
            else : dstat,pval=np.nan,np.nan
            bins = [i for i in np.linspace(np.nanmin(V_map[index]),np.nanmax(V_map[index]),50)]
            counthalf1,_ = np.histogram(half1_list, bins=bins)
            counthalf2,_ = np.histogram(half2_list, bins=bins)
            output.extend([dstat,pval,countgoodpix])
        if measure!='None':PAdiff = findPAdiff([jpa], [PA])[0]
        # output.append(PAdiff)
        
        # Plot results
        if plot:
            plt.subplot(3,col,col*b+2)
            if not applybinning:
                size = max(max(abs(X[index])),max(abs(Y[index])))
                side = np.array([-size, size])
                plt.plot(side,side/np.tan(np.radians(PA)), 'limegreen', linewidth=3) # PA
                plt.plot(side,-side*np.tan(np.radians(PA)), 'k--', linewidth=3) # Zero-velocity line
                if index[0].shape[0]>2:display_pixels(X[index], Y[index], V_map[index], vmin=np.nanmin(V_map[index]),vmax=np.nanmax(V_map[index]),colorbar=False)
                if measure=='KS':
                    if dap=='SPX':plt.title(f'kPA={PA:.1f}$^\circ$,PAdiff={PAdiff:.1f}$^\circ$ unbinned\n(Re<{re_criterion},SNR>={snthreshold},ma,+e)\nhas {countgoodpix} good pixels',size=10)
                    if dap=='VOR10':plt.title(f'kPA={PA:.1f}$^\circ$(Re<{re_criterion},ma,+e)\nhas {countgoodpix} good pixels',size=10)
                if measure!='KS':plt.title(f'kPA={PA:.1f}+/-{angErr:.2f}$^\circ$,PAdiff={PAdiff:.1f}$^\circ$, V_offset={vSyst}, unbinned\n(Re<{re_criterion},SNR>={snthreshold},ma,+e)\nhas {countgoodpix} good pixels',size=10)
            if applybinning:
                size = max(max(abs(xBin)),max(abs(yBin)))
                side = np.array([-size, size])
                plt.plot(side,side/np.tan(np.radians(PA)), 'limegreen', linewidth=3) # PA
                plt.plot(side,-side*np.tan(np.radians(PA)), 'k--', linewidth=3) # Zero-velocity line
                display_bins_generators(xBin,yBin,V_mapBin, X[index], Y[index], vmin=np.nanmin(V_mapBin),vmax=np.nanmax(V_mapBin),colorbar=True)
                if measure=='KS':
                    if dap=='SPX':plt.title(f'kPA={PA:.1f}$^\circ$,PAdiff={PAdiff:.1f}$^\circ$ binned\n(Re<{re_criterion},SNR>={snthreshold},ma,+e)\nhas {countgoodpix} good pixels',size=10)
                    if dap=='VOR10':plt.title(f'kPA={PA:.1f}$^\circ$(Re<{re_criterion},ma,+e)\nhas {countgoodpix} good pixels',size=10)
                if measure!='KS':plt.title(f'kPA={PA:.1f}+/-{angErr:.2f}$^\circ$,PAdiff={PAdiff:.1f}$^\circ$, V_offset={vSyst:.2f}, binned\n(Re<{re_criterion},SNR>={snthreshold},ma,+e)\nhas {countgoodpix} good pixels',size=10)
            plt.gca().invert_xaxis()

            plt.subplot(3,col,col*b+3)
            if index[0].shape[0]>2:plt.plot([(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)],counthalf1/counthalf1.sum())
            if index[0].shape[0]>2:plt.plot([(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)],counthalf2/counthalf2.sum())
            # plt.title(f'D statistic:{dstat:.4f},p-value:{pval:.2E}')
            plt.title(f'D statistic:{dstat:.4f}')
            plt.xlabel('velocity (km/s)')

        if plot and b==2:
            maxarcsec_rad = data[(plateifu,'nomial FOV radius (arcsec)','float')]*np.sqrt(3)#/2*2
            rad_list=[e+0.1 for e in np.arange(0.5,30,0.5) if e<maxarcsec_rad]

            rad_ks,pa_ks,paerr_ks = [],[],[]
            rad_kin,pa_kin,paerr_kin,k51_kin,k51err_kin = [],[],[],[],[]
            rad_capellari,pa_capellari,paerr_capellari = [],[],[]
            for arcsec_criterion in rad_list:
                if not binning : 
                    index = ((snr>=snthreshold)&(V_ivar>0.)&(r<arcsec_criterion/p50)&(mask_any==0)).nonzero()
                    index_kin = ((snr>=snthreshold)&(V_ivar>0.)&(mask_any==0)).nonzero()
                if binning : 
                    index = ((r<arcsec_criterion/p50)&(mask_any!=1)).nonzero()
                    index_kin = ((mask_any!=1)).nonzero()

                # try:
                if not binning : 
                    PA=fit_kpa_ks(X[index],Y[index],V_map[index],arcsec_criterion/p50,binning=binning)[0]
                if binning : 
                    PA=fit_kpa_ks(xBin,yBin,V_mapBin,arcsec_criterion/p50,binning=binning)[0]
                rad_ks.append(arcsec_criterion)
                pa_ks.append(PA)
                # paerr_ks.append(PAerr)
                # except Exception:pass

                # try:
                pix = arcsec_criterion/0.5
                if not binning:k = kinemetry(xbin=X[index_kin], ybin=Y[index_kin], moment=V_map[index_kin], error=np.sqrt(1/V_ivar[index_kin]), scale=0.5, radius=np.full(3,pix), cover=0.0, plot=False)
                if binning:k = kinemetry(xbin=xBin, ybin=yBin, moment=V_mapBin, scale=0.5, radius=np.full(3,pix), cover=0.0, plot=False)
                PA,angErr = -k.pa[0]%180,min(k.er_pa[0],90)
                k0 = k.cf[:,0]
                k1 = np.sqrt(k.cf[:,1]**2 + k.cf[:,2]**2)
                k5 = np.sqrt(k.cf[:,5]**2 + k.cf[:,6]**2)
                k51 = k5/k1
                erk1 = (np.sqrt( (k.cf[:,1]*k.er_cf[:,1])**2 + (k.cf[:,2]*k.er_cf[:,2])**2 ))/k1
                erk5 = (np.sqrt( (k.cf[:,5]*k.er_cf[:,5])**2 + (k.cf[:,6]*k.er_cf[:,6])**2 ))/k5
                erk51 = ( np.sqrt( ((k5/k1) * erk1)**2 + erk5**2  ) )/k1 
                rad_kin.append(k.rad[0])  
                pa_kin.append(PA)
                paerr_kin.append(angErr)
                k51_kin.append(k51[0])
                k51err_kin.append(erk51[0])
                # except Exception:pass

                try:
                    if not binning : 
                        PA,angErr,_=fit_kinematic_pa(X[index],Y[index],V_map[index], debug=False, nsteps=361,quiet=False, plot=False)
                    if binning : 
                        PA,angErr,_=fit_kinematic_pa(xBin,yBin,V_mapBin, debug=False, nsteps=361,quiet=False, plot=False)
                    PA=-PA%180.      
                    rad_capellari.append(arcsec_criterion)  
                    paerr_capellari.append(angErr)
                    pa_capellari.append(PA)
                except Exception:pass

            plt.subplot(3,col,2*col)
            plt.errorbar(rad_kin,pa_kin,3*np.array(paerr_kin),capsize=3,fmt='.g',label='Kinemetry')
            plt.plot(rad_ks,pa_ks,color='k',label='KS')
            plt.errorbar(rad_capellari,pa_capellari,paerr_capellari,fmt='.b',label='Capellari',capsize=3,linewidth=1)
            plt.legend()
            for criterion in re_criterion_list: 
                plt.axvline(criterion*p50,color='k',linestyle='--')
            plt.title('kPA')
            plt.xlabel('arcsec')
            plt.xticks(rad_kin[::3],fontsize=8, rotation=0)
            plt.ylabel('degrees')
            plt.yticks(np.linspace(0, 180, 7))
            # plt.ylim(0,200)

            plt.subplot(3,col,3*col)
            plt.errorbar(rad_kin,k51_kin,k51err_kin,capsize=3,fmt='.g')
            for criterion in re_criterion_list: 
                plt.axvline(criterion*p50,color='k',linestyle='--')
            plt.title('k5/k1')
            plt.xlabel('arcsec')
            plt.xticks(rad_kin[::3],fontsize=8, rotation=0)
            # plt.ylim(0,1)
            plt.yscale('log')
            plt.axhline(0.04,color='blue',linestyle='--')
            
            plt.tight_layout()

    output.extend([index])
    return output

if __name__ == '__main__':
    data = filehandling("500_e6.csv")
    print(kPA('11755-3701', data,re_criterion_list=[1,0.5,0.3],plot=False,source='STAR',dap='VOR10',para='NSA',measure='KS',binning=False,snthreshold=0))
    plt.show()
    # for et in data.extract('plateifu',selection='hasjpa'):
    #     print(kPA(et, data, source='STAR',re_criterion_list=[0.3,0.5,1],snthreshold=3,dap='SPX',measure=False))
    #     plt.show()

        # if not binning:PA,angErr,vSyst = fit_kinematic_pa(X[index],Y[index],V_map[index], debug=False, nsteps=361,quiet=False, plot=False)
        # if binning:PA,angErr,vSyst = fit_kinematic_pa(xBin,yBin,V_mapBin, debug=False, nsteps=361,quiet=False, plot=False)
        # PA=-PA%180.        

# 11755-3701
# 10214-3703
# 12700-1901
# 11011-9102
# 8438-12703
# 8614-12701

# # extract data
# plt.subplot(3,col,3*col)
# plt.axis('off')
# plt.axis([0, 10, 0, 10])
# plt.text(0, 10, f'{plateifu}', fontsize=10)
# sm = data[(plateifu,'NSA_ELPETRO_MASS','float')]
# plt.text(0, 9, f'stellar mass {sm:.2E} $M_\odot/h^2$', fontsize=10)
# hm = data[(plateifu,'Mh_L','float')]
# plt.text(0, 8, f'halo mass {10**hm:.2E} $M_\odot/h$', fontsize=10)
# eps = data[(plateifu,'ym_epsilon','float')]
# plt.text(0, 7, f'YM epsilon {eps:.4f}', fontsize=10)
# lam = data[(plateifu,'ym_lambda1','float')]
# plt.text(0, 6, f'YM lamda {lam:.4f}', fontsize=10)
# if (lam<(eps*0.25+0.1) and eps<0.4):plt.text(0, 5, 'is fast rotator', fontsize=10)
# else:plt.text(0, 5, 'is not fast rotator', fontsize=10)
# objgp = data[(plateifu,'objgp(arcsec)','float')]
# r180 = data[(plateifu,'r_200','float')]
# plt.text(0, 4, f'distance to gp center {objgp/r180:.4f} r180', fontsize=10)
# nog = data[(plateifu,'300kpc 1500 km/s','float')]
# plt.text(0, 3, f'{int(nog)} galaxies in 300 kpc 1500 km/s', fontsize=10)
# near = data[(plateifu,'5th nearest (mpc) by cylindrical method cut = 1000 km/s','float')]
# plt.text(0, 2, f'5th nearest is {near:.4f} mpc away', fontsize=10)
# bcg = data[(plateifu,'BCG','float')]
# plt.text(0, 1, f'BCG {str(bcg==1)}', fontsize=10)
# mmg = data[(plateifu,'MMG','float')]
# plt.text(0, 0, f'MMG {str(mmg==1)}', fontsize=10)

# These are what YM did
# if source=='STAR':
#     if binning : 
#         index = ((r<re_criterion)&(mask[1]!=1)).nonzero()
#         # index_ks = ((mask[1]!=1)).nonzero()
#     if not binning : 
#         index = ((V_ivar>0.)&(r<re_criterion)&(mask[1]!=1)).nonzero()
#         # index_ks = ((V_ivar>0.)&(mask[1]!=1)).nonzero()
# if source=='HA' or source=='O3':
#     if binning:
#         index = ((r<re_criterion)&(mask[1]!=1)).nonzero()
#         # index_ks = ((mask[1]!=1)).nonzero()
#     if not binning:
#         index = ((F_map/Ferr>snthreshold)&(V_ivar>0.)&(r<re_criterion)&(mask[1]!=1)).nonzero()
#         # index_ks = ((F_map/Ferr>snthreshold)&(V_ivar>0.)&(mask[1]!=1)).nonzero()

# # plot png result from kinemetry
# if data[(plateifu,'kinemetry','int')]==1:
#     ax = plt.subplot(3,col,2*col)
#     pos = ax.get_position()
#     posnew = [pos.x0*0.94, pos.y0*0.78,  pos.width*1.7, pos.height*1.7 ]
#     ax.set_position(posnew)    
#     plt.imshow(mpimg.imread(f'/Volumes/SDrive/yenting_pa_alignment/results/kinemetry/{dap}/{plateifu}.png'))
#     plt.margins(0.1,0.1)
#     plt.axis('off')