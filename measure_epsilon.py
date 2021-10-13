import numpy as np
import math as math
import matplotlib.pyplot as plt
from astropy.table import Table, Column, MaskedColumn
from astropy.io import ascii
from astropy.io import fits
import sys
from LambdaR_2D_General import Derive_LR_VS_Profiles
from pynmap_master.pynmap.misc_functions import characterise_ima2d

def interval(x,list):
    found = 0
    for i in range(len(list)-1):
        if x>list[i] and x<list[i+1]:
            found = 1
            return i
    if found == 0:
        return len(list)-1

#parameters
maskk=False
what=0
numofgalaxies = 200
eps0=0.
pa0=0.
re_criterion_list=[0.3,0.5,1]


plate,epspt2,epspt5,eps1,lambdapt2profile,lambdapt5profile,lambda1profile,lambdapt2default,lambdapt5default,lambda1default=[],[],[],[],[],[],[],[],[],[]
for platenumber in np.arange(what,what+numofgalaxies,1):
    list_objects = ascii.read("D:\\yenting\\codes\\extracted_mpl11_my.txt")
    list_plate = list_objects['plateifu']
    list_mastarhc2 = list_objects['GAU']
    halpha_channel=[23,18] #Ha-6564
    #halpha_channel=[16,13] #Oiii-5008
    list_ra = list_objects['objra']
    list_dec = list_objects['objdec']
    list_z = list_objects['z_1']
    list_theta =list_objects['ELPETRO_TH50_R']
    list_Re = []
    for b in range(len(re_criterion_list)):
        list_Re.append(re_criterion_list[b]*list_theta[platenumber]) #in arcsec
    #print(list_ra[platenumber],list_dec[platenumber],list_plate[platenumber])

    #==================================MaNGA=======================================#
    # open MaNGA Maps file and extract x,y, flux, and stellar velocity
    if list_mastarhc2[platenumber]==0: galaxy=fits.open('D:\\yenting\\MaNGA\\all_map_DAP\\manga-%s-MAPS-HYB10-MILESHC-MASTARHC2.fits' %(list_plate[platenumber]))
    else: galaxy=fits.open('D:\\yenting\\MaNGA\\all_map_DAP\\manga-%s-MAPS-HYB10-GAU-MILESHC.fits' %(list_plate[platenumber]))
    X = galaxy['SPX_SKYCOO'].data[0]
    Y = galaxy['SPX_SKYCOO'].data[1]
    F = galaxy['SPX_MFLUX'].data
    V = galaxy['STELLAR_VEL'].data
    # V = galaxy['EMLINE_GVEL'].data[halpha_channel[list_mastarhc2[platenumber]]]
    S = galaxy['STELLAR_SIGMA'].data
    # S = galaxy['EMLINE_GSIGMA'].data[halpha_channel[list_mastarhc2[platenumber]]]
    mV = galaxy['STELLAR_VEL_MASK'].data
    # mV = galaxy['EMLINE_GVEL_MASK'].data[halpha_channel[list_mastarhc2[platenumber]]]
    mS = galaxy['STELLAR_SIGMA_MASK'].data
    # mS = galaxy['EMLINE_GSIGMA_MASK'].data[halpha_channel[list_mastarhc2[platenumber]]]
    mask = mS > 1e9

    #eps and pa profile
    if maskk==True:r, lumprof, eps, pa, re_50, eps_50, pa_50 = characterise_ima2d(X, Y, F, mask = mask)
    else:r, lumprof, eps, pa, re_50, eps_50, pa_50 = characterise_ima2d(X, Y, F)
    #print(eps, pa)
    r_index = [interval(list_Re[b],r) for b in range(len(list_Re))]
    
    # 1 go where you use the PA and Eps profiles but beware: at small radii they are noisy
    if maskk==True:result = Derive_LR_VS_Profiles(X[~mask], Y[~mask], F[~mask], V[~mask],S[~mask], R_EpsPA=r, Eps=eps, PA=pa)
    else:result = Derive_LR_VS_Profiles(X, Y, F, V,S, R_EpsPA=r, Eps=eps, PA=pa)
    #print(result.LambdaR)
    effrad_index = [interval(list_Re[b],result.EffRadius) for b in range(len(list_Re))]
    
    # 2 go where you use a single value for PA and Eps
    if maskk==True:result2 = Derive_LR_VS_Profiles(X[~mask], Y[~mask], F[~mask], V[~mask], S[~mask], Eps=eps0, PA=pa0)
    else:result2 = Derive_LR_VS_Profiles(X, Y, F, V, S, Eps=eps0, PA=pa0)
    #print(result2.LambdaR)
    effrad2_index = [interval(list_Re[b],result2.EffRadius) for b in range(len(list_Re))]

    #print(r_index,effrad_index,effrad2_index)
    #plot
    fig,axs=plt.subplots(3,2,constrained_layout=True,figsize=(10, 10))
    fig.suptitle('left:set pa=%g and eps=%g ; right: use the PA and Eps profiles ; %gRe = %.2f,%gRe = %.2f,%gRe = %.2f (arcsec)' %(pa0,eps0,re_criterion_list[0],list_Re[0],re_criterion_list[1],list_Re[1],re_criterion_list[2],list_Re[2]))
    axs[0][0].plot(r, pa)
    axs[0][0].set_title('pa = %.2f(%gRe),%.2f(%gRe),%.2f(%gRe)\nhas mean = %.2f'%(pa[r_index[0]] if r_index[0]<len(pa) else 9999,re_criterion_list[0],pa[r_index[1]] if r_index[1]<len(pa) else 9999,re_criterion_list[1],pa[r_index[2]] if r_index[2]<len(pa) else 9999 ,re_criterion_list[2],np.mean(pa)))
    axs[0][0].set_xlabel('radius(arcsec)')
    axs[0][0].set_ylabel('pa')
    axs[0][1].plot(r,eps)
    axs[0][1].set_title('eps = %.2f(%gRe),%.2f(%gRe),%.2f(%gRe)\nhas mean = %f'%(eps[r_index[0]] if r_index[0]<len(eps) else 9999,re_criterion_list[0],eps[r_index[1]] if r_index[1]<len(eps) else 9999,re_criterion_list[1],eps[r_index[2]] if r_index[2]<len(eps) else 9999 ,re_criterion_list[2],np.mean(eps)))
    axs[0][1].set_xlabel('radius(arcsec)')
    axs[0][1].set_ylabel('eps')

    # You can see the two results
    axs[1][0].plot(result2.EffRadius, result2.VS)
    axs[1][0].set_title('V/Sig = %.2f(%gRe),%.2f(%gRe),%.2f(%gRe)\nhas mean = %f'%(result2.VS[effrad2_index[0]] if effrad2_index[0]<len(result2.VS) else 9999,re_criterion_list[0],result2.VS[effrad2_index[1]] if effrad2_index[1]<len(result2.VS) else 9999,re_criterion_list[1],result2.VS[effrad2_index[2]] if effrad2_index[2]<len(result2.VS) else 9999 ,re_criterion_list[2],np.mean(result2.VS)))
    axs[1][0].set_xlabel('radius(arcsec)')
    axs[1][0].set_ylabel('V/Sig')
    axs[1][1].plot(result.EffRadius, result.VS)
    axs[1][1].set_title('V/Sig = %.2f(%gRe),%.2f(%gRe),%.2f(%gRe)\nhas mean = %f'%(result.VS[effrad_index[0]] if effrad_index[0]<len(result.VS) else 9999,re_criterion_list[0],result.VS[effrad_index[1]] if effrad_index[1]<len(result.VS) else 9999,re_criterion_list[1],result.VS[effrad_index[2]] if effrad_index[2]<len(result.VS) else 9999 ,re_criterion_list[2],np.nanmean(result.VS)))
    axs[1][1].set_xlabel('radius(arcsec)')
    axs[1][1].set_ylabel('V/Sig')
    # or
    axs[2][0].plot(result2.EffRadius, result2.LambdaR)
    axs[2][0].set_title('Lambda = %.2f(%gRe),%.2f(%gRe),%.2f(%gRe)\nhas mean = %f'%(result2.LambdaR[effrad2_index[0]] if effrad2_index[0]<len(result2.LambdaR) else 9999,re_criterion_list[0],result2.LambdaR[effrad2_index[1]] if effrad2_index[1]<len(result2.LambdaR) else 9999,re_criterion_list[1],result2.LambdaR[effrad2_index[2]] if effrad2_index[2]<len(result2.LambdaR) else 9999 ,re_criterion_list[2],np.mean(result2.LambdaR)))
    axs[2][0].set_xlabel('radius(arcsec)')
    axs[2][0].set_ylabel('Lambda')
    axs[2][1].plot(result.EffRadius, result.LambdaR)
    axs[2][1].set_title('Lambda = %.2f(%gRe),%.2f(%gRe),%.2f(%gRe)\nhas mean = %f'%(result.LambdaR[effrad_index[0]] if effrad_index[0]<len(result.LambdaR) else 9999,re_criterion_list[0],result.LambdaR[effrad_index[1]] if effrad_index[1]<len(result.LambdaR) else 9999,re_criterion_list[1],result.LambdaR[effrad_index[2]] if effrad_index[2]<len(result.LambdaR) else 9999 ,re_criterion_list[2],np.nanmean(result.LambdaR)))
    axs[2][1].set_xlabel('radius(arcsec)')
    axs[2][1].set_ylabel('Lambda')
    plt.savefig('D:\\yenting\\results\\measureEps&Lambda\\stellar_e4(unmasked)\\ifu=%s_ra=%4.3f_dec=%4.3f_manga_result.png' %(list_plate[platenumber],list_ra[platenumber],list_dec[platenumber]))
    print('finish ifu=%s'%(list_plate[platenumber]))
    # plt.show()
    plt.clf()
    
    epspt2.append(round(eps[r_index[0]],3) if r_index[0]<len(eps) else 9999)
    epspt5.append(round(eps[r_index[1]],3) if r_index[1]<len(eps) else 9999)
    eps1.append(round(eps[r_index[2]],3) if r_index[2]<len(eps) else 9999)
    lambdapt2profile.append(round(result.LambdaR[effrad_index[0]],3) if effrad_index[0]<len(result.LambdaR) else 9999)
    lambdapt5profile.append(round(result.LambdaR[effrad_index[1]],3) if effrad_index[1]<len(result.LambdaR) else 9999)
    lambda1profile.append(round(result.LambdaR[effrad_index[2]],3) if effrad_index[2]<len(result.LambdaR) else 9999)
    lambdapt2default.append(round(result2.LambdaR[effrad_index[0]],3) if effrad_index[0]<len(result2.LambdaR) else 9999)
    lambdapt5default.append(round(result2.LambdaR[effrad_index[1]],3) if effrad_index[1]<len(result2.LambdaR) else 9999)
    lambda1default.append(round(result2.LambdaR[effrad_index[2]],3) if effrad_index[2]<len(result2.LambdaR) else 9999)
    table = Table([epspt2,lambdapt2profile,lambdapt2default,epspt5,lambdapt5profile,lambdapt5default,eps1,lambda1profile,lambda1default], names=['epspt2','lambdapt2(profile)','lambdapt2(default)','epspt5','lambdapt5(profile)','lambdapt5(default)','eps1','lambda1(profile)','lambda1(default)'])
    ascii.write(table, 'D:\\yenting\\results\\measureEps&Lambda\\stellar_e4(unmasked)\\Eps&Lambda.txt',overwrite=True)
    