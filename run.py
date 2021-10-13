import numpy as np
from astropy.io import fits
from astropy.table import Table, Column, MaskedColumn
from astropy.io import ascii
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import matplotlib.backend_bases as bb
from matplotlib.widgets import Slider, RangeSlider, Button, RadioButtons
import math
import csv
from measurejPA_class import jPA
from measurekPA_class import kPA
from datahandling import filehandling,match,comparePA,comparePAdiff,comparePAdiff_scatter,findPAdiff,newPAdiff,rewrite
# import montecarlo_class as mc
from montecarlo_analytic import analyticmc
from tqdm import tqdm as tqdm
from plotbin.display_pixels import display_pixels
import scipy.stats as stats
from itertools import zip_longest

#############################################################################
# initialize
#############################################################################
data = filehandling("500_e6.csv")
plateifu = np.array(data.extract('plateifu'))
my = data.extract('plateifu','hasjpa')
hasjpa = data.extract('hasjpa')
jPA_val = data.extract('jPA',tofloat=True)

#############################################################################
# compare with yenmei
#############################################################################

# kPA_source = 'STAR'
# kPA_range = '1'
# ympa = data.extract('PA_'+kPA_source+'_'+kPA_range+'RE',tofloat=True)
# tbcpa = data.extract('PA_app',tofloat=True)
# ympaerr = data.extract('PA_'+kPA_source+'_ERR_'+kPA_range+'RE',tofloat=True)
# tbcpaerr = data.extract('err_PA_app',tofloat=True)
# err = [90,25,20,15,10]
# comparePA,comparePAdiff(plateifu,ympa,ympaerr,tbcpa,tbcpaerr,err,err)

# str1 = 'lambda1(profile)(masked)'
# str2 = 'ym_lambda1'
# a = data.extract(str1,'hasjpa',tofloat=True)
# b = data.extract(str2,'hasjpa',tofloat=True)
# assert len(a)==len(b)
# an = [a[i] for i,_ in enumerate(a) if a[i]< 10 and b[i]<10]
# bn = [b[i] for i,_ in enumerate(b) if a[i]< 10 and b[i]<10]
# plt.scatter(an,bn)
# plt.plot([0,0.8],[0,0.8],color='red')
# plt.title(str1+' vs '+str2)
# plt.show()
# for i in my:
#     try:
#         print(i,plateifu.index(i)-2,ympa[plateifu.index(i)]) 
#         # pltgivenkPA(plateifu.index(i)-2,float(ympa[plateifu.index(i)]),float(ympaerr[plateifu.index(i)]))   
#         kPA()
#     except Exception:
#         pass

#############################################################################
# compare two monte carlo results
#############################################################################

# plt.figure(figsize=(10,10))
# ax1=plt.subplot(211)
# mc.to3d([30],modelingsize=300,ax=ax1)
# analyticmc([30],ax=ax1)
# # ax1.set_title(r'$\mathit{\phi=30^\circ}$', fontsize=20)
# # ax1.spines["right"].set_visible(False)
# # ax1.spines["top"].set_visible(False)
# plt.text(0.80, 0.90, r'$\mathit{\phi=30^\circ}$', fontsize=25, transform=ax1.transAxes)
# ax1.set_xlabel(r'$\mathit{\alpha(deg)}$', fontsize=20)
# ax1.set_ylabel(r'$\mathit{\mathbb{P}(\alpha|\phi)}$', fontsize=20)
# plt.setp(ax1, xticks=np.arange(0, 90, 10),yticks=np.arange(0, 0.041, 0.01))
# ax2=plt.subplot(212)
# mc.to3d([60],modelingsize=300,ax=ax2)
# analyticmc([60],ax=ax2)
# # ax2.set_title(r'$\mathit{\phi=60^\circ}$', fontsize=20)
# plt.text(0.80, 0.90, r'$\mathit{\phi=60^\circ}$', fontsize=25, transform=ax2.transAxes)
# ax2.set_xlabel(r'$\mathit{\alpha(deg)}$', fontsize=20)
# ax2.set_ylabel(r'$\mathit{\mathbb{P}(\alpha|\phi)}$', fontsize=20)
# plt.setp(ax2, xticks=np.arange(0, 91, 10),yticks=np.arange(0, 0.041, 0.01))
# # ax2.yaxis.set_visible(False)
# # plt.savefig(savedir+'allmc.png')
# plt.tight_layout()
# plt.show()
# plt.clf()

#############################################################################
# interactive data selection
#############################################################################

binsnum = 5
binsinit = np.linspace(0, 90, binsnum+1)
apply = [True,True,False,False,False,False,False,False,False,False,False]
result = None

def kinemetry_result(kPA_range):
    rad_data = filehandling("/Volumes/SDrive/yenting_pa_alignment/results/kinemetry/VOR10_1/result_rad.csv")
    pa_data = filehandling("/Volumes/SDrive/yenting_pa_alignment/results/kinemetry/VOR10_1/result_pa.csv")
    paerr_data = filehandling("/Volumes/SDrive/yenting_pa_alignment/results/kinemetry/VOR10_1/result_paerr.csv")
    k51_data = filehandling("/Volumes/SDrive/yenting_pa_alignment/results/kinemetry/VOR10_1/result_k51.csv")
    k51err_data = filehandling("/Volumes/SDrive/yenting_pa_alignment/results/kinemetry/VOR10_1/result_k51err.csv")
    kinemetry_rad,kinemetry_pa,kinemetry_paerr,kinemetry_k51,kinemetry_k51err=[],[],[],[],[]
    print('reading kinemetry data...')
    for plateifu in tqdm(plateifu):
        if plateifu!='' and data[(plateifu,'kinemetry','int')]==1:
            rad_kinemetry = [rad_data[(plateifu,str(i),'float')] for i in range(34)]
            pa_kinemetry = [pa_data[(plateifu,str(i),'float')]%180 for i in range(34)]
            paerr_kinemetry = [paerr_data[(plateifu,str(i),'float')] if paerr_data[(plateifu,str(i),'float')]<=180 else 180 for i in range(34)]
            k51_kinemetry = [k51_data[(plateifu,str(i),'float')] for i in range(34)]
            k51err_kinemetry = [k51err_data[(plateifu,str(i),'float')] for i in range(34)]
            index_range = abs(np.asarray(rad_kinemetry)-kPA_range*data[(plateifu,'NSA_ELPETRO_TH50_R','float')]).argmin()
            maxind_kinemetry = abs(np.asarray(rad_kinemetry)-data[(plateifu,'nomial FOV radius (arcsec)','float')]*np.sqrt(3)/2).argmin()
        if plateifu!='' and data[(plateifu,'kinemetry','int')]==1 and index_range<=maxind_kinemetry:
            kinemetry_rad.append(rad_kinemetry[index_range])
            kinemetry_pa.append(pa_kinemetry[index_range])
            kinemetry_paerr.append(paerr_kinemetry[index_range])
            kinemetry_k51.append(k51_kinemetry[index_range])
            kinemetry_k51err.append(k51err_kinemetry[index_range])
        else:
            kinemetry_rad.append(np.nan)
            kinemetry_pa.append(np.nan)
            kinemetry_paerr.append(np.nan)
            kinemetry_k51.append(np.nan)
            kinemetry_k51err.append(np.nan)

    return kinemetry_rad,kinemetry_pa,kinemetry_paerr,kinemetry_k51,kinemetry_k51err

def findscore(datanum,totdatanum,errmax,kPAerr,errlim,k51max,k51,k51lim,rad,kPAkin,kPAks):
    # return 2*datanum/totdatanum-errmax/30-k51max/0.20
    # return 3*datanum/totdatanum-errmax/30-k51max/0.20-np.nanmean(abs(kPAkin-kPAks))/30
    return datanum/totdatanum-errmax/errlim-k51max/k51lim-np.nanmean(abs(kPAkin-kPAks))/90

def optimize(row,kPA_range,binsnum,frombuffer=True,plot='35',ana=True ,directory='1011_2'):
    if not frombuffer:
        rad,kPAkin,kPAerr,k51,k51err=[],[],[],[],[]
        kPAks,pn=[],[]
        for i,e in enumerate(tqdm(plateifu)):
            if e!='':
                kslist=kPA(e, data,re_criterion_list=[float(kPA_range)],plot=False, source='STAR',dap='VOR10',para='NSA',measure='KS',binning=False,snthreshold=0)
                kinlist=kPA(e, data,re_criterion_list=[float(kPA_range)],plot=False, source='STAR',dap='VOR10',para='NSA',measure='Kinemetry',binning=False,snthreshold=0)
                kPAks.append(kslist[1]) 
                pn.append(kslist[5]) 
                rad.append(kinlist[1])
                kPAkin.append(kinlist[2])
                kPAerr.append(kinlist[3])
                k51.append(kinlist[4])
                k51err.append(kinlist[5])
            else:
                kPAks.append(np.nan) 
                pn.append(np.nan) 
                rad.append(np.nan)
                kPAkin.append(np.nan)
                kPAerr.append(np.nan)
                k51.append(np.nan)
                k51err.append(np.nan)
        table = Table([plateifu,rad,kPAkin,kPAerr,k51,k51err,kPAks,pn], names=['plateifu','rad','kPAkin','kPAerr','k51','k51err','kPAks','pn'])
        ascii.write(table, f'buffer{kPA_range}.csv',overwrite=True)

    if frombuffer:
        bufferdata=filehandling(f'buffer{kPA_range}.csv')
        rad=bufferdata.extract('rad',tofloat=True)
        kPAkin=bufferdata.extract('kPAkin',tofloat=True)
        kPAerr=bufferdata.extract('kPAerr',tofloat=True)
        k51=bufferdata.extract('k51',tofloat=True)
        k51err=bufferdata.extract('k51err',tofloat=True)
        kPAks=bufferdata.extract('kPAks',tofloat=True)
        pn=bufferdata.extract('pn',tofloat=True)

    
    errmaxlist,k51maxlist,scorelist=[],[],[]
    errlim=np.nanpercentile(kPAerr,100)
    # k51lim=np.nanpercentile(k51,100)
    k51lim=np.nanmax([e for i,e in enumerate(k51) if e/k51err[i]>3])
    for errmax in np.linspace(0,90,50):
        for k51max in np.linspace(0,1,50):
            PAdiff = findPAdiff(jPA_val,kPAkin)
            totdatanum = np.count_nonzero(~np.isnan(PAdiff))
            PAdiff = newPAdiff(kPAerr,0,errmax,PAdiff)
            PAdiff = newPAdiff(k51,0,k51max,PAdiff)
            PAdiff = newPAdiff(pn,30,np.inf,PAdiff)
            datanum = np.count_nonzero(~np.isnan(PAdiff))
            score = findscore(datanum,totdatanum,errmax,np.array(kPAerr),errlim,k51max,np.array(k51),k51lim,rad,np.array(kPAkin),np.array(kPAks))
            errmaxlist.append(errmax)
            k51maxlist.append(k51max)
            scorelist.append(score)
    opterr=errmaxlist[np.argmax(scorelist)]
    optk51=k51maxlist[np.argmax(scorelist)]

    PAdiff = findPAdiff(jPA_val,kPAkin)
    if apply[0]:PAdiff = newPAdiff(kPAerr,0,opterr,PAdiff)
    if apply[1]:PAdiff = np.array(newPAdiff(k51,0,optk51,PAdiff))
    datanum = np.count_nonzero(~np.isnan(PAdiff))
    write = csv.writer(open(f'/Volumes/SDrive/yenting_pa_alignment/results/for_yt/{directory}/all{kPA_range}.csv', 'w'))
    for values in zip(plateifu[~np.isnan(PAdiff)],PAdiff[~np.isnan(PAdiff)]):
        write.writerow(values)

    if '1' in plot:
        global fig1,axes1
        if row==0:fig1,axes1=plt.subplots(3,3,figsize=(4*3,10))
        axes1[0][row].scatter(kPAerr,k51)
        # plt.subplot(3,3,1*3+row-2).errorbar(kPAerr,k51,yerr=k51err,fmt='o')
        axes1[0][row].set_xlabel('kPA error')
        axes1[0][row].set_ylabel(r'$k_5/k_1$')
        axes1[1][row].hist(kPAerr)
        # axes1[0][row].axvline(np.nanpercentile(kPAerr,95),c='k')
        axes1[1][row].axvline(errlim,c='k')
        axes1[1][row].set_title('kPA err distribution')
        axes1[2][row].hist(k51)
        # axes1[2][row].axvline(np.nanpercentile(k51,95),c='k')
        axes1[2][row].axvline(k51lim,c='k')
        axes1[2][row].set_title('k5/k1 distribution')
        fig1.tight_layout()

    if '2' in plot:
        global fig2,axes2
        if row==0:fig2,axes2=plt.subplots(3,3,figsize=(4*3,10),sharey='row',sharex='col')
        axes2[0][row].set_title(f'{kPA_range} $R_e$', fontsize=15)
        axes2[1][row].text(0.40, 0.95, 'kPA error', transform=axes2[1][row].transAxes)
        axes2[0][0].text(-0.05,0.45,r'$k_5/k_1$',rotation=90, transform=axes2[0][row].transAxes)
        axes2[0][0].set_ylabel('score', fontsize=15)
        # axes2[0][0].yaxis.set_label_coords(-0.05,0.5)
        axes2[0][row].clabel(axes2[0][row].tricontour(errmaxlist,k51maxlist,scorelist,30), inline=True, fontsize=10)
        axes2[0][row].scatter(opterr,optk51)
        axes2[0][row].annotate(f'max at ({opterr:.2f},{optk51:.2f})',(opterr+.1,optk51+.03))
        axes2[2][row].text(0.40, 0.95, 'PA difference', transform=axes2[2][row].transAxes)
        axes2[1][0].set_ylabel('number of galaxies', fontsize=15)
        axes2[1][row].text(0.50, 0.85, '%d galaxies'%(datanum), fontsize=15, transform=axes2[1][row].transAxes)
        hist, bins = np.histogram([x for x in PAdiff if not np.isnan(x)],bins=binsinit)
        axes2[1][row].bar(bins[:-1]+45/binsnum,hist, width=90/binsnum,color='green',alpha=0.6)
        axes2[2][0].set_ylabel('3D angle probability', fontsize=15)
        if ana:analyticmc(PAdiff,ax=axes2[2][row])
        axes2[2][row].set_xlabel('angle',fontsize=13)
        fig2.tight_layout()
        fig2.subplots_adjust(hspace=0,wspace=0)

    '''
    DAP data
    35:
    ELG17   = 'OIII-5008'
    ELG24   = 'Ha-6564 '
    ELG12   = 'Hdel-4102'
    46:
    SPI22   = 'HDeltaA '
    SPI24   = 'HDeltaF '
    SPI45	= 'Dn4000  '
    '''
    h=0.7
    lumdist = np.array(data.extract('LDIST_Z',tofloat=True))
    stellarmass = np.log10(np.array(data.extract('NSA_ELPETRO_MASS',tofloat=True))*h**2)
    stellar_sigma_1re = data.extract('STELLAR_SIGMA_1RE',tofloat=True)
    o3_gflux_1re = np.array([float(e.split(',')[16]) if e!='' else np.nan for e in data.extract('EMLINE_GFLUX_1RE')])
    o3_lum_1re = np.array([np.log10(4*np.pi*f*1e-17*(d/0.7*(3.085678*1e24))**2) for (f,d) in zip(o3_gflux_1re,lumdist)])
    dn4000_specindex_1re = np.array([float(e.split(',')[44]) if e!='' else np.nan for e in data.extract('SPECINDEX_1RE')])
    hd_specindex_1re = np.array([float(e.split(',')[21]) if e!='' else np.nan for e in data.extract('SPECINDEX_1RE')]) #HDeltaA
    sfr_1re = data.extract('SFR_1RE',tofloat=True)
    ha_gew_1re = np.array([float(e.split(',')[23]) if e!='' else np.nan for e in data.extract('EMLINE_GEW_1RE')]) #######################!!!!!!!!!!!!!!!!!!!!!!!!!
    stellarmass_r,stellar_sigma_r,o3_lum_r,dn4000_specindex_r,hd_specindex_r,sfr_r = [],[],[],[],[],[]
    for i,e in enumerate(tqdm(plateifu)):
        if e!='':
            print('ji')
            index = kPA(e, data,re_criterion_list=[float(kPA_range)],plot=False, source='STAR',dap='VOR10',para='NSA',measure='None',binning=False,snthreshold=0)[-1]
            if e not in ['8479-12705','8479-9101','9183-3701']:
                stellarmass_r.append(np.log10(np.sum(10**fits.open(f"/Volumes/SDrive/yenting_pa_alignment/MaNGA/Pipe3D/manga-{e}.Pipe3D.cube.fits.gz")['SSP'].data[19][index]))) #stellar mass density dust corrected in 'm_Sun/spaxels^2'
            else: stellarmass_r.append(np.nan)
            hdu = fits.open(f"/Volumes/SDrive/yenting_pa_alignment/MaNGA/VOR10-MILESHC-MASTARSSP/manga-{e}-MAPS-VOR10-MILESHC-MASTARSSP.fits.gz")
            stellar_sigma_r.append(np.average((hdu['STELLAR_VEL'].data[index]-np.average(hdu['STELLAR_VEL'].data[index], weights=hdu['BIN_MFLUX'].data[index]))**2, weights=hdu['BIN_MFLUX'].data[index]))
            o3_gflux_r = np.sum(hdu['EMLINE_GFLUX'].data[16][index]) #Oiii-5008
            o3_lum_r.append(np.log10(np.sum(4*np.pi*hdu['EMLINE_GFLUX'].data[16][index]*1e-17*(lumdist[i]/0.7*(3.085678*1e24))**2)))
            dn4000_specindex_r.append(np.mean(hdu['SPECINDEX'].data[44][index])) #Dn4000
            hd_specindex_r.append(np.mean(hdu['SPECINDEX'].data[21][index])) #HDeltaA
            sfr_r.append(np.sum(10**(np.log10(4*np.pi*hdu['EMLINE_GFLUX'].data[23][index]*1e-17*(lumdist[i]/0.7*(3.085678*1e24))**2)-41.27))) #Ha-6564
        else: 
            stellarmass_r.append(np.nan)
            stellar_sigma_r.append(np.nan)
            o3_lum_r.append(np.nan)
            dn4000_specindex_r.append(np.nan)
            hd_specindex_r.append(np.nan)
            sfr_r.append(np.nan)
    
    ssfr = np.array([np.log10(sfr/sm) for (sfr,sm) in zip(sfr_1re,data.extract('NSA_ELPETRO_MASS',tofloat=True))])
    surface_mass_density = np.array([np.log10(mass/np.pi/re**2) for (mass,re) in zip(data.extract('NSA_ELPETRO_MASS',tofloat=True),data.extract('NSA_ELPETRO_TH50_R',tofloat=True))])
    blackholemass = np.array([10**(8.13+4.02*np.log10(sig/200))/1e8 for (lum,sig) in zip(o3_lum_1re,stellar_sigma_1re)])
    oer = np.array([np.log10(10**(lum)/(1.28*1e46*bhm)) for (lum,bhm) in zip(o3_lum_1re,blackholemass)])
    ssfr_r = np.array([np.log10(sfr/sm) for (sfr,sm) in zip(sfr_r,10**(np.array(stellarmass_r)))])
    surface_mass_density_r = np.array([np.log10(mass/np.pi/re**2) for (mass,re) in zip(10**(np.array(stellarmass_r)),data.extract('NSA_ELPETRO_TH50_R',tofloat=True))])
    blackholemass_r = np.array([10**(8.13+4.02*np.log10(sig/200))/1e8 for (lum,sig) in zip(o3_lum_r,stellar_sigma_r)])
    oer_r = np.array([np.log10(10**(lum)/(1.28*1e46*bhm)) for (lum,bhm) in zip(o3_lum_r,blackholemass_r)])
    
    nsa_sersic_n = np.array([e if e!=-9999.0 else np.nan for e in data.extract('NSA_SERSIC_N',tofloat=True)])
    radio_lum = np.array(data.extract('lum',tofloat=True))
    radio_morphology = np.array(data.extract('radio_morphology',tofloat=True))

    halomass = data.extract('Mh_L',tofloat=True)
    objgp = data.extract('objgp(arcsec)',tofloat=True)
    r180 = data.extract('r_200',tofloat=True)
    distance = [o/r if not (np.isnan(o) or np.isnan(r)) else np.nan for (o,r) in zip(objgp,r180)]
    numofgal = data.extract('300kpc 1500 km/s',tofloat=True)
    nearest = data.extract('5th nearest (mpc) by cylindrical method cut = 1000 km/s',tofloat=True)
    bcg = [x if x!=2 else 0.0 for x in data.extract('BCG',tofloat=True)]
    mmg = [x if x!=2 else 0.0 for x in data.extract('MMG',tofloat=True)]

    
    if '3' in plot:
        global fig3,axes3
        if row==0:fig3,axes3=plt.subplots(4,3,figsize=(3.5*3,10),sharey='row',sharex='col')
        cut = np.nanmedian([e for i,e in enumerate(stellarmass) if not np.isnan(PAdiff[i])])
        PAdiff_original=np.copy(PAdiff)
        PAdiff_copy1=np.copy(PAdiff)
        PAdiff_copy1 = newPAdiff(stellarmass,0,cut+1e-3,PAdiff_copy1)
        datanum = np.count_nonzero(~np.isnan(PAdiff_copy1))
        hist, bins = np.histogram([x for x in PAdiff_copy1 if not np.isnan(x)],bins=binsinit)
        axes3[0][row].set_title(f'{kPA_range} $R_e$', fontsize=15)
        axes3[0][row].bar(bins[:-1]+45/binsnum,hist, width=90/binsnum,color='green',alpha=0.6)
        # axes3[0][row].set_title(f'low stellar mass < {cut:.2f}\n{datanum} galaxies')
        axes3[0][row].text(0.95,0.95, f'low stellar mass < {cut:.2f}\n{datanum} galaxies', fontsize=13, transform=axes3[0][row].transAxes, horizontalalignment='left',verticalalignment='bottom')
        axes3[0][0].set_ylabel('number of galaxies',fontsize=13)
        axes3[1][0].set_ylabel('3D angle probability',fontsize=13)
        if ana:analyticmc(PAdiff_copy1,ax=axes3[1][row])
        
        PAdiff_copy2=np.copy(PAdiff)
        PAdiff_copy2 = newPAdiff(stellarmass,cut+1e-3,np.inf,PAdiff_copy2)
        datanum = np.count_nonzero(~np.isnan(PAdiff_copy2))
        hist, bins = np.histogram([x for x in PAdiff_copy2 if not np.isnan(x)],bins=binsinit)
        axes3[2][row].bar(bins[:-1]+45/binsnum,hist, width=90/binsnum,color='green',alpha=0.6)
        # axes3[2][row].set_title(f'high stellar mass > {cut:.2f}\n{datanum} galaxies')
        axes3[2][row].text(0.95,0.95, f'high stellar mass > {cut:.2f}\n{datanum} galaxies', fontsize=13, transform=axes3[2][row].transAxes, horizontalalignment='left',verticalalignment='bottom')
        axes3[2][0].set_ylabel('number of galaxies',fontsize=13)
        axes3[3][row].set_xlabel('angle',fontsize=13)
        axes3[3][0].set_ylabel('3D angle probability',fontsize=13)
        if ana:analyticmc(PAdiff_copy2,ax=axes3[3][row])
        fig3.tight_layout()
        fig3.subplots_adjust(hspace=0,wspace=0)

        ## for debugging
        # if row==0:
        #     global fig3_1,axes3_1
        #     fig3_1,axes3_1=plt.subplots(3,3,figsize=(4*3,7))
        # comparePAdiff(PAdiff_copy1, PAdiff_copy2, 'low stellar mass', 'high stellar mass', f'{np.count_nonzero(~np.isnan(PAdiff_copy1))}, {np.count_nonzero(~np.isnan(PAdiff_copy2))}', stellarmass, ax=axes3_1[0][row])
        # comparePAdiff(PAdiff_original, PAdiff_copy2, 'low stellar mass', 'high stellar mass', f'{np.count_nonzero(~np.isnan(PAdiff_original))}, {np.count_nonzero(~np.isnan(PAdiff_copy2))}', stellarmass, ax=axes3_1[1][row])
        # comparePAdiff(PAdiff_copy1, PAdiff_original, 'low stellar mass', 'high stellar mass', f'{np.count_nonzero(~np.isnan(PAdiff_copy1))}, {np.count_nonzero(~np.isnan(PAdiff_original))}', stellarmass, ax=axes3_1[2][row])
        
        if kPA_range=='0.3':
            global fig3_1,axes3_1
            fig3_1,axes3_1=plt.subplots(3,3,figsize=(5*3,10))
            write = csv.writer(open(f'/Volumes/SDrive/yenting_pa_alignment/results/for_yt/{directory}/stellarmass{kPA_range}.csv', 'w'))
            write.writerow(['low stellar mass','PA difference','high stellar mass','PA difference'])
            for values in zip_longest(*[plateifu[~np.isnan(PAdiff_copy1)],PAdiff_copy1[~np.isnan(PAdiff_copy1)],plateifu[~np.isnan(PAdiff_copy2)],PAdiff_copy2[~np.isnan(PAdiff_copy2)]]):
                write.writerow(values)
            comparePAdiff(PAdiff_copy1, PAdiff_copy2, 'low stellar mass', 'high stellar mass', 'sersic index n', nsa_sersic_n, ax=axes3_1[0][0])
            comparePAdiff(PAdiff_copy1, PAdiff_copy2, 'low stellar mass', 'high stellar mass', 'log(sSFR)',ssfr, ax=axes3_1[0][1])
            comparePAdiff(PAdiff_copy1, PAdiff_copy2, 'low stellar mass', 'high stellar mass', 'Halpha EW',ha_gew_1re, ax=axes3_1[0][2])
            comparePAdiff(PAdiff_copy1, PAdiff_copy2, 'low stellar mass', 'high stellar mass', 'stellar surface mass density',surface_mass_density, ax=axes3_1[1][0])
            comparePAdiff(PAdiff_copy1, PAdiff_copy2, 'low stellar mass', 'high stellar mass', 'OER',oer, ax=axes3_1[1][1])
            comparePAdiff(PAdiff_copy1, PAdiff_copy2, 'low stellar mass', 'high stellar mass', 'radio_morphology',radio_morphology, ax=axes3_1[1][2],binnum=4, setticks=[])
            axes3_1[1][2].text(0.4,-0.05, 'NaN', fontsize=13, horizontalalignment='center',verticalalignment='top')
            axes3_1[1][2].text(1.15,-0.05, 'FR1', fontsize=13, horizontalalignment='center',verticalalignment='top')
            axes3_1[1][2].text(1.85,-0.05, 'FR2', fontsize=13, horizontalalignment='center',verticalalignment='top')
            axes3_1[1][2].text(2.6,-0.05, 'NAT', fontsize=13, horizontalalignment='center',verticalalignment='top')
            comparePAdiff_scatter(PAdiff_copy1, PAdiff_copy2, 'low stellar mass', 'high stellar mass', '', dn4000_specindex_1re, hd_specindex_1re, r'D$_n$(4000)', r'HDelta$_A$', axes3_1[2][0])
            comparePAdiff_scatter(PAdiff_copy1, PAdiff_copy2, 'low stellar mass', 'high stellar mass', '', o3_lum_1re, radio_lum, r'log(O[III] luminosity)', r'log($P_{1.4GHz}$)', axes3_1[2][1])
            axes3_1[2][2].axis('off')
            fig3_1.tight_layout()

    if '4' in plot:
        global fig4,axes4
        if row==0:fig4,axes4=plt.subplots(4,3,figsize=(4*3,10),sharey='row',sharex='col')
        cut = np.nanmedian([e for i,e in enumerate(halomass) if not np.isnan(PAdiff[i])])
        PAdiff_copy=np.copy(PAdiff)
        PAdiff_copy = newPAdiff(halomass,0,cut,PAdiff_copy)
        datanum = np.count_nonzero(~np.isnan(PAdiff_copy))
        hist, bins = np.histogram([x for x in PAdiff_copy if not np.isnan(x)],bins=binsinit)
        axes4[0][row].set_title(f'{kPA_range} $R_e$', fontsize=15)
        axes4[0][row].bar(bins[:-1]+45/binsnum,hist, width=90/binsnum,color='green',alpha=0.6)
        # axes4[0][row].set_title(f'low halo mass < {cut:.2f}\n{datanum} galaxies')
        axes4[0][row].text(0.95,0.95, f'low halo mass < {cut:.2f}\n{datanum} galaxies', fontsize=13, transform=axes4[0][row].transAxes, horizontalalignment='right',verticalalignment='top')
        axes4[0][0].set_ylabel('number of galaxies',fontsize=13)
        axes4[1][0].set_ylabel('3D angle probability',fontsize=13)
        if ana:analyticmc(PAdiff_copy,ax=axes4[1][row])
        
        PAdiff_copy=np.copy(PAdiff)
        PAdiff_copy = newPAdiff(halomass,cut,np.inf,PAdiff_copy)
        datanum = np.count_nonzero(~np.isnan(PAdiff_copy))
        hist, bins = np.histogram([x for x in PAdiff_copy if not np.isnan(x)],bins=binsinit)
        axes4[2][row].bar(bins[:-1]+45/binsnum,hist, width=90/binsnum,color='green',alpha=0.6)
        # axes4[2][row].set_title(f'high halo mass > {cut:.2f}\n{datanum} galaxies')
        axes4[2][row].text(0.95,0.95, f'high halo mass > {cut:.2f}\n{datanum} galaxies', fontsize=13, transform=axes4[2][row].transAxes, horizontalalignment='right',verticalalignment='top')
        axes4[2][0].set_ylabel('number of galaxies',fontsize=13)
        axes4[3][row].set_xlabel('angle',fontsize=13)
        axes4[3][0].set_ylabel('3D angle probability',fontsize=13)
        if ana:analyticmc(PAdiff_copy,ax=axes4[3][row])
        fig4.tight_layout()
        fig4.subplots_adjust(hspace=0,wspace=0)

    if '5' in plot:
        global fig5,axes5
        if row==0:fig5,axes5=plt.subplots(4,3,figsize=(4*3,10),sharey='row',sharex='col')
        cut = np.nanmedian([e for i,e in enumerate(stellar_sigma_1re) if not np.isnan(PAdiff[i])])
        PAdiff_copy1=np.copy(PAdiff)
        PAdiff_copy1 = newPAdiff(stellar_sigma_1re,0,cut,PAdiff_copy1)
        datanum = np.count_nonzero(~np.isnan(PAdiff_copy1))
        hist, bins = np.histogram([x for x in PAdiff_copy1 if not np.isnan(x)],bins=binsinit)
        axes5[0][row].set_title(f'{kPA_range} $R_e$', fontsize=15)
        axes5[0][row].bar(bins[:-1]+45/binsnum,hist, width=90/binsnum,color='green',alpha=0.6)
        # axes5[0][row].set_title(f'low stellar velocity\ndispersion < {cut:.2f}\n{datanum} galaxies')
        axes5[0][row].text(0.95,0.95, f'low stellar velocity\ndispersion < {cut:.2f}km/s\n{datanum} galaxies', fontsize=13, transform=axes5[0][row].transAxes, horizontalalignment='right',verticalalignment='top')
        axes5[0][0].set_ylabel('number of galaxies',fontsize=13)
        axes5[1][0].set_ylabel('3D angle probability',fontsize=13)
        if ana:analyticmc(PAdiff_copy1,ax=axes5[1][row])
        
        PAdiff_copy2=np.copy(PAdiff)
        PAdiff_copy2 = newPAdiff(stellar_sigma_1re,cut,np.inf,PAdiff_copy2)
        datanum = np.count_nonzero(~np.isnan(PAdiff_copy2))
        hist, bins = np.histogram([x for x in PAdiff_copy2 if not np.isnan(x)],bins=binsinit)
        axes5[2][row].bar(bins[:-1]+45/binsnum,hist, width=90/binsnum,color='green',alpha=0.6)
        # axes5[2][row].set_title(f'high stellar velocity\ndispersion > {cut:.2f}\n{datanum} galaxies')
        axes5[2][row].text(0.95,0.95, f'high stellar velocity\ndispersion > {cut:.2f}km/s\n{datanum} galaxies', fontsize=13, transform=axes5[2][row].transAxes, horizontalalignment='right',verticalalignment='top')
        axes5[2][0].set_ylabel('number of galaxies',fontsize=13)
        axes5[3][row].set_xlabel('angle',fontsize=13)
        axes5[3][0].set_ylabel('3D angle probability',fontsize=13)
        if ana:analyticmc(PAdiff_copy2,ax=axes5[3][row])
        fig5.tight_layout()
        fig5.subplots_adjust(hspace=0,wspace=0)

        if kPA_range=='0.3':
            global fig5_1,axes5_1
            fig5_1,axes5_1=plt.subplots(3,3,figsize=(5*3,10))
            write = csv.writer(open(f'/Volumes/SDrive/yenting_pa_alignment/results/for_yt/{directory}/stellarveldispersion{kPA_range}.csv', 'w'))
            write.writerow(['low stellar velocity dispersion','PA difference','high stellar velocity dispersion','PA difference'])
            for values in zip_longest(*[plateifu[~np.isnan(PAdiff_copy1)],PAdiff_copy1[~np.isnan(PAdiff_copy1)],plateifu[~np.isnan(PAdiff_copy2)],PAdiff_copy2[~np.isnan(PAdiff_copy2)]]):
                write.writerow(values)
            write.writerow('low stellar velocity dispersion')
            write.writerow(plateifu[~np.isnan(PAdiff_copy1)])
            write.writerow(PAdiff_copy1[~np.isnan(PAdiff_copy1)])
            write.writerow('high stellar velocity dispersion')
            write.writerow(plateifu[~np.isnan(PAdiff_copy2)])
            write.writerow(PAdiff_copy2[~np.isnan(PAdiff_copy2)])
            comparePAdiff(PAdiff_copy1, PAdiff_copy2, 'low velocity dispersion', 'high velocity dispersion', 'sersic index n', nsa_sersic_n, ax=axes5_1[0][0])
            comparePAdiff(PAdiff_copy1, PAdiff_copy2, 'low velocity dispersion', 'high velocity dispersion', 'log(sSFR)',ssfr, ax=axes5_1[0][1])
            comparePAdiff(PAdiff_copy1, PAdiff_copy2, 'low velocity dispersion', 'high velocity dispersion', 'Halpha EW',ha_gew_1re, ax=axes5_1[0][2])
            comparePAdiff(PAdiff_copy1, PAdiff_copy2, 'low velocity dispersion', 'high velocity dispersion', 'stellar surface mass density',surface_mass_density, ax=axes5_1[1][0])
            comparePAdiff(PAdiff_copy1, PAdiff_copy2, 'low velocity dispersion', 'high velocity dispersion', 'OER',oer, ax=axes5_1[1][1])
            comparePAdiff(PAdiff_copy1, PAdiff_copy2, 'low velocity dispersion', 'high velocity dispersion', 'radio_morphology',radio_morphology, ax=axes5_1[1][2],binnum=4, setticks=[])
            axes5_1[1][2].text(0.4,-0.05, 'NaN', fontsize=13, horizontalalignment='center',verticalalignment='top')
            axes5_1[1][2].text(1.15,-0.05, 'FR1', fontsize=13, horizontalalignment='center',verticalalignment='top')
            axes5_1[1][2].text(1.85,-0.05, 'FR2', fontsize=13, horizontalalignment='center',verticalalignment='top')
            axes5_1[1][2].text(2.6,-0.05, 'NAT', fontsize=13, horizontalalignment='center',verticalalignment='top')
            comparePAdiff_scatter(PAdiff_copy1, PAdiff_copy2, 'low velocity dispersion', 'high velocity dispersion', '', dn4000_specindex_1re, hd_specindex_1re, r'D$_n$(4000)', r'HDelta$_A$', axes5_1[2][0])
            comparePAdiff_scatter(PAdiff_copy1, PAdiff_copy2, 'low velocity dispersion', 'high velocity dispersion', '', o3_lum_1re, radio_lum, r'log(O[III] luminosity)', r'log($P_{1.4GHz}$)', axes5_1[2][1])
            axes5_1[2][2].axis('off')
            fig5_1.tight_layout()

    if '6' in plot and kPA_range=='0.3':
        global fig6,axes6
        fig6,axes6=plt.subplots(4,3,figsize=(4*3,10),sharey='row',sharex='col')
        cut = np.nanmedian([e for i,e in enumerate(distance) if not np.isnan(PAdiff[i])])
        PAdiff_copy=np.copy(PAdiff)
        PAdiff_copy = newPAdiff(distance,cut,np.inf,PAdiff_copy)
        datanum = np.count_nonzero(~np.isnan(PAdiff_copy))
        hist, bins = np.histogram([x for x in PAdiff_copy if not np.isnan(x)],bins=binsinit)
        axes6[0][0].set_title(f'{kPA_range} $R_e$', fontsize=15)
        axes6[0][0].bar(bins[:-1]+45/binsnum,hist, width=90/binsnum,color='green',alpha=0.6)
        axes6[0][0].text(0.95,0.95, 'distance to group\ncenter>%.2f$R_{180}$\n%d galaxies'%(cut,datanum), fontsize=13, transform=axes6[0][0].transAxes, horizontalalignment='right',verticalalignment='top')
        axes6[0][0].set_ylabel('number of galaxies',fontsize=13)
        axes6[1][0].set_ylabel('3D angle probability',fontsize=13)
        if ana:analyticmc(PAdiff_copy,ax=axes6[1][0])

        PAdiff_copy=np.copy(PAdiff)
        PAdiff_copy = newPAdiff(distance,0,cut,PAdiff_copy)
        datanum = np.count_nonzero(~np.isnan(PAdiff_copy))
        hist, bins = np.histogram([x for x in PAdiff_copy if not np.isnan(x)],bins=binsinit)
        axes6[2][0].bar(bins[:-1]+45/binsnum,hist, width=90/binsnum,color='green',alpha=0.6)
        axes6[2][0].text(0.95,0.95, 'distance to group\ncenter<%.2f$R_{180}$\n%d galaxies'%(cut,datanum), fontsize=13, transform=axes6[2][0].transAxes, horizontalalignment='right',verticalalignment='top')
        axes6[2][0].set_ylabel('number of galaxies',fontsize=13)
        axes6[3][0].set_ylabel('3D angle probability',fontsize=13)
        axes6[3][0].set_xlabel('angle',fontsize=13)
        if ana:analyticmc(PAdiff_copy,ax=axes6[3][0])

        cut = np.nanmedian([e for i,e in enumerate(numofgal) if not np.isnan(PAdiff[i])])
        PAdiff_copy=np.copy(PAdiff)
        PAdiff_copy = newPAdiff(numofgal,0,cut,PAdiff_copy)
        datanum = np.count_nonzero(~np.isnan(PAdiff_copy))
        hist, bins = np.histogram([x for x in PAdiff_copy if not np.isnan(x)],bins=binsinit)
        axes6[0][1].set_title(f'{kPA_range} $R_e$', fontsize=15)
        axes6[0][1].bar(bins[:-1]+45/binsnum,hist, width=90/binsnum,color='green',alpha=0.6)
        axes6[0][1].text(0.95,0.95, f'number of galaxies in\ncylindrical region < {cut:n}\n{datanum} galaxies', fontsize=13, transform=axes6[0][1].transAxes, horizontalalignment='right',verticalalignment='top')
        if ana:analyticmc(PAdiff_copy,ax=axes6[1][1])

        PAdiff_copy=np.copy(PAdiff)
        PAdiff_copy = newPAdiff(numofgal,cut,np.inf,PAdiff_copy)
        datanum = np.count_nonzero(~np.isnan(PAdiff_copy))
        hist, bins = np.histogram([x for x in PAdiff_copy if not np.isnan(x)],bins=binsinit)
        axes6[2][1].bar(bins[:-1]+45/binsnum,hist, width=90/binsnum,color='green',alpha=0.6)
        axes6[2][1].text(0.95,0.95, f'number of galaxies in\ncylindrical region > {cut:n}\n{datanum} galaxies', fontsize=13, transform=axes6[2][1].transAxes, horizontalalignment='right',verticalalignment='top')
        axes6[3][1].set_xlabel('angle',fontsize=13)
        if ana:analyticmc(PAdiff_copy,ax=axes6[3][1])

        cut = np.nanmedian([e for i,e in enumerate(nearest) if not np.isnan(PAdiff[i])])
        PAdiff_copy=np.copy(PAdiff)
        PAdiff_copy = newPAdiff(nearest,cut,np.inf,PAdiff_copy)
        datanum = np.count_nonzero(~np.isnan(PAdiff_copy))
        hist, bins = np.histogram([x for x in PAdiff_copy if not np.isnan(x)],bins=binsinit)
        axes6[0][2].set_title(f'{kPA_range} $R_e$', fontsize=15)
        axes6[0][2].bar(bins[:-1]+45/binsnum,hist, width=90/binsnum,color='green',alpha=0.6)
        axes6[0][2].text(0.95,0.95, f'distance to 5th nearest\ngalaxy > {cut:.2f}kpc\n{datanum} galaxies', fontsize=13, transform=axes6[0][2].transAxes, horizontalalignment='right',verticalalignment='top')
        if ana:analyticmc(PAdiff_copy,ax=axes6[1][2])

        PAdiff_copy=np.copy(PAdiff)
        PAdiff_copy = newPAdiff(nearest,0,cut,PAdiff_copy)
        datanum = np.count_nonzero(~np.isnan(PAdiff_copy))
        hist, bins = np.histogram([x for x in PAdiff_copy if not np.isnan(x)],bins=binsinit)
        axes6[2][2].bar(bins[:-1]+45/binsnum,hist, width=90/binsnum,color='green',alpha=0.6)
        axes6[2][2].text(0.95,0.95, f'distance to 5th nearest\ngalaxy < {cut:.2f}kpc\n{datanum} galaxies', fontsize=13, transform=axes6[2][2].transAxes, horizontalalignment='right',verticalalignment='top')
        axes6[3][2].set_xlabel('angle',fontsize=13)
        if ana:analyticmc(PAdiff_copy,ax=axes6[3][2])
        fig6.tight_layout()
        fig6.subplots_adjust(hspace=0,wspace=0)

    if '7' in plot:
        global fig7,axes7
        if row==0:fig7,axes7=plt.subplots(4,3,figsize=(4*3,10),sharey='row',sharex='col')
        if kPA_range=='0.3':
            global fig7_1,axes7_1
            fig7_1,axes7_1=plt.subplots(2,2,figsize=(3.5*2,4),sharey='row',sharex='col')
        PAdiff_copy=np.copy(PAdiff)
        PAdiff_copy = newPAdiff(bcg,-0.5,0.5,PAdiff_copy)
        hist, bins = np.histogram([x for x in PAdiff_copy if not np.isnan(x)],bins=binsinit)
        axes7[0][row].set_title(f'{kPA_range} $R_e$', fontsize=15)
        axes7[0][0].set_ylabel('number of galaxies',fontsize=13)
        axes7[0][row].bar(bins[:-1]+45/binsnum,hist, width=90/binsnum,color='green',alpha=0.6)
        axes7[0][row].text(0.95,0.95, f'not BCG', fontsize=13, transform=axes7[0][row].transAxes, horizontalalignment='right',verticalalignment='top')
        if kPA_range=='0.3' and ana:
            axes7_1[0][0].text(0.05,0.95, f'not BCG', fontsize=13, transform=axes7_1[0][0].transAxes, horizontalalignment='left',verticalalignment='top')
            analyticmc(PAdiff_copy,ax=axes7_1[0][0])   

        PAdiff_copy=np.copy(PAdiff)
        PAdiff_copy = newPAdiff(bcg,0.5,1.5,PAdiff_copy)
        hist, bins = np.histogram([x for x in PAdiff_copy if not np.isnan(x)],bins=binsinit)
        axes7[1][0].set_ylabel('number of galaxies',fontsize=13)
        axes7[1][row].bar(bins[:-1]+45/binsnum,hist, width=90/binsnum,color='green',alpha=0.6)
        axes7[1][row].text(0.95,0.95, f'BCG', fontsize=13, transform=axes7[1][row].transAxes, horizontalalignment='right',verticalalignment='top')
        if kPA_range=='0.3' and ana:
            axes7_1[1][0].text(0.05,0.95, f'BCG', fontsize=13, transform=axes7_1[1][0].transAxes, horizontalalignment='left',verticalalignment='top')
            analyticmc(PAdiff_copy,ax=axes7_1[1][0])   

        PAdiff_copy=np.copy(PAdiff)
        PAdiff_copy = newPAdiff(mmg,-0.5,0.5,PAdiff_copy)
        hist, bins = np.histogram([x for x in PAdiff_copy if not np.isnan(x)],bins=binsinit)
        axes7[2][0].set_ylabel('number of galaxies',fontsize=13)
        axes7[2][row].bar(bins[:-1]+45/binsnum,hist, width=90/binsnum,color='green',alpha=0.6)
        axes7[2][row].text(0.95,0.95, f'not MMG', fontsize=13, transform=axes7[2][row].transAxes, horizontalalignment='right',verticalalignment='top')
        if kPA_range=='0.3' and ana:
            axes7_1[0][1].text(0.05,0.95, f'not MMG', fontsize=13, transform=axes7_1[0][1].transAxes, horizontalalignment='left',verticalalignment='top')
            analyticmc(PAdiff_copy,ax=axes7_1[0][1])   

        PAdiff_copy=np.copy(PAdiff)
        PAdiff_copy = newPAdiff(mmg,0.5,1.5,PAdiff_copy)
        hist, bins = np.histogram([x for x in PAdiff_copy if not np.isnan(x)],bins=binsinit)
        axes7[3][0].set_ylabel('number of galaxies',fontsize=13)
        axes7[3][row].bar(bins[:-1]+45/binsnum,hist, width=90/binsnum,color='green',alpha=0.6)
        axes7[3][row].text(0.95,0.95, f'MMG', fontsize=13, transform=axes7[3][row].transAxes, horizontalalignment='right',verticalalignment='top')
        axes7[3][row].set_xlabel('angle',fontsize=13)
        if kPA_range=='0.3' and ana:
            axes7_1[1][1].text(0.05,0.95, f'MMG', fontsize=13, transform=axes7_1[1][1].transAxes, horizontalalignment='left',verticalalignment='top')
            analyticmc(PAdiff_copy,ax=axes7_1[1][1]) 
            fig7_1.suptitle(r'3D PDF at 0.3$R_e$')  
            fig7_1.supxlabel('angle')  
            fig7_1.supylabel('probability')  
            # fig7_1.tight_layout()
            fig7_1.subplots_adjust(hspace=0,wspace=0)
        fig7.tight_layout()
        fig7.subplots_adjust(hspace=0,wspace=0)


        
    if '8' in plot:
        global fig8,axes8
        if row==0:fig8,axes8=plt.subplots(3,9,figsize=(21,10),sharex='col')
        axes8[row][0].scatter(PAdiff,stellarmass)
        axes8[0][0].set_title('stellarmass')
        axes8[row][1].scatter(PAdiff,halomass)
        axes8[0][1].set_title('halomass')
        axes8[row][2].scatter(PAdiff,distance)
        axes8[row][2].set_ylim(0,sorted(distance)[-3])
        axes8[0][2].set_title('distance')
        axes8[row][3].scatter(PAdiff,numofgal)
        axes8[0][3].set_title('numofgal')
        axes8[row][4].scatter(PAdiff,nearest)
        axes8[0][4].set_title('nearest')
        axes8[row][5].scatter(PAdiff,bcg)
        axes8[0][5].set_title('bcg')
        axes8[row][6].scatter(PAdiff,mmg)
        axes8[0][6].set_title('mmg')
        axes8[row][7].scatter(PAdiff,pn)
        axes8[0][7].set_title('pixel number')
        axes8[row][8].scatter(PAdiff,stellar_sigma_1re)
        axes8[0][8].set_title('velocity dispersion')
        fig8.tight_layout()
        fig8.subplots_adjust(hspace=0,wspace=0)
    
    if row==2:
        if '1' in plot:fig1.savefig(f'/Volumes/SDrive/yenting_pa_alignment/results/for_yt/{directory}/figure1.png')
        if '2' in plot:fig2.savefig(f'/Volumes/SDrive/yenting_pa_alignment/results/for_yt/{directory}/figure2.png')
        if '3' in plot:fig3.savefig(f'/Volumes/SDrive/yenting_pa_alignment/results/for_yt/{directory}/figure3.png')
        if '3' in plot:fig3_1.savefig(f'/Volumes/SDrive/yenting_pa_alignment/results/for_yt/{directory}/figure3_1.png')
        if '4' in plot:fig4.savefig(f'/Volumes/SDrive/yenting_pa_alignment/results/for_yt/{directory}/figure4.png')
        if '5' in plot:fig5.savefig(f'/Volumes/SDrive/yenting_pa_alignment/results/for_yt/{directory}/figure5.png')
        if '5' in plot:fig5_1.savefig(f'/Volumes/SDrive/yenting_pa_alignment/results/for_yt/{directory}/figure5_1.png')
        if '6' in plot:fig6.savefig(f'/Volumes/SDrive/yenting_pa_alignment/results/for_yt/{directory}/figure6.png')
        if '7' in plot:fig7.savefig(f'/Volumes/SDrive/yenting_pa_alignment/results/for_yt/{directory}/figure7.png')
        if '7' in plot:fig7_1.savefig(f'/Volumes/SDrive/yenting_pa_alignment/results/for_yt/{directory}/figure7_1.png')
        if '8' in plot:fig8.savefig(f'/Volumes/SDrive/yenting_pa_alignment/results/for_yt/{directory}/figure8.png')


    return opterr, optk51

# manual run
for row,kPA_range in enumerate(['1.0','0.5','0.3']):
    optimize(row,kPA_range,binsnum)
plt.show()


def plotting(kPA_source,kPA_range,apply,k=None):
    plateifu = np.array(data.extract('plateifu'))
    if k!=None : apply=[i==k for i in range(9)]
    # find subsamples
    # kPA = data.extract('PA_'+kPA_source+'_'+kPA_range+'RE',tofloat=True)
    # kPAerr = data.extract('PA_'+kPA_source+'_ERR_'+kPA_range+'RE',tofloat=True)
    # kPAdstat = data.extract('PA_'+kPA_source+'_DSTAT_'+kPA_range+'RE',tofloat=True)
    # pn = data.extract('PA_'+kPA_source+'_PXL_'+kPA_range+'RE',tofloat=True)
    global PAdiff,kPA,kPAerr,k51,k51err
    kPA,kPAerr,k51,k51err=kinemetry_result(float(kPA_range))
    PAdiff = findPAdiff(jPA_val,kPA)
    if apply[0]:PAdiff = newPAdiff(kPAerr,0,90,PAdiff)
    if apply[1]:PAdiff = newPAdiff(k51,0,1,PAdiff)

    # plot histogram
    fig = plt.figure(figsize=(18,8))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    datanum = np.count_nonzero(~np.isnan(PAdiff))
    ax1.set_title('%s in %s Re has %d data'%(kPA_source,kPA_range,datanum))
    ax2.set_title('3D angle probability')
    fig.subplots_adjust(bottom=0.3, left=0.25)
    hist, bins = np.histogram([x for x in PAdiff if not np.isnan(x)],bins=binsinit)
    b = ax1.bar(bins[:-1]+45/binsnum,hist, width=90/binsnum)

    # Define axes area
    if apply[0]:kPAerr_slider_ax  = fig.add_axes([0.25, 0.17, 0.6, 0.03])
    if apply[1]:k51_slider_ax  = fig.add_axes([0.25, 0.10, 0.6, 0.03])

    # Sliders
    if apply[0]:kPAerr_slider = Slider(kPAerr_slider_ax, 'kPA error max', 0, 90, valinit=50)
    if apply[1]:k51_slider = Slider(k51_slider_ax, 'k51 max', 0, 1, valinit=0.1)

    def sliders_on_changed(val):
        global PAdiff,hist,bins
        # kPA = data.extract('PA_'+kPA_source+'_'+kPA_range+'RE',tofloat=True)
        # kPAerr = data.extract('PA_'+kPA_source+'_ERR_'+kPA_range+'RE',tofloat=True)
        # kPAdstat = data.extract('PA_'+kPA_source+'_DSTAT_'+kPA_range+'RE',tofloat=True)
        # pn = data.extract('PA_'+kPA_source+'_PXL_'+kPA_range+'RE',tofloat=True)
        # kPA,kPAerr,k51,k51err=kinemetry_result(float(kPA_range))
        PAdiff = findPAdiff(jPA_val,kPA)
        totdatanum = np.count_nonzero(~np.isnan(PAdiff))
        if apply[0]:PAdiff = newPAdiff(kPAerr,0,kPAerr_slider.val,PAdiff)
        if apply[1]:PAdiff = newPAdiff(k51,0,k51_slider.val,PAdiff)
        
        datanum = np.count_nonzero(~np.isnan(PAdiff))
        score = findscore(datanum,totdatanum,kPAerr_slider.val,kPAerr,k51_slider.val,k51)
        ax1.set_title('%s in %s Re has %d data\nscore=%f'%(kPA_source,kPA_range,datanum,score))
        hist, bins = np.histogram([x for x in PAdiff if not np.isnan(x)],bins=binsinit)
        [bar.set_height(hist[i]) for i, bar in enumerate(b)]
        [bar.set_x(bins[i]) for i, bar in enumerate(b)]
        fig.canvas.draw_idle()

    if apply[0]:kPAerr_slider.on_changed(sliders_on_changed)
    if apply[1]:k51_slider.on_changed(sliders_on_changed)

    # Buttons
    post_button_ax = fig.add_axes([0.1, 0.70, 0.1, 0.06])
    post_button = Button(post_button_ax, 'POST', hovercolor='red')
    def post_button_on_clicked(mouse_event):
        global result
        ax2.cla()
        result=analyticmc(PAdiff)
        ax2.set_title('3D angle probability')
        ax2.plot(np.degrees(np.linspace(0.0, np.pi/2, 300)),result)
        # plt.show()
    post_button.on_clicked(post_button_on_clicked)

    export_button_ax = fig.add_axes([0.1, 0.80, 0.1, 0.06])
    export_button = Button(export_button_ax, 'Export', hovercolor='red')
    def export_button_on_clicked(mouse_event):
        PAdiffar = np.array(PAdiff)
        write = csv.writer(open('result.csv', 'w'))
        write.writerow(['What are the rows? 1. plateifus that match the criterion 2. their PA difference, and if calculated, 3. 3d degrees from 0~90 4. corresponding probabilities'])
        write.writerow(plateifu[~np.isnan(PAdiffar)])
        write.writerow(PAdiffar[~np.isnan(PAdiffar)])
        if result is not None:
            write.writerow(np.degrees(np.linspace(0.0, np.pi/2, 300)))
            write.writerow(result)
    export_button.on_clicked(export_button_on_clicked)

    if k!=None:
        criterion_dict = [
            {'kpaerr max':np.nanmedian(np.array(kPAerr)[~np.isnan(np.array(PAdiff))])},
        ]
        a = criterion_dict[k].items()
        for (name,criterion) in criterion_dict[k].items():
            print(name,kPA_source,kPA_range)
            if k==0:kPAerr_slider.set_val(criterion)
                
            post_button_on_clicked(bb.MouseEvent('button_release_event',bb.FigureCanvasBase(fig),650,850))
            PAdiffar = np.array(PAdiff)
            # with open(f'/Volumes/SDrive/yenting/results/for yt/715/{name} for {kPA_source} in {kPA_range} Re.csv', 'w') as f:
            #     write = csv.writer(f)
            #     write.writerow(['What are the rows? 1. plateifus that match the criterion 2. their PA difference, 3. 3d degrees from 0~90 4. corresponding probabilities'])
            #     write.writerow(plateifu[~np.isnan(PAdiffar)])
            #     write.writerow(PAdiffar[~np.isnan(PAdiffar)])
            #     write.writerow(np.degrees(np.linspace(0.0, np.pi/2, 300)))
            #     write.writerow(result)
            # plt.savefig(f'/Volumes/SDrive/yenting/results/for yt/715/{name} for {kPA_source} in {kPA_range} Re.jpg')

    else:plt.show()

# manual run
# plotting(kPA_source,kPA_range,apply)

# auto run
# for k in range(-1,9):
#     for kPA_source in ['STAR','HA','O3']:
#         for kPA_range in ['1','0.5','0.3']:
#             plotting(kPA_source,kPA_range,apply=[i==k for i in range(9)],k=k)

# stellarmass = np.log(data.extract('NSA_ELPETRO_MASS',tofloat=True))
# halomass = data.extract('Mh_L',tofloat=True)
# eps = data.extract('ym_epsilon',tofloat=True)
# lam = data.extract('ym_lambda1',tofloat=True)
# fastrotator = [float(l<(e*0.25+0.1) and e<0.4) if not (np.isnan(e) or np.isnan(l)) else np.nan for (e,l) in zip(eps,lam)]
# objgp = data.extract('objgp(arcsec)',tofloat=True)
# r180 = data.extract('r_200',tofloat=True)
# distance = [o/r if not (np.isnan(o) or np.isnan(r)) else np.nan for (o,r) in zip(objgp,r180)]
# numofgal = data.extract('300kpc 1500 km/s',tofloat=True)
# nearest = data.extract('5th nearest (mpc) by cylindrical method cut = 1000 km/s',tofloat=True)
# bcg = [x if x!=2 else 0.0 for x in data.extract('BCG',tofloat=True)]
# mmg = [x if x!=2 else 0.0 for x in data.extract('MMG',tofloat=True)]
# dstat = data.extract('PA_'+kPA_source+'_DSTAT_'+kPA_range+'RE',tofloat=True)
# pn = data.extract('PA_'+kPA_source+'_PXL_'+kPA_range+'RE',tofloat=True)
# bins_dict = {
#     'sm':[np.nanmin(np.array(stellarmass)[~np.isnan(np.array(PAdiff))]),np.nanmedian(np.array(stellarmass)[~np.isnan(np.array(PAdiff))]),np.nanmax(np.array(stellarmass)[~np.isnan(np.array(PAdiff))])],
#     'hm':[np.nanmin(np.array(halomass)[~np.isnan(np.array(PAdiff))]),np.nanmedian(np.array(halomass)[~np.isnan(np.array(PAdiff))]),np.nanmax(np.array(halomass)[~np.isnan(np.array(PAdiff))])],
#     'dr':[np.nanmin(np.array(distance)[~np.isnan(np.array(PAdiff))]),np.nanmedian(np.array(distance)[~np.isnan(np.array(PAdiff))]),np.nanmax(np.array(distance)[~np.isnan(np.array(PAdiff))])],
#     'nog':[np.nanmin(np.array(numofgal)[~np.isnan(np.array(PAdiff))]),np.nanmedian(np.array(numofgal)[~np.isnan(np.array(PAdiff))]),np.nanmax(np.array(numofgal)[~np.isnan(np.array(PAdiff))])],
#     'nea':[np.nanmin(np.array(nearest)[~np.isnan(np.array(PAdiff))]),np.nanmedian(np.array(nearest)[~np.isnan(np.array(PAdiff))]),np.nanmax(np.array(nearest)[~np.isnan(np.array(PAdiff))])],
#     'pn':[np.nanmin(np.array(pn)[~np.isnan(np.array(PAdiff))]),np.nanmedian(np.array(pn)[~np.isnan(np.array(PAdiff))]),np.nanmax(np.array(pn)[~np.isnan(np.array(PAdiff))])],
# }


#############################################################################
# match yt's data
#############################################################################
# my = match.xify(plateifu,hasjpa)
# additional = filehandling("D:\\yenting\\list.csv")
# yt = additional.extract(-1,['plateifu'])[0]
# print(len(yt),yt)
# both = match.inxy(yt,my)
# inytnotmy = match.inxnoty(yt,my)
# inmynotyt = match.inxnoty(my,yt)
# print(len(both))
# print(plateifu)
# print(match.xify_criterion(plateifu,yt))
# print(inmynotyt)
# print(inmynotyt)
# print(len(my),len(yt),len(both),len(inytnotmy),len(inmynotyt))
# table = Table([match.xify_criterion(plateifu,yt)], names=['yt_selection'])
# ascii.write(table, 'D:\\yenting\\possible galaxies\\result.txt',overwrite=True)


#############################################################################
# measure kPA
#############################################################################
# manga = []
# for plateifu in my:
#     try:
#         output = kPA(plateifu, data,re_criterion_list=[1,0.5,0.3],plot=True,source='STAR',dap='VOR10',para='NSA',measure='KS',binning=False,snthreshold=0)
#         # plt.savefig(f'/Volumes/SDrive/yenting_pa_alignment/results/measurekPA/temp/ifu={plateifu}_manga_result.png')
#         plt.show()
#         plt.close()
#     except Exception:
#         output=[plateifu,0,0,0,0,0,0,0,0,0,0,0,0]
#     manga.append(output)
    # table = Table(np.array(manga).T.tolist(), names=['plateifu','pa1','dstat1','pval1','measurement1','papt5','dstatpt5','pvalpt5','measurementpt5','papt3','dstatpt3','pvalpt3','measurementpt3'])
    # table = Table(np.array(manga).T.tolist(), names=['plateifu','pa1','err1','voffset1','measurement1','dstat1','pval1','papt5','errpt5','voffsetpt5','measurement5','dstatpt5','pvalpt5','papt3','errpt3','voffsetpt3','measurement3','dstatpt3','pvalpt3'])
    # table = Table(np.array(manga).T.tolist(), names=['plateifu','dstat1','pval1','dstatpt5','pvalpt5','dstatpt3','pvalpt3'])
    # ascii.write(table, '/Volumes/SDrive/yenting_pa_alignment/results/measurekPA/temp/result.csv',overwrite=True)

# print(kPA('10216-9101', data, source='stellar',re_criterion_list=[0.5,1,1.5],snthreshold=3))
# plt.show()
# kPA('8095-6101', data, source='stellar',re_criterion_list=[1,1.5,2],snthreshold=0,dap='SPX')
# plt.show()


#############################################################################
# measure jPA
#############################################################################
# resultf, resultv, resultn = [],[],[]
# # for plateifu in data.plateifus[2:]:
# for plateifu in ['11755-3701']:
#     plt.figure(figsize=(18,10))
#     plt.suptitle('final decision : '+data[(plateifu,'my_decision')])
#     perc,deg,top,redshift=data[(plateifu,'perc','float')],45,data[(plateifu,'top','int')],data[(plateifu,'z_1','float')]
#     try:
#         resultf.append(jPA(plateifu, data=data, source='FIRST', perc = perc, deg = deg ,crop=0, redshift=redshift,top=top,ax=plt.subplot(231),measure=False))
#     except Exception:
#         resultf.append([plateifu,0,0,0])
#     try:
#         resultv.append(jPA(plateifu, data=data, source='VLASS', perc = perc, deg = deg ,crop=0, redshift=redshift,top=top,ax=plt.subplot(232),measure=False))
#     except Exception:
#         resultv.append([plateifu,0,0,0])
#     try:
#         resultn.append(jPA(plateifu, data=data, source='NVSS', perc = perc, deg = deg ,crop=0, redshift=redshift,top=top,ax=plt.subplot(233),measure=False))
#         # plt.title('(NVSS)')
#         plt.gca().add_artist(plt.Rectangle((-54,-54), 108, 108,fill=False,color='blue'))
#     except Exception:
#         resultn.append([plateifu,0,0,0])
#     try:
#         plt.subplot(234)
#         first_cutout = mpimg.imread('/Volumes/SDrive/yenting_pa_alignment/radio_and_optical_images/potential/%s.gif'%(plateifu))
#         plt.imshow(first_cutout)
#         plt.title('FIRST cutout')
#         plt.axis('off')
#     except Exception:
#         pass
#     try:
#         plt.subplot(235)
#         optical = mpimg.imread('/Volumes/SDrive/yenting_pa_alignment/radio_and_optical_images/optical2/name=%s_optical.jpg'%(plateifu))
#         plt.imshow(optical)
#         plt.title('SDSS optical image')
#         plt.axis('off')
#     except Exception:
#         pass
#     assert data[(plateifu,'my_decision')] in ['FIRST','VLASS','NVSS','FIRST CATALOG','GIVENPA','DELETE'],'invalid decision'
#     try:
#         if data[(plateifu,'my_decision')] in ['GIVENPA','FIRST CATALOG']:
#             resultn.append(jPA(plateifu, data=data, source='FIRST', perc = perc, deg = deg ,crop=0, redshift=redshift,top=top,ax=plt.subplot(236),measure=False,givenpa=data[(plateifu,'jPA','float')]))
#         elif data[(plateifu,'my_decision')] in ['FIRST','VLASS','NVSS']:
#             resultn.append(jPA(plateifu, data=data, source=data[(plateifu,'my_decision','str')], perc = perc, deg = deg ,crop=0, redshift=data[(plateifu,'z_1','float')],top=top,ax=plt.subplot(236),measure=False))
#         # plt.title('(NVSS)')
#     except Exception:
#         resultn.append([plateifu,0,0,0])
#     # plt.tight_layout()
#     # plt.savefig(f'/Volumes/SDrive/yenting_pa_alignment/results/measurejPA/e16/{plateifu}.png')
#     plt.show()
#     plt.close()
#     # table = Table(np.array(resultf).T.tolist(), names=['plateifu','jpa','err_prop','err_overall'])
#     # ascii.write(table, '/Volumes/SDrive/yenting_pa_alignment/results/measurejPA/e14/result_first.csv',overwrite=True)
#     # table = Table(np.array(resultv).T.tolist(), names=['plateifu','jpa','err_prop','err_overall'])
#     # ascii.write(table, '/Volumes/SDrive/yenting_pa_alignment/results/measurejPA/e14/result_vlass.csv',overwrite=True)
#     # table = Table(np.array(resultn).T.tolist(), names=['plateifu','jpa','err_prop','err_overall'])
#     # ascii.write(table, '/Volumes/SDrive/yenting_pa_alignment/results/measurejPA/e14/result_nvss.csv',overwrite=True)

