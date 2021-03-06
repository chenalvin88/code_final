import numpy as np
from astropy.io import fits
from astropy.table import Table, Column, MaskedColumn
from astropy.io import ascii
from astropy.stats import sigma_clip
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from mpl_toolkits.axes_grid1 import make_axes_locatable
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
import os
from astropy.cosmology import WMAP9 as cosmo
import astropy.units as u

#############################################################################
# initialize
#############################################################################
data = filehandling("100_e9.csv")
data_all = filehandling("all_e8.csv")
data_nsa = filehandling('nsa_linRG.csv')
data_environment_nsa = filehandling('environment_nsaRQs.csv')
plateifu = np.array(data.extract('plateifu'))
my = data.extract('plateifu','hasjpa')
hasjpa = data.extract('hasjpa')
jPA_val = data.extract('jPA',tofloat=True)

#############################################################################
# find Re distribution
#############################################################################

re_arcsec = data.extract('NSA_ELPETRO_TH50_R',tofloat=True)
redshift = data.extract('z_1',tofloat=True)
re_kpc = (cosmo.kpc_proper_per_arcmin(redshift).to('kpc/arcsec')*re_arcsec*u.arcsec).value
# plt.subplot(211).hist(re_arcsec)
# plt.subplot(211).set_title('Re (arcsec)')
# plt.subplot(211).text(0.7,0.9,r'$\mu$='+f'{np.nanmean(re_arcsec):.2f}, m={np.nanmedian(re_arcsec):.2f}',transform=plt.gca().transAxes, horizontalalignment='left',verticalalignment='top')
# plt.subplot(212).hist(re_kpc)
# plt.subplot(212).set_title('Re (kpc)')
# plt.subplot(212).text(0.7,0.9,r'$\mu$='+f'{np.nanmean(re_kpc):.2f}, m={np.nanmedian(re_kpc):.2f}',transform=plt.gca().transAxes, horizontalalignment='left',verticalalignment='top')
# plt.tight_layout()
# plt.show()

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
hand = [40,0.15]

def optimize_hand(kPAerr,k51,k51err,kPAkin,kPAks,pn):
    global opterr,optk51,errlim,k51lim,errmaxlist,k51maxlist,scorelist
    opterr=hand[0]
    optk51=hand[1]
    PAdiff_out = findPAdiff(jPA_val,kPAkin)
    PAdiff_out = newPAdiff(kPAerr,0,opterr,PAdiff_out)
    PAdiff_out = np.array(newPAdiff(k51,0,optk51,PAdiff_out))
    return PAdiff_out

def optimize_score(kPAerr,k51,k51err,kPAkin,kPAks,pn):
    global opterr,optk51,errlim,k51lim,errmaxlist,k51maxlist,scorelist
    errmaxlist,k51maxlist,scorelist=[],[],[]
    errlim=np.nanpercentile(kPAerr,100)
    # k51lim=np.nanpercentile(k51,80)
    k51lim=np.nanmax([e for i,e in enumerate(k51) if e/k51err[i]>3])
    for errmax in np.linspace(0,90,50):
        for k51max in np.linspace(0,1,50):
            PAdiff = findPAdiff(jPA_val,kPAkin)
            totdatanum = np.count_nonzero(~np.isnan(PAdiff))
            PAdiff = newPAdiff(kPAerr,0,errmax,PAdiff)
            PAdiff = newPAdiff(k51,0,k51max,PAdiff)
            PAdiff = newPAdiff(pn,30,np.inf,PAdiff)
            datanum = np.count_nonzero(~np.isnan(PAdiff))
            score = datanum/totdatanum-errmax/errlim-k51max/k51lim-np.nanmean(abs(kPAkin-kPAks))/90
            errmaxlist.append(errmax)
            k51maxlist.append(k51max)
            scorelist.append(score)
    opterr=errmaxlist[np.argmax(scorelist)]
    optk51=k51maxlist[np.argmax(scorelist)]
    PAdiff_out = findPAdiff(jPA_val,kPAkin)
    PAdiff_out = newPAdiff(kPAerr,0,opterr,PAdiff_out)
    PAdiff_out = np.array(newPAdiff(k51,0,optk51,PAdiff_out))
    return PAdiff_out

def optimize_ranking_ellipse(kPAerr,k51,k51err,kPAkin,kPAks,pn):
    global errlim,k51lim,kPAerr_ranking,k51_ranking,rad_ranking
    datanum = 100
    errlim=np.nanpercentile(kPAerr,100)
    k51lim=np.nanpercentile(k51,100)
    # k51lim=np.nanmax([e for i,e in enumerate(k51) if e/k51err[i]>3])
    print('finding better samples by ranking')
    for rad_ranking in (np.arange(0.04,1,1e-6)):
        kPAerr_ranking,k51_ranking = [],[]
        index = []
        for ind,(i,j) in enumerate(zip(kPAerr,k51)):
            if i**2/(rad_ranking*errlim)**2+j**2/(rad_ranking*k51lim)**2<1:
                kPAerr_ranking.append(i)
                k51_ranking.append(j)
                index.append(ind)
        if len(index)==datanum:
            print(f'for {datanum:d} data, found rad_ranking={rad_ranking:.6f}')
            break
        assert len(index)<=datanum, 'not dense enough resolution'
    kPAkin_new = [e if i in index else np.nan for i,e in enumerate(kPAkin)]
    jPA_val_new = [e if i in index else np.nan for i,e in enumerate(jPA_val)]
    PAdiff = np.array(findPAdiff(jPA_val_new,kPAkin_new))
    return PAdiff

def optimize_ranking_max(kPAerr,k51,k51err,kPAkin,kPAks,pn):
    global kPAerr_ranking,k51_ranking
    kPAerr_ranking,k51_ranking = np.copy(kPAerr),np.copy(k51)
    index = []
    datanum = 100
    for i in range(int(datanum/2)):
        index.append(np.nanargmax(kPAerr_ranking))
        kPAerr_ranking[np.nanargmax(kPAerr_ranking)]=np.nan
        index.append(np.nanargmax(k51_ranking))
        k51_ranking[np.nanargmax(k51_ranking)]=np.nan
    kPAkin_new = [e if i in index else np.nan for i,e in enumerate(kPAkin)]
    jPA_val_new = [e if i in index else np.nan for i,e in enumerate(jPA_val)]
    PAdiff = np.array(findPAdiff(jPA_val_new,kPAkin_new))
    return PAdiff

def optimize_ranking_clipping(kPAerr,k51,k51err,kPAkin,kPAks,pn):
    global kPAerr_ranking,k51_ranking
    sigma = 3
    index = np.full_like(kPAerr,True,dtype=bool)
    converge,count = np.inf,0
    while True:
        for i in range(len(kPAerr)):
            if kPAerr[i]>np.nanmean(kPAerr[index])+sigma*np.nanstd(kPAerr[index]):
                index[i]=False
        count+=1
        if converge==np.count_nonzero(index):break
        converge=np.count_nonzero(index)
        for i in range(len(kPAerr)):
            if k51[i]>np.nanmean(k51[index])+sigma*np.nanstd(k51[index]):
                index[i]=False
        count+=1
        if converge==np.count_nonzero(index):break
        converge=np.count_nonzero(index)
    print(f'takes {count} times until convergence of {converge} galaxies')
    kPAerr_ranking = [e if index[i] else np.nan for i,e in enumerate(kPAerr)]
    k51_ranking = [e if index[i] else np.nan for i,e in enumerate(k51)]
    kPAkin_new = [e if index[i] else np.nan for i,e in enumerate(kPAkin)]
    jPA_val_new = [e if index[i] else np.nan for i,e in enumerate(jPA_val)]
    PAdiff = np.array(findPAdiff(jPA_val_new,kPAkin_new))
    return PAdiff

def optimize_ranking_percentage(kPAerr,k51,k51err,kPAkin,kPAks,pn):
    global kPAerr_ranking,k51_ranking
    perc = 90
    index = np.full_like(kPAerr,True,dtype=bool)
    datanum = 100
    while True:
        for i in range(len(kPAerr)):
            if kPAerr[i]>np.nanpercentile(kPAerr[index],perc):
                index[i]=False
        if np.count_nonzero(index)<=datanum:break
        for i in range(len(kPAerr)):
            if k51[i]>np.nanpercentile(k51[index],perc):
                index[i]=False
        if np.count_nonzero(index)<=datanum:break
    kPAerr_ranking = [e if index[i] else np.nan for i,e in enumerate(kPAerr)]
    k51_ranking = [e if index[i] else np.nan for i,e in enumerate(k51)]
    kPAkin_new = [e if index[i] else np.nan for i,e in enumerate(kPAkin)]
    jPA_val_new = [e if index[i] else np.nan for i,e in enumerate(jPA_val)]
    PAdiff = np.array(findPAdiff(jPA_val_new,kPAkin_new))
    return PAdiff

def radius_dependent_parameters_all():
    print('finding radius dependent parameters in all dap catalog')
    stellarmass_all, o3_lum_all, ha_lum_all, dn4000_specindex_all, hd_specindex_all, sfr_all = np.full((3,len(data_all.extract('plateifu'))),np.nan),np.full((3,len(data_all.extract('plateifu'))),np.nan),np.full((3,len(data_all.extract('plateifu'))),np.nan),np.full((3,len(data_all.extract('plateifu'))),np.nan),np.full((3,len(data_all.extract('plateifu'))),np.nan),np.full((3,len(data_all.extract('plateifu'))),np.nan)
    mangaid_all = np.full_like(data_all.extract('plateifu'), np.nan)
    h=0.7
    lumdist = np.array(data_all.extract('LDIST_NSA_Z',tofloat=True))/h
    for j,plateifu in enumerate(tqdm(data_all.extract('plateifu')[0:])):
        if data_all[plateifu,'DAPDONE']=='TRUE' and data_all[plateifu,'DAPQUAL','float']==0 and data_all[plateifu,'MANGAID'] not in mangaid_all: #dapdone=True, dapqual=0, mangaid not repeat
            plate = plateifu.split('-')[0]
            ifu = plateifu.split('-')[1]
            # path_pipe3d = f"./bufferforfile/manga-{plateifu}.Pipe3D.cube.fits.gz"
            # url_pipe3d = f"--user='sdss' --password='2.5-meters' http://data.sdss.org/sas/mangawork/manga/sandbox/pipe3d/v3_1_1/3.1.1/{plate}/manga-{plateifu}.Pipe3D.cube.fits.gz"
            # os.system(f"wget -O {path_pipe3d} {url_pipe3d}")
            path_hyb10 = f"./bufferforfile/manga-{plateifu}-MAPS-HYB10-MILESHC-MASTARSSP.fits.gz"
            url_hyb10 = f"--user='sdss' --password='2.5-meters' https://data.sdss.org/sas/mangawork/manga/spectro/analysis/MPL-11/HYB10-MILESHC-MASTARSSP/{plate}/{ifu}/manga-{plateifu}-MAPS-HYB10-MILESHC-MASTARSSP.fits.gz"
            os.system(f"wget -O {path_hyb10} {url_hyb10}")
            mangaid_all[j]=data_all[plateifu,'MANGAID']
            for i,kPA_range in enumerate(['1.0','0.5','0.3']):
                index = kPA(plateifu, data_all,re_criterion_list=[float(kPA_range)],plot=False, source='STAR',dap='HYB10',para='NSA',measure='None',binning=False,snthreshold=0,mangadir='./bufferforfile')[-1]
                # stellarmass_all[i][j] = np.log10(np.sum(10**fits.open(path_pipe3d)['SSP'].data[19][index]))
                hdu = fits.open(path_hyb10)
                # o3_gflux_map = hdu['EMLINE_GFLUX'].data[16][index]*1e-17 #Oiii-5008 1E-17 erg/s/cm^2/spaxel
                # o3_lum_all[i][j]=np.log10(np.sum(4*np.pi*o3_gflux_map*(lumdist[j]*(3.085678*1e24))**2))
                # dn4000_specindex_all[i][j]=np.median(hdu['SPECINDEX'].data[44][index])
                # hd_specindex_all[i][j]=np.median(hdu['SPECINDEX'].data[21][index])
                ha_gflux_map = hdu['EMLINE_GFLUX'].data[23][index]*1e-17 #Ha-6564 1E-17 erg/s/cm^2/spaxel
                ha_lum_all[i][j]=np.log10(np.sum(4*np.pi*ha_gflux_map*(lumdist[j]*(3.085678*1e24))**2))
                sfr_all[i][j]=np.sum(10**(np.log10(4*np.pi*ha_gflux_map*(lumdist[j]*(3.085678*1e24))**2)-41.27))
            # os.system(f"rm {path_pipe3d}")
            os.system(f"rm {path_hyb10}")
        # table = Table([data_all.extract('plateifu'),mangaid_all,stellarmass_all[0],o3_lum_all[0],dn4000_specindex_all[0],hd_specindex_all[0],sfr_all[0],stellarmass_all[1],o3_lum_all[1],dn4000_specindex_all[1],hd_specindex_all[1],sfr_all[1],stellarmass_all[2],o3_lum_all[2],dn4000_specindex_all[2],hd_specindex_all[2],sfr_all[2]], names=['plateifu','mangaid','stellarmass_1re','o3_lum_1re','dn4000_specindex_1re','hd_specindex_1re','sfr_1re','stellarmass_0.5re','o3_lum_0.5re','dn4000_specindex_0.5re','hd_specindex_0.5re','sfr_0.5re','stellarmass_0.3re','o3_lum_0.3re','dn4000_specindex_0.3re','hd_specindex_0.3re','sfr_0.3re'])
        table = Table([mangaid_all,ha_lum_all[0],ha_lum_all[1],ha_lum_all[2]], names=['mangaid','ha_lum_1re','ha_lum_0.5re','ha_lum_0.3re'])
        ascii.write(table, f'all_control_ha.csv',overwrite=True)
# radius_dependent_parameters_all()
        
def find_control_index(separate_criterion='mass',control_group_size=5,input_data=data_all,plateifu=plateifu,selection='all'):
    control_ind = np.full((len(plateifu),control_group_size),np.nan,dtype=object)
    control = np.full((len(plateifu),control_group_size),np.nan,dtype=object)
    print('finding control group')
    all_plateifu = input_data.extract('plateifu',selection=selection)
    if separate_criterion=='mass':separate_criterion_all = input_data.extract('NSA_ELPETRO_MASS',tofloat=True,selection=selection)
    if separate_criterion=='svd':separate_criterion_all = input_data.extract('STELLAR_SIGMA_1RE',tofloat=True,selection=selection)
    sersic_n_all = input_data.extract('NSA_SERSIC_N',tofloat=True,selection=selection)
    re_all = input_data.extract('NSA_ELPETRO_TH50_R',tofloat=True,selection=selection)
    # mangaid_all = input_data.extract('identified_mangaid')
    mangaid_all = all_plateifu
    for i,e in enumerate((plateifu)):
        tuning = 1
        startover = 0
        count = 0
        lim = [[0.8,1.2],[0.8,1.2],[0.8,1.2]]
        if e!='':
            if separate_criterion=='mass':r = data[(e,'NSA_ELPETRO_MASS','float')]
            if separate_criterion=='svd':r = data[(e,'STELLAR_SIGMA_1RE','float')]
            p = data[(e,'NSA_SERSIC_N','float')]
            q = data[(e,'NSA_ELPETRO_TH50_R','float')]
        while e!='' and not np.isnan(r):
            for j,c in enumerate(all_plateifu):
                if c!='' and mangaid_all[j]!='nan':
                    s = separate_criterion_all[j]
                    t = sersic_n_all[j]
                    u = re_all[j]
                    if s>lim[0][0]*r and s<lim[0][1]*r and u>lim[1][0]*q and u<lim[1][1]*q and t>lim[2][0]*p and t<lim[2][1]*p\
                    and c not in plateifu and j not in control_ind[i]:
                        control_ind[i][count] = j
                        control[i][count] = c
                        count+=1
                if j==len(all_plateifu)-1:
                    startover+=1
                    # lim[0] = [0.9-startover//tuning*0.05,1.1+startover//tuning*0.05]
                    # lim[1] = [0.9-startover*0.05,1.1+startover*0.05]
                    # lim[2] = [0.9-startover*0.05,1.1+startover*0.05]
                if count==control_group_size:break
            if count==control_group_size or startover==tuning:break
        if e!='' and count!=control_group_size:
            if separate_criterion=='mass':print(f'separated by {separate_criterion}, {e}, has {np.log10(r/0.7**2),p,q}')
            if separate_criterion=='svd':print(f'separated by {separate_criterion}, {e}, has {r,p,q}')
            print(f'started over {startover} time, has {count} control sample, lim are {lim[0], lim[1], lim[2]}')
    return control_ind,control


def main(row,kPA_range,binsnum,fixed_range=False,frombuffer=True,plot='23',ana=False,find_dependencies=True,control=True,optimize='scoring',directory='1201'):
    assert optimize in ['hand','scoring','ranking_ellipse','ranking_max','ranking_clipping','ranking_percentage']
    if not frombuffer:
        rad,kPAkin,kPAerr,k51,k51err=np.full_like(plateifu,np.nan,dtype=object),np.full_like(plateifu,np.nan,dtype=object),np.full_like(plateifu,np.nan,dtype=object),np.full_like(plateifu,np.nan,dtype=object),np.full_like(plateifu,np.nan,dtype=object)
        kPAks,pn=np.full_like(plateifu,np.nan,dtype=object),np.full_like(plateifu,np.nan,dtype=object)
        for i,e in enumerate(tqdm(plateifu)):
            if e!='':
                ##########################
                if fixed_range:kPA_range=2*re_kpc[i]
                ##########################
                try:kslist=kPA(e, data,re_criterion_list=[float(kPA_range)],plot=False, source='STAR',dap='HYB10',para='NSA',measure='KS',binning=False,snthreshold=0)
                except Exception : kslist=[np.nan]*6
                try:kinlist=kPA(e, data,re_criterion_list=[float(kPA_range)],plot=False, source='STAR',dap='HYB10',para='NSA',measure='Kinemetry',binning=False,snthreshold=0)
                except Exception : kinlist=[np.nan]*6
                kPAks[i]=kslist[1]
                pn[i]=kslist[5]
                rad[i]=kinlist[1]
                kPAkin[i]=kinlist[2]
                kPAerr[i]=kinlist[3]
                k51[i]=kinlist[4]
                k51err[i]=kinlist[5]
        table = Table([plateifu,rad,kPAkin,kPAerr,k51,k51err,kPAks,pn], names=['plateifu','rad','kPAkin','kPAerr','k51','k51err','kPAks','pn'])
        if fixed_range:ascii.write(table, f'buffer_fixed.csv',overwrite=True)
        else:ascii.write(table, f'buffer{kPA_range}.csv',overwrite=True)

    if frombuffer:
        if fixed_range:bufferdata=filehandling(f'./bufferforfile/kpa/buffer_fixed.csv')
        else:bufferdata=filehandling(f'./bufferforfile/kpa/buffer{kPA_range}.csv')
        rad=bufferdata.extract('rad',tofloat=True)
        kPAkin=bufferdata.extract('kPAkin',tofloat=True)
        kPAerr=bufferdata.extract('kPAerr',tofloat=True)
        k51=bufferdata.extract('k51',tofloat=True)
        k51err=bufferdata.extract('k51err',tofloat=True)
        kPAks=bufferdata.extract('kPAks',tofloat=True)
        pn=bufferdata.extract('pn',tofloat=True)

    if optimize=='hand':PAdiff = optimize_hand(np.array(kPAerr),np.array(k51),np.array(k51err),np.array(kPAkin),np.array(kPAks),np.array(pn))
    if optimize=='scoring':PAdiff = optimize_score(np.array(kPAerr),np.array(k51),np.array(k51err),np.array(kPAkin),np.array(kPAks),np.array(pn))
    if optimize=='ranking_ellipse':PAdiff = optimize_ranking_ellipse(np.array(kPAerr),np.array(k51),np.array(k51err),np.array(kPAkin),np.array(kPAks),np.array(pn))
    if optimize=='ranking_max':PAdiff = optimize_ranking_max(np.array(kPAerr),np.array(k51),np.array(k51err),np.array(kPAkin),np.array(kPAks),np.array(pn))
    if optimize=='ranking_clipping':PAdiff = optimize_ranking_clipping(np.array(kPAerr),np.array(k51),np.array(k51err),np.array(kPAkin),np.array(kPAks),np.array(pn))
    if optimize=='ranking_percentage':PAdiff = optimize_ranking_percentage(np.array(kPAerr),np.array(k51),np.array(k51err),np.array(kPAkin),np.array(kPAks),np.array(pn))

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
        if optimize=='scoring':axes1[1][row].axvline(errlim,c='k')
        axes1[1][row].set_title('kPA err distribution')
        axes1[2][row].hist(k51)
        # axes1[2][row].axvline(np.nanpercentile(k51,95),c='k')
        if optimize=='scoring':axes1[2][row].axvline(k51lim,c='k')
        axes1[2][row].set_title('k5/k1 distribution')
        fig1.tight_layout()

    if '2' in plot:
        global fig2,axes2
        if row==0:fig2,axes2=plt.subplots(3,3,figsize=(4*3,10),sharey='row',sharex='col')
        axes2[0][row].set_title(f'{kPA_range} $R_e$', fontsize=15)
        axes2[1][row].text(0.40, 0.95, 'kPA error', transform=axes2[1][row].transAxes)
        axes2[0][0].text(-0.05,0.45,r'$k_5/k_1$',rotation=90, transform=axes2[0][row].transAxes)
        # axes2[0][row].scatter(kPAerr,k51)
        # for i,txt in enumerate(findPAdiff(jPA_val,kPAkin)):
        #     axes2[0][row].annotate(f'{txt:.0f}', (kPAerr[i], k51[i]))
        # axes2[0][row].tricontour(kPAerr,k51,PAdiff)
        divider = make_axes_locatable(axes2[0][row])
        sc = axes2[0][row].scatter(kPAerr,k51,c=findPAdiff(jPA_val,kPAkin),cmap='Blues',s=15)
        cax = divider.append_axes('right', size='5%', pad=0.)
        fig2.colorbar(sc, cax=cax, orientation='vertical')
        if optimize=='scoring':
            axes2[0][0].set_ylabel('optimizing metric', fontsize=15)
            axes2[0][row].axvline(opterr)
            axes2[0][row].axhline(optk51)
            # axes2[0][0].yaxis.set_label_coords(-0.05,0.5)
            # axes2[0][row].clabel(axes2[0][row].tricontour(errmaxlist,k51maxlist,scorelist,30), inline=True, fontsize=10)
            # axes2[0][row].scatter(opterr,optk51)
            # axes2[0][row].annotate(f'max at ({opterr:.2f},{optk51:.2f})',(opterr+.1,optk51+.03))
        if optimize=='hand':
            axes2[0][0].set_ylabel('set limit by hand', fontsize=15)
            axes2[0][row].axvline(hand[0])
            axes2[0][row].axhline(hand[1])
        if optimize=='ranking_ellipse':
            axes2[0][0].set_ylabel('ranking', fontsize=15)
            axes2[0][row].add_artist(Ellipse((0,0),2*rad_ranking*errlim,2*rad_ranking*k51lim,facecolor='none',edgecolor='red'))   
            axes2[0][row].scatter(kPAerr_ranking,k51_ranking,s=3,c='r',alpha=0.5)
        if optimize in ['ranking_max','ranking_clipping','ranking_percentage']:
            axes2[0][0].set_ylabel('ranking', fontsize=15)
            axes2[0][row].scatter(kPAerr_ranking,k51_ranking,s=3,c='r',alpha=0.5)
        axes2[2][row].text(0.40, 0.95, 'PA difference', transform=axes2[2][row].transAxes)
        axes2[1][0].set_ylabel('number of galaxies', fontsize=15)
        axes2[1][row].text(0.50, 0.85, '%d galaxies'%(datanum), fontsize=15, transform=axes2[1][row].transAxes)
        hist, bins = np.histogram([x for x in PAdiff if not np.isnan(x)],bins=binsinit)
        axes2[1][row].bar(bins[:-1]+45/binsnum,hist, width=90/binsnum,color='green',alpha=0.6)
        axes2[1][row].set_xlim(0,90)
        axes2[2][0].set_ylabel('3D angle probability', fontsize=15)
        if ana:analyticmc(PAdiff,ax=axes2[2][row],rticks=row==2)
        axes2[2][row].set_xlabel('angle',fontsize=13)
        fig2.tight_layout()
        fig2.subplots_adjust(hspace=0,wspace=0)
        cax.yaxis.set_ticks_position('left')

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
    stellar_sigma_1re = np.array(data.extract('STELLAR_SIGMA_1RE',tofloat=True)) # km/s
    ha_gew_1re = np.array([float(e.split(',')[23]) if e!='' else np.nan for e in data.extract('EMLINE_GEW_1RE')]) #######################!!!!!!!!!!!!!!!!!!!!!!!!!
    nsa_sersic_n = np.array([e if e!=-9999.0 else np.nan for e in data.extract('NSA_SERSIC_N',tofloat=True)])
    radio_lum = np.array(data.extract('lum',tofloat=True))
    radio_morphology = np.array([e if e==1 or e==2 else np.nan for e in data.extract('radio_morphology',tofloat=True)])
    # radio_morphology = np.array([e if (e==1 or e==2) and data.extract('conf',tofloat=True)[i]==1 else np.nan for i,e in enumerate(data.extract('radio_morphology',tofloat=True))])
    
    lumdist = np.array(data.extract('LDIST_NSA_Z',tofloat=True))/h # h-1 Mpc
    stellarmass = np.log10(np.array(data.extract('NSA_ELPETRO_MASS',tofloat=True))/h**2) # h-2 solar masses
    o3_gflux_1re = np.array([float(e.split(',')[16]) if e!='' else np.nan for e in data.extract('EMLINE_GFLUX_1RE')])*1e-17 # 10-17 erg/s/cm2
    o3_lum_1re = np.array([np.log10(4*np.pi*f*(d*(3.085678*1e24))**2) for (f,d) in zip(o3_gflux_1re,lumdist)])
    dn4000_specindex_1re = np.array([float(e.split(',')[44]) if e!='' else np.nan for e in data.extract('SPECINDEX_1RE')])
    hd_specindex_1re = np.array([float(e.split(',')[21]) if e!='' else np.nan for e in data.extract('SPECINDEX_1RE')]) #HDeltaA
    sfr_1re = np.array(data.extract('SFR_1RE',tofloat=True))/h**2 # h-2 Msun/yr
    print(f'finding parameters in {kPA_range} Re')
    plateifu_all = data_all.extract('plateifu')
    stellarmass_all = np.array(data_all.extract(f'stellarmass_{kPA_range}re',tofloat=True))
    o3_lum_all = np.array(data_all.extract(f'o3_lum_{kPA_range}re',tofloat=True))
    ha_lum_all = np.array(data_all.extract(f'ha_lum_{kPA_range}re',tofloat=True))
    dn4000_specindex_all = np.array(data_all.extract(f'dn4000_specindex_{kPA_range}re',tofloat=True))
    hd_specindex_all = np.array(data_all.extract(f'hd_specindex_{kPA_range}re',tofloat=True))
    sfr_all = np.array(data_all.extract(f'sfr_{kPA_range}re',tofloat=True))
    numofgal_all = np.array([count-1 if dist0<1 else count for count,dist0 in zip(data_all.extract('300kpc1500cz_count',tofloat=True),data_all.extract('0th_kpc_1500cz',tofloat=True))])
    nearest_all = np.array([dist5 if dist0<1 else dist4 for dist4,dist5,dist0 in zip(data_all.extract('4th_kpc_1500cz',tofloat=True),data_all.extract('5th_kpc_1500cz',tofloat=True),data_all.extract('0th_kpc_1500cz',tofloat=True))])
    nearest_density_all = 5/nearest_all**2
    gema_overdensity_all = np.array(data_all.extract('gema_overdensity',tofloat=True))
    mangaid_all = data_all.extract('identified_mangaid')
    stellarmass_r,o3_lum_r,ha_lum_r,dn4000_specindex_r,hd_specindex_r,sfr_r,numofgal,nearest_density,gema_overdensity = np.full_like(plateifu,np.nan,dtype=float),np.full_like(plateifu,np.nan,dtype=float),np.full_like(plateifu,np.nan,dtype=float),np.full_like(plateifu,np.nan,dtype=float),np.full_like(plateifu,np.nan,dtype=float),np.full_like(plateifu,np.nan,dtype=float),np.full_like(plateifu,np.nan,dtype=float),np.full_like(plateifu,np.nan,dtype=float),np.full_like(plateifu,np.nan,dtype=float)
    for i,e in enumerate(tqdm(plateifu)):
        i_inall = plateifu_all.index(e)
        if e!='' and mangaid_all[i_inall]!='nan':
            stellarmass_r[i]=stellarmass_all[i_inall] # stellar mass density dust corrected in 'm_Sun/spaxels^2'
            o3_lum_r[i]=o3_lum_all[i_inall]
            ha_lum_r[i]=ha_lum_all[i_inall]
            dn4000_specindex_r[i]=dn4000_specindex_all[i_inall] #Dn4000
            hd_specindex_r[i]=hd_specindex_all[i_inall] #HDeltaA
            sfr_r[i]=sfr_all[i_inall]
            numofgal[i]=numofgal_all[i_inall]
            nearest_density[i]=nearest_density_all[i_inall]
            gema_overdensity[i]=gema_overdensity_all[i_inall]
    
    ssfr = np.array([np.log10(sfr/sm) for (sfr,sm) in zip(sfr_1re,data.extract('NSA_ELPETRO_MASS',tofloat=True))])
    surface_mass_density = np.array([np.log10(mass/np.pi/re**2) for (mass,re) in zip(data.extract('NSA_ELPETRO_MASS',tofloat=True),data.extract('NSA_ELPETRO_TH50_R',tofloat=True))])
    blackholemass = np.array([10**(8.13+4.02*np.log10(sig/200))/1e8 for (lum,sig) in zip(o3_lum_1re,stellar_sigma_1re)])
    oer = np.array([np.log10(10**(lum)/(1.28*1e46*bhm)) for (lum,bhm) in zip(o3_lum_1re,blackholemass)])
    ssfr_r = np.array([np.log10(sfr/sm) for (sfr,sm) in zip(sfr_r,10**(np.array(stellarmass_r)))])
    surface_mass_density_r = np.array([np.log10(mass/np.pi/(float(kPA_range)*re)**2) for (mass,re) in zip(10**(np.array(stellarmass_r)),data.extract('NSA_ELPETRO_TH50_R',tofloat=True))])
    blackholemass_r = np.array([10**(8.13+4.02*np.log10(sig/200))/1e8 for (lum,sig) in zip(o3_lum_r,stellar_sigma_1re)])
    oer_r = np.array([np.log10(10**(lum)/(1.28*1e46*bhm)) for (lum,bhm) in zip(o3_lum_r,blackholemass_r)])

    halomass = data.extract('Mh_L',tofloat=True)
    objgp = data.extract('objgp(arcsec)',tofloat=True)
    r180 = data.extract('r_200',tofloat=True)
    distance = [o/r if not (np.isnan(o) or np.isnan(r)) else np.nan for (o,r) in zip(objgp,r180)]
    bcg = [x if x!=2 else 0.0 for x in data.extract('BCG',tofloat=True)]
    mmg = [x if x!=2 else 0.0 for x in data.extract('MMG',tofloat=True)]

    if find_dependencies:
        # manga control group
        stellarmass_control,o3_lum_control,ha_lum_control,dn4000_specindex_control,hd_specindex_control,sfr_control,numofgal_control,nearest_density_control,gema_overdensity_control = np.full_like(plateifu,np.nan,dtype=object),np.full_like(plateifu,np.nan,dtype=object),np.full_like(plateifu,np.nan,dtype='float64'),np.full_like(plateifu,np.nan,dtype=object),np.full_like(plateifu,np.nan,dtype=object),np.full_like(plateifu,np.nan,dtype=object),np.full_like(plateifu,np.nan,dtype='float64'),np.full_like(plateifu,np.nan,dtype='float64'),np.full_like(plateifu,np.nan,dtype='float64')
        control_index,_ = find_control_index(separate_criterion='mass')
        for i,e in enumerate(plateifu):
            control_index_i = np.array(list(filter(lambda a: ~np.isnan(a), control_index[i])))
            if not np.isnan(control_index[i][0]) and len(control_index_i)>2:
                stellarmass_control[i] = np.log10(np.mean(10**np.array(stellarmass_all[control_index_i])))
                o3_lum_control[i] = np.log10(np.mean(10**np.array(o3_lum_all[control_index_i])))
                ha_lum_control[i] = np.log10(np.mean(10**np.array(ha_lum_all[control_index_i])))
                dn4000_specindex_control[i] = np.mean(dn4000_specindex_all[control_index_i])
                hd_specindex_control[i] = np.mean(hd_specindex_all[control_index_i])
                sfr_control[i] = np.mean(sfr_all[control_index_i])
                numofgal_control[i] = np.mean(numofgal_all[control_index_i])
                nearest_density_control[i] = np.median(nearest_density_all[control_index_i])
                gema_overdensity_control[i] = np.median(gema_overdensity_all[control_index_i])
        ssfr_control = np.array([np.log10(sfr/sm) for (sfr,sm) in zip(sfr_control,10**(np.array(stellarmass_control)))])
        surface_mass_density_control = np.array([np.log10(mass/np.pi/(float(kPA_range)*re)**2) for (mass,re) in zip(10**(np.array(stellarmass_control)),data.extract('NSA_ELPETRO_TH50_R',tofloat=True))])
        blackholemass_control = np.array([10**(8.13+4.02*np.log10(sig/200))/1e8 for (lum,sig) in zip(o3_lum_control,stellar_sigma_1re)])
        oer_control = np.array([np.log10(10**(lum)/(1.28*1e46*bhm)) for (lum,bhm) in zip(o3_lum_control,blackholemass_control)])
        
        # nsa control group
        plateifu_nsa = np.array(data_environment_nsa.extract('plateifu'))
        numofgal_nsa = np.array([count-1 if dist0<1 else count for count,dist0 in zip(data_environment_nsa.extract('300kpc1500cz_count',tofloat=True),data_environment_nsa.extract('0th_kpc_1500cz',tofloat=True))])
        nearest_nsa = np.array([dist5 if dist0<1 else dist4 for dist4,dist5,dist0 in zip(data_environment_nsa.extract('4th_kpc_1500cz',tofloat=True),data_environment_nsa.extract('5th_kpc_1500cz',tofloat=True),data_environment_nsa.extract('0th_kpc_1500cz',tofloat=True))])
        nearest_density_nsa = 5/nearest_nsa**2
        numofgal_nsa_control,nearest_density_nsa_control=np.full_like(plateifu,np.nan,dtype='float64'),np.full_like(plateifu,np.nan,dtype='float64')
        control_index,control_element = find_control_index(separate_criterion='mass',control_group_size=50,input_data=data_nsa,selection='selection')
        control_element = control_element.astype('float')
        for i,e in enumerate(plateifu):
            control_i = np.array(list(filter(lambda a: ~np.isnan(a), control_element[i])))
            if not np.isnan(control_element[i][0]) and len(control_i)>2:
                numofgal_nsa_control[i] = np.mean([n for p,n in zip(plateifu_nsa,numofgal_nsa) if float(p) in control_i])
                nearest_density_nsa_control[i] = np.median([d for p,d in zip(plateifu_nsa,nearest_density_nsa) if float(p) in control_i])
    

    def dependencies(name,data,axis):
        if name=='stellar mass':comparePAdiff(PAdiff_copy1, PAdiff_copy2, f'low {name}', f'high {name}', f'{name} from pipe 3d', data, ax=axis[0][0],show_pval=False)
        else:comparePAdiff(PAdiff_copy1, PAdiff_copy2, f'low {name}', f'high {name}', f'{name}', data, ax=axis[0][0],show_pval=False)
        comparePAdiff(PAdiff_copy1, PAdiff_copy2, f'low {name}', f'high {name}', 'sersic index n', nsa_sersic_n, ax=axis[0][1],show_pval=False)
        comparePAdiff(PAdiff_copy1, PAdiff_copy2, f'low {name}', f'high {name}', 'radio_morphology',radio_morphology, ax=axis[2][0],binnum=4, setticks=[],show_pval=False)
        axis[2][0].text(1.15,-0.05, 'FR1', fontsize=13, horizontalalignment='center',verticalalignment='top')
        axis[2][0].text(1.85,-0.05, 'FR2', fontsize=13, horizontalalignment='center',verticalalignment='top')
        unique, counts = np.unique(radio_morphology[(~np.isnan(PAdiff_copy1))], return_counts=True)
        bin1 = dict(zip(unique, counts))
        unique, counts = np.unique(radio_morphology[(~np.isnan(PAdiff_copy2))], return_counts=True)
        bin2 = dict(zip(unique, counts))
        comparePAdiff(PAdiff_copy1, PAdiff_copy2, f'low {name}', f'high {name}', 'log(overdensity from GEMA)',np.log10(gema_overdensity), ax=axis[3][2])
        if control:
            comparePAdiff(PAdiff_copy1, PAdiff_copy2, f'low {name}', f'high {name}', 'log(Ha luminosity) - log(Ha luminosity)_control',ha_lum_r-ha_lum_control, ax=axis[0][2])
            comparePAdiff(PAdiff_copy1, PAdiff_copy2, f'low {name}', f'high {name}', 'log(sSFR)-log(sSFR_control)',ssfr_r-ssfr_control, ax=axis[1][0])
            comparePAdiff(PAdiff_copy1, PAdiff_copy2, f'low {name}', f'high {name}', 'stellar surface mass density / ssmd_control',surface_mass_density_r/surface_mass_density_control, ax=axis[1][1])
            comparePAdiff(PAdiff_copy1, PAdiff_copy2, f'low {name}', f'high {name}', 'OER/OER_control',oer_r/oer_control, ax=axis[1][2])
            comparePAdiff_scatter(PAdiff_copy1, PAdiff_copy2, f'low {name}', f'high {name}', '', dn4000_specindex_r/dn4000_specindex_control,hd_specindex_r-hd_specindex_control, r'D$_n$(4000)/D$_n$(4000)_control', r'HDelta$_A$-HDelta$_A$_control', axis[2][1])
            comparePAdiff_scatter(PAdiff_copy1, PAdiff_copy2, f'low {name}', f'high {name}', '', o3_lum_r-o3_lum_control, radio_lum, 'log(O[III] luminosity) - log(O[III] luminosity)_control', r'log($P_{1.4GHz}$)', axis[2][2])
            comparePAdiff(PAdiff_copy1, PAdiff_copy2, f'low {name}', f'high {name}', '(number of galaxies - nog_control) / nog_control',(numofgal-numofgal_nsa_control)/numofgal_nsa_control, ax=axis[3][0], std_dev=True)
            axis[3][1].set_xlim(-2,20)
            comparePAdiff(PAdiff_copy1, PAdiff_copy2, f'low {name}', f'high {name}', '(surface density - sd_control) / sd_control',(nearest_density-nearest_density_nsa_control)/nearest_density_nsa_control, ax=axis[3][1], binsize=0.8, std_dev=True)
            # axis[3][2].scatter((numofgal-numofgal_control)/numofgal_control,(nearest_density-nearest_density_control)/nearest_density_control)
        if not control:
            comparePAdiff(PAdiff_copy1, PAdiff_copy2, f'low {name}', f'high {name}', 'log(Ha luminosity)',ha_lum_r, ax=axis[0][2])
            comparePAdiff(PAdiff_copy1, PAdiff_copy2, f'low {name}', f'high {name}', 'log(sSFR)',ssfr_r, ax=axis[1][0])
            comparePAdiff(PAdiff_copy1, PAdiff_copy2, f'low {name}', f'high {name}', 'stellar surface mass density',surface_mass_density_r, ax=axis[1][1])
            comparePAdiff(PAdiff_copy1, PAdiff_copy2, f'low {name}', f'high {name}', 'OER',oer_r, ax=axis[1][2])
            comparePAdiff_scatter(PAdiff_copy1, PAdiff_copy2, f'low {name}', f'high {name}', '', dn4000_specindex_r,hd_specindex_r, r'D$_n$(4000)', r'HDelta$_A$', axis[2][1])
            comparePAdiff_scatter(PAdiff_copy1, PAdiff_copy2, f'low {name}', f'high {name}', '', o3_lum_r, radio_lum, 'log(O[III] luminosity)', r'log($P_{1.4GHz}$)', axis[2][2])
            comparePAdiff(PAdiff_copy1, PAdiff_copy2, f'low {name}', f'high {name}', 'number of galaxies',numofgal, ax=axis[3][0])
            comparePAdiff(PAdiff_copy1, PAdiff_copy2, f'low {name}', f'high {name}', 'surface density',nearest_density, ax=axis[3][1])
    
    if '3' in plot:
        global fig3,axes3
        if row==0:fig3,axes3=plt.subplots(4,3,figsize=(4*3,10),sharey='row',sharex='col')
        cut = np.nanmedian([e for i,e in enumerate(stellarmass) if not np.isnan(PAdiff[i])])
        PAdiff_copy1=np.copy(PAdiff)
        PAdiff_copy1 = newPAdiff(stellarmass,0,cut,PAdiff_copy1)
        datanum1 = np.count_nonzero(~np.isnan(PAdiff_copy1))
        hist, bins = np.histogram([x for x in PAdiff_copy1 if not np.isnan(x)],bins=binsinit)
        axes3[0][row].set_title(f'{kPA_range} $R_e$', fontsize=15)
        axes3[0][row].bar(bins[:-1]+45/binsnum,hist, width=90/binsnum,color='green',alpha=0.6)
        # axes3[0][row].set_title(f'low stellar mass\n< {cut:.2f}\n{datanum1} galaxies')
        axes3[0][row].text(0.95,0.95, r'log($M_*$)'+f' < {cut:.2f}\n{datanum1} galaxies', fontsize=13, transform=axes3[0][row].transAxes, horizontalalignment='right',verticalalignment='top')
        axes3[0][0].set_ylabel('number of galaxies',fontsize=13)
        axes3[1][0].set_ylabel('3D angle probability',fontsize=13)
        if ana:analyticmc(PAdiff_copy1,ax=axes3[1][row],rticks=row==2)
        
        PAdiff_copy2=np.copy(PAdiff)
        PAdiff_copy2 = newPAdiff(stellarmass,cut,np.inf,PAdiff_copy2)
        datanum2 = np.count_nonzero(~np.isnan(PAdiff_copy2))
        hist, bins = np.histogram([x for x in PAdiff_copy2 if not np.isnan(x)],bins=binsinit)
        axes3[2][row].bar(bins[:-1]+45/binsnum,hist, width=90/binsnum,color='green',alpha=0.6)
        # axes3[2][row].set_title(f'high stellar mass\n> {cut:.2f}\n{datanum2} galaxies')
        axes3[2][row].text(0.95,0.95, r'log($M_*$)'+f' > {cut:.2f}\n{datanum2} galaxies', fontsize=13, transform=axes3[2][row].transAxes, horizontalalignment='right',verticalalignment='top')
        axes3[2][0].set_ylabel('number of galaxies',fontsize=13)
        axes3[3][row].set_xlabel('angle',fontsize=13)
        axes3[3][0].set_ylabel('3D angle probability',fontsize=13)
        if ana:analyticmc(PAdiff_copy2,ax=axes3[3][row],rticks=row==2)
        fig3.tight_layout()
        fig3.subplots_adjust(hspace=0,wspace=0)

        ## for debugging
        # if row==0:
        #     global fig3_1,axes3_1
        #     fig3_1,axes3_1=plt.subplots(3,3,figsize=(4*3,7))
        # comparePAdiff(PAdiff_copy1, PAdiff_copy2, 'low stellar mass', 'high stellar mass', f'{np.count_nonzero(~np.isnan(PAdiff_copy1))}, {np.count_nonzero(~np.isnan(PAdiff_copy2))}', stellarmass, ax=axes3_1[0][row])
        # comparePAdiff(PAdiff_original, PAdiff_copy2, 'low stellar mass', 'high stellar mass', f'{np.count_nonzero(~np.isnan(PAdiff_original))}, {np.count_nonzero(~np.isnan(PAdiff_copy2))}', stellarmass, ax=axes3_1[1][row])
        # comparePAdiff(PAdiff_copy1, PAdiff_original, 'low stellar mass', 'high stellar mass', f'{np.count_nonzero(~np.isnan(PAdiff_copy1))}, {np.count_nonzero(~np.isnan(PAdiff_original))}', stellarmass, ax=axes3_1[2][row])
        
        if find_dependencies:
            global fig3_1,axes3_1
            fig3_1,axes3_1=[[],[],[]],[[],[],[]]
            fig3_1[row],axes3_1[row]=plt.subplots(4,3,figsize=(5*3,13))
            write = csv.writer(open(f'/Volumes/SDrive/yenting_pa_alignment/results/for_yt/{directory}/stellarmass{kPA_range}.csv', 'w'))
            write.writerow(['low stellar mass','PA difference','stellar mass','high stellar mass','PA difference','stellar mass'])
            for values in zip_longest(*[plateifu[~np.isnan(PAdiff_copy1)],PAdiff_copy1[~np.isnan(PAdiff_copy1)],stellarmass_r[~np.isnan(PAdiff_copy1)],plateifu[~np.isnan(PAdiff_copy2)],PAdiff_copy2[~np.isnan(PAdiff_copy2)]],stellarmass_r[~np.isnan(PAdiff_copy2)]):
                write.writerow(values)
            # name='stellar mass'
            # axis=axes3_1[row]
            dependencies(name='stellar mass',data=stellarmass_r,axis=axes3_1[row])
            fig3_1[row].tight_layout()

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
        if ana:analyticmc(PAdiff_copy,ax=axes4[1][row],rticks=row==2)
        
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
        if ana:analyticmc(PAdiff_copy,ax=axes4[3][row],rticks=row==2)
        fig4.tight_layout()
        fig4.subplots_adjust(hspace=0,wspace=0)

        if find_dependencies:
            global fig4_1,axes4_1
            fig4_1,axes4_1=[[],[],[]],[[],[],[]]
            fig4_1[row],axes4_1[row]=plt.subplots(4,3,figsize=(5*3,13))
            write = csv.writer(open(f'/Volumes/SDrive/yenting_pa_alignment/results/for_yt/{directory}/halomass{kPA_range}.csv', 'w'))
            write.writerow(['low halo mass','PA difference','high halo mass','PA difference'])
            for values in zip_longest(*[plateifu[~np.isnan(PAdiff_copy1)],PAdiff_copy1[~np.isnan(PAdiff_copy1)],plateifu[~np.isnan(PAdiff_copy2)],PAdiff_copy2[~np.isnan(PAdiff_copy2)]]):
                write.writerow(values)
            dependencies(name='halo mass',data=halomass,axis=axes4_1[row])
            fig4_1[row].tight_layout()

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
        axes5[0][row].text(0.95,0.95, r'$\sigma$'+f' < {cut:.1f}km/s\n{datanum} galaxies', fontsize=13, transform=axes5[0][row].transAxes, horizontalalignment='right',verticalalignment='top')
        axes5[0][0].set_ylabel('number of galaxies',fontsize=13)
        axes5[1][0].set_ylabel('3D angle probability',fontsize=13)
        if ana:analyticmc(PAdiff_copy1,ax=axes5[1][row],rticks=row==2)
        
        PAdiff_copy2=np.copy(PAdiff)
        PAdiff_copy2 = newPAdiff(stellar_sigma_1re,cut,np.inf,PAdiff_copy2)
        datanum = np.count_nonzero(~np.isnan(PAdiff_copy2))
        hist, bins = np.histogram([x for x in PAdiff_copy2 if not np.isnan(x)],bins=binsinit)
        axes5[2][row].bar(bins[:-1]+45/binsnum,hist, width=90/binsnum,color='green',alpha=0.6)
        axes5[2][row].text(0.95,0.95, r'$\sigma$'+f' > {cut:.1f}km/s\n{datanum} galaxies', fontsize=13, transform=axes5[2][row].transAxes, horizontalalignment='right',verticalalignment='top')
        axes5[2][0].set_ylabel('number of galaxies',fontsize=13)
        axes5[3][row].set_xlabel('angle',fontsize=13)
        axes5[3][0].set_ylabel('3D angle probability',fontsize=13)
        if ana:analyticmc(PAdiff_copy2,ax=axes5[3][row],rticks=row==2)
        fig5.tight_layout()
        fig5.subplots_adjust(hspace=0,wspace=0)

        if find_dependencies:
            global fig5_1,axes5_1
            fig5_1,axes5_1=[[],[],[]],[[],[],[]]
            fig5_1[row],axes5_1[row]=plt.subplots(4,3,figsize=(5*3,13))
            write = csv.writer(open(f'/Volumes/SDrive/yenting_pa_alignment/results/for_yt/{directory}/stellarveldispersion{kPA_range}.csv', 'w'))
            write.writerow(['low stellar velocity dispersion','PA difference','high stellar velocity dispersion','PA difference'])
            for values in zip_longest(*[plateifu[~np.isnan(PAdiff_copy1)],PAdiff_copy1[~np.isnan(PAdiff_copy1)],plateifu[~np.isnan(PAdiff_copy2)],PAdiff_copy2[~np.isnan(PAdiff_copy2)]]):
                write.writerow(values)
            dependencies(name='stellar velocity dispersion',data=stellar_sigma_1re,axis=axes5_1[row])
            fig5_1[row].tight_layout()

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

        cut = np.nanmedian([e for i,e in enumerate(nearest_density) if not np.isnan(PAdiff[i])])
        PAdiff_copy=np.copy(PAdiff)
        PAdiff_copy = newPAdiff(nearest_density,cut,np.inf,PAdiff_copy)
        datanum = np.count_nonzero(~np.isnan(PAdiff_copy))
        hist, bins = np.histogram([x for x in PAdiff_copy if not np.isnan(x)],bins=binsinit)
        axes6[0][2].set_title(f'{kPA_range} $R_e$', fontsize=15)
        axes6[0][2].bar(bins[:-1]+45/binsnum,hist, width=90/binsnum,color='green',alpha=0.6)
        axes6[0][2].text(0.95,0.95, f'surface density > {cut:.2E}\n{datanum} galaxies', fontsize=13, transform=axes6[0][2].transAxes, horizontalalignment='right',verticalalignment='top')
        if ana:analyticmc(PAdiff_copy,ax=axes6[1][2])

        PAdiff_copy=np.copy(PAdiff)
        PAdiff_copy = newPAdiff(nearest_density,0,cut,PAdiff_copy)
        datanum = np.count_nonzero(~np.isnan(PAdiff_copy))
        hist, bins = np.histogram([x for x in PAdiff_copy if not np.isnan(x)],bins=binsinit)
        axes6[2][2].bar(bins[:-1]+45/binsnum,hist, width=90/binsnum,color='green',alpha=0.6)
        axes6[2][2].text(0.95,0.95, f'surface density < {cut:.2E}\n{datanum} galaxies', fontsize=13, transform=axes6[2][2].transAxes, horizontalalignment='right',verticalalignment='top')
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
        axes8[row][4].scatter(PAdiff,nearest_density)
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
    
    if find_dependencies:
        if '3' in plot:fig3_1[row].savefig(f'/Volumes/SDrive/yenting_pa_alignment/results/for_yt/{directory}/figure3_1_{kPA_range}Re.png')
        if '4' in plot:fig4_1[row].savefig(f'/Volumes/SDrive/yenting_pa_alignment/results/for_yt/{directory}/figure4_1_{kPA_range}Re.png')
        if '5' in plot:fig5_1[row].savefig(f'/Volumes/SDrive/yenting_pa_alignment/results/for_yt/{directory}/figure5_1_{kPA_range}Re.png')
    if row==2:
        if '1' in plot:fig1.savefig(f'/Volumes/SDrive/yenting_pa_alignment/results/for_yt/{directory}/figure1.png')
        if '2' in plot:fig2.savefig(f'/Volumes/SDrive/yenting_pa_alignment/results/for_yt/{directory}/figure2.png')
        if '3' in plot:fig3.savefig(f'/Volumes/SDrive/yenting_pa_alignment/results/for_yt/{directory}/figure3.png')
        if '4' in plot:fig4.savefig(f'/Volumes/SDrive/yenting_pa_alignment/results/for_yt/{directory}/figure4.png')
        if '5' in plot:fig5.savefig(f'/Volumes/SDrive/yenting_pa_alignment/results/for_yt/{directory}/figure5.png')
        if '6' in plot:fig6.savefig(f'/Volumes/SDrive/yenting_pa_alignment/results/for_yt/{directory}/figure6.png')
        if '7' in plot:fig7.savefig(f'/Volumes/SDrive/yenting_pa_alignment/results/for_yt/{directory}/figure7.png')
        if '7' in plot:fig7_1.savefig(f'/Volumes/SDrive/yenting_pa_alignment/results/for_yt/{directory}/figure7_1.png')
        if '8' in plot:fig8.savefig(f'/Volumes/SDrive/yenting_pa_alignment/results/for_yt/{directory}/figure8.png')


# manual run
if __name__ == '__main__':
    for row,kPA_range in enumerate(['1.0','0.5','0.3']):
        main(row,kPA_range,binsnum)
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
#         output = kPA(plateifu, data,re_criterion_list=[1,0.5,0.3],plot=True,source='STAR',dap='HYB10',para='NSA',measure='KS',binning=False,snthreshold=0)
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

