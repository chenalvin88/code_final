import numpy as np
from astropy.io import fits
from astropy.table import Table, Column, MaskedColumn
from astropy.io import ascii
from datahandling import filehandling
from astropy.cosmology import WMAP9 as cosmo
from astropy import units as u
from astropy.coordinates import SkyCoord
from tqdm.auto import tqdm
import multiprocessing
import itertools
import matplotlib.pyplot as plt
from run import find_control_index

#parameters 
kpccriterion = 300 #kpc
czcriterion = 1500 #c*z(km/sec)

#constants
c = 299792.458 #km/sec
zcriterion = czcriterion/c
#print(zcriterion)


tested = 'nsa' # or 'manga'
if tested=='manga':
    data = filehandling('all_e8.csv')
    list_ra = data.extract('OBJRA',tofloat=True)
    list_dec = data.extract('OBJDEC',tofloat=True)
    list_z = data.extract('Z',tofloat=True)
    list_mangaid = data.extract('identified_mangaid')
    list_identified = [e for i,e in enumerate(list_mangaid) if e!='nan' and e!='']
if tested == 'nsa':
    data = filehandling('nsa_linRG.csv')
    control_ind,control = find_control_index(separate_criterion='mass',control_group_size=50,input_data=data,selection='selection')
    control = control.astype('float')
    control = sorted(list(np.ravel(control[~np.isnan(control)])))
    control = np.array(list(dict.fromkeys(control)),dtype=int).astype(str)
    list_ra = data.extract('NSA_RA',tofloat=True)
    list_dec = data.extract('NSA_DEC',tofloat=True)
    list_z = data.extract('NSA_Z',tofloat=True)
    list_mangaid = data.extract('plateifu')
    list_identified = [e for i,e in enumerate(list_mangaid) if e!='nan' and e!='' and e in control]

vagc=fits.open('/Volumes/SDrive/yenting_pa_alignment/files/post_catalog.dr72all0.fits')
ra = vagc[1].data['RA']
dec = vagc[1].data['DEC']
z = vagc[1].data['z']
#print(ra[0],dec[0],z[0])

def find(id,density_method):
    assert density_method in ['fixed','nth','both']
    i = list_mangaid.index(id)
    kpc_per_arcmin = cosmo.kpc_proper_per_arcmin(list_z[i]).value
    diff = kpccriterion/kpc_per_arcmin/60 #criterion from kpc to arcmin to deg
    count = 0
    dist_list = []
    distance = [np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf]
    redshift = [np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf]
    if density_method in ['nth','both']: accelerate_deg = 5
    else : accelerate_deg = diff
    for j in (range(len(ra))):
        # if abs(list_ra[i]-ra[j])<accelerate_deg and abs(list_dec[i]-dec[j])<accelerate_deg and abs(list_z[i]-z[j])<zcriterion:
        if abs(list_ra[i]-ra[j])<accelerate_deg and abs(list_dec[i]-dec[j])<accelerate_deg and abs(list_z[i]-z[j])<zcriterion:
            separation = np.degrees(np.arccos(np.sin(np.radians(list_dec[i]))*np.sin(np.radians(dec[j]))+np.cos(np.radians(list_dec[i]))*np.cos(np.radians(dec[j]))*np.cos(abs(np.radians(list_ra[i]-ra[j])))))
            a = separation*kpc_per_arcmin*60
            if density_method in ['nth','both'] and a<max(distance):
                distance[np.argmax(distance)] = a
                redshift[np.argmax(distance)] = (z[j]-list_z[i])*c
            if density_method in ['fixed','both'] and separation<diff:
                # print(a)
                dist_list.append(a)
                count+=1
    distance = sorted(distance)
    redshift = [x for _, x in sorted(zip(distance, redshift))]
    # print('galaxy',list_plate[i],'has',count,'galaxies in a surrounding of',kpccriterion,'kpc and within c*z =',czcriterion,'km/sec')
    return [id,count,str(dist_list),str(distance),str(redshift)]

# for i in range(2,100):
if __name__ == '__main__':
    density_method = 'both'
    assert density_method in ['fixed','nth','both']
    plateifu, count, distance, nth_distance, nth_redshift = [],[],[],[],[]
    for id in tqdm(list_identified):
        out = find(id,density_method=density_method)
        plateifu.append(out[0])
        count.append(out[1])
        distance.append(out[2])
        nth_distance.append(out[3])
        nth_redshift.append(out[4])
        if density_method == 'both':
            table = Table([plateifu,count,distance,nth_distance,nth_redshift], names=['plateifu','count in aperture','distance(kpc)','nth distance','nth redsuift difference'])
            ascii.write(table, f'./bufferforfile/{kpccriterion}kpc{czcriterion}cz_6nearest_nsaRQs.csv',overwrite=True)
        if density_method == 'fixed':
            table = Table([plateifu,count,distance], names=['plateifu','count in aperture','distance(kpc)'])
            ascii.write(table, f'./bufferforfile/{kpccriterion}kpc{czcriterion}cz.csv',overwrite=True)
        if density_method == 'nth':
            table = Table([plateifu,nth_distance,nth_redshift], names=['plateifu','nth distance','nth redsuift difference'])
            ascii.write(table, f'./bufferforfile/{czcriterion}cz_6nearest.csv',overwrite=True)
