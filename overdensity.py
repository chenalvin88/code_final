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

#parameters 
kpccriterion = 300 #kpc
czcriterion = 1000 #c*z(km/sec)

#constants
c = 299792.458 #km/sec
zcriterion = czcriterion/c
#print(zcriterion)

data = filehandling('nsa_vagc_linRG.csv')
# list_plate = data.extract('plateifu')
# list_ra = data.extract('OBJRA',tofloat=True)
# list_dec = data.extract('OBJDEC',tofloat=True)
list_ra = data.extract('nsa_ra',tofloat=True)
list_dec = data.extract('nsa_dec',tofloat=True)
list_z = data.extract('nsa_z',tofloat=True)
list_mangaid = data.extract('identified_mangaid')
list_identified = [e for i,e in enumerate(data.extract('identified_mangaid')) if e!='nan' and e!='']

vagc=fits.open('/Volumes/SDrive/yenting_pa_alignment/files/post_catalog.dr72all0.fits')
ra = vagc[1].data['RA']
dec = vagc[1].data['DEC']
z = vagc[1].data['z']
#print(ra[0],dec[0],z[0])

def find_old(manga_id):
    diff = kpccriterion*u.kpc/cosmo.kpc_proper_per_arcmin(list_z[i]) #range of kpccriterion kpc
    # print(list_plate[i],diff,cosmo.kpc_proper_per_arcmin(list_z[i]))
    # object = SkyCoord(ra=list_ra[i]*u.degree, dec=list_dec[i]*u.degree, distance=cosmo.comoving_distance(list_z[i]), frame='icrs')
    object = SkyCoord(ra=list_ra[i]*u.degree, dec=list_dec[i]*u.degree, frame='icrs')
    count = 0
    dist_list = []
    for j in range(len(ra)):
        # surrounding = SkyCoord(ra=ra[j]*u.degree, dec=dec[j]*u.degree, distance=cosmo.comoving_distance(z[j]), frame='icrs')
        surrounding = SkyCoord(ra=ra[j]*u.degree, dec=dec[j]*u.degree, frame='icrs')
        # c1.separation_3d(c2) 
        asdf = diff.to(u.deg).value
        if abs(list_ra[i]-ra[j])<asdf and abs(list_dec[i]-dec[j])<asdf and abs(list_z[i]-z[j])<zcriterion and object.separation(surrounding)<diff:
            a = (object.separation(surrounding)*cosmo.kpc_proper_per_arcmin(list_z[i])).to(u.kpc)
            print(a)
            dist_list.append(a.value)
            count+=1
    print('galaxy',list_plate[i],'has',count,'galaxies in a surrounding of',kpccriterion,'kpc and within c*z =',czcriterion,'km/sec')
    return [list_plate[i],count,str(dist_list)]

def find(i,density_method):
    assert density_method in ['fixed','nth','both']
    kpc_per_arcmin = cosmo.kpc_proper_per_arcmin(list_z[i]).value
    diff = kpccriterion/kpc_per_arcmin/60 #criterion from kpc to arcmin to deg
    count = 0
    dist_list = []
    distance = [np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf]
    redshift = [np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf]
    if density_method=='nth': accelerate_deg = 5
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
    return [list_plate[i],count,str(dist_list),str(distance),str(redshift)]

# for i in range(2,100):
if __name__ == '__main__':
    density_method = 'both'
    assert density_method in ['fixed','nth','both']
    plateifu, count, distance, nth_distance, nth_redshift = [],[],[],[],[]
    for id in tqdm(list_identified):
        i = list_mangaid.index(id)
        out = find(i,density_method=density_method)
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
