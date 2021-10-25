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

#parameters 
kpccriterion = 300 #kpc
czcriterion = 1500 #c*z(km/sec)

#constants
c = 299792.458 #km/sec
zcriterion = czcriterion/c
#print(zcriterion)

data = filehandling('all_e8.csv')
list_plate = data.extract('plateifu')
list_ra = data.extract('OBJRA',tofloat=True)
list_dec = data.extract('OBJDEC',tofloat=True)
list_z = data.extract('Z',tofloat=True)
list_identified = [e for i,e in enumerate(data.extract('identified_mangaid')) if e!='nan']

vagc=fits.open('/Volumes/SDrive/yenting_pa_alignment/files/post_catalog.dr72all0.fits')
ra = vagc[1].data['RA']
dec = vagc[1].data['DEC']
z = vagc[1].data['z']
#print(ra[0],dec[0],z[0])

def find(i):
    diff = kpccriterion*u.kpc/cosmo.kpc_proper_per_arcmin(list_z[i]) #range of kpccriterion kpc
    print(list_plate[i],diff)
    # object = SkyCoord(ra=list_ra[i]*u.degree, dec=list_dec[i]*u.degree, distance=cosmo.comoving_distance(list_z[i]), frame='icrs')
    object = SkyCoord(ra=list_ra[i]*u.degree, dec=list_dec[i]*u.degree, frame='icrs')
    count = 0
    dist_list = []
    for j in range(len(ra)):
        # surrounding = SkyCoord(ra=ra[j]*u.degree, dec=dec[j]*u.degree, distance=cosmo.comoving_distance(z[j]), frame='icrs')
        surrounding = SkyCoord(ra=ra[j]*u.degree, dec=dec[j]*u.degree, frame='icrs')
        # c1.separation_3d(c2) 
        asdf = diff.to(u.deg).value
        if list_ra[i]-ra[j]<asdf and list_dec[i]-dec[j]<asdf and abs(list_z[i]-z[j])<zcriterion and object.separation(surrounding)<diff:
            a = (object.separation(surrounding)*cosmo.kpc_proper_per_arcmin(list_z[i])).to(u.kpc)
            print(a)
            dist_list.append(a.value)
            count+=1
    print('galaxy',list_plate[i],'has',count,'galaxies in a surrounding of',kpccriterion,'kpc and within c*z =',czcriterion,'km/sec')
    return [list_plate[i],count,str(dist_list)]

# for i in range(2,100):
if __name__ == '__main__':
    plateifu, result, distance = [],[],[]
    for i in range(0,len(list_identified),10):
        a_pool = multiprocessing.Pool()
        out = a_pool.map(find,list_identified[i:i+10])
        out = list(map(list, itertools.zip_longest(*out, fillvalue=None)))
        plateifu.extend(out[0])
        result.extend(out[1])
        distance.extend(out[2])
        table = Table([plateifu,result,distance], names=['plateifu','count','distance(kpc)'])
        ascii.write(table, './bufferforfile/300kpc1500cz.csv',overwrite=True)
        print('hi')