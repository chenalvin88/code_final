import numpy as np
from astropy.io import fits
from astropy.table import Table, Column, MaskedColumn
from astropy.io import ascii
from montecarlo_e1 import discriminant
from astropy.cosmology import WMAP9 as cosmo
from astropy import units as u
from astropy.coordinates import SkyCoord
from tqdm.auto import tqdm


#parameters 
kpccriterion = 200 #kpc
czcriterion = 1000 #c*z(km/sec)
what = 0
numberofgalaxies = 200

#constants
c = 299792.458 #km/sec
zcriterion = czcriterion/c
#print(zcriterion)

list_objects = ascii.read("D:\\yenting\\codes\\extracted_mpl11_my.txt")
list_plate = list_objects['plateifu']
list_ra = list_objects['objra']
list_dec = list_objects['objdec']
list_z = list_objects['z_1']

vagc=fits.open('D:\\yenting\\post_catalog.dr72all0.fits')
vagc.info()
ra = vagc[1].data['RA']
dec = vagc[1].data['DEC']
z = vagc[1].data['z']
#print(ra[0],dec[0],z[0])

plateifu, result, distance = [],[],[]
for i in np.arange(what,what+numberofgalaxies,1):
    print(list_plate[i])
    diff = kpccriterion*u.kpc/cosmo.kpc_proper_per_arcmin(list_z[i]) #range of kpccriterion kpc
    print(diff)
    # object = SkyCoord(ra=list_ra[i]*u.degree, dec=list_dec[i]*u.degree, distance=cosmo.comoving_distance(list_z[i]), frame='icrs')
    object = SkyCoord(ra=list_ra[i]*u.degree, dec=list_dec[i]*u.degree, frame='icrs')
    count = 0
    temp = []
    for j in tqdm(range(len(ra))):
        # surrounding = SkyCoord(ra=ra[j]*u.degree, dec=dec[j]*u.degree, distance=cosmo.comoving_distance(z[j]), frame='icrs')
        surrounding = SkyCoord(ra=ra[j]*u.degree, dec=dec[j]*u.degree, frame='icrs')
        # c1.separation_3d(c2) 
        
        if list_ra[i]-ra[j]<diff.to(u.deg).value and list_dec[i]-dec[j]<diff.to(u.deg).value and abs(list_z[i]-z[j])<zcriterion and object.separation(surrounding)<diff:
            a = (object.separation(surrounding)*cosmo.kpc_proper_per_arcmin(list_z[i])).to(u.kpc)
            print(a)
            temp.append(a.value)
            count+=1
    print('galaxy',list_plate[i],'has',count,'galaxies in a surrounding of',kpccriterion,'kpc and within c*z =',czcriterion,'km/sec')
    
    plateifu.append(list_plate[i])
    result.append(count)
    distance.append(str(temp))
    table = Table([plateifu,result,distance], names=['plateifu','count','distance(kpc)'])
    ascii.write(table, 'D:\\yenting\\results\\countsurroundings\\fixed aperture e1\\200kpc1000cz1.txt',overwrite=True)