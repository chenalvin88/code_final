import numpy as np
from astropy.io import fits
from astropy.table import Table, Column, MaskedColumn
from astropy.io import ascii
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import WMAP9 as cosmo
from tqdm.auto import tqdm

#parameters 
czcriterion = 1500 #c*z(km/sec)
n = 5 # nth nearest neighbor
what = 0
numberofgalaxies = 200

#constants
c = 299792.458 #km/sec
H0 = 0.7*100 #km*sec^-1*Mpc^-1
zcriterion = czcriterion/c

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

plateifu, result, result_all, redshift_all = [],[],[],[]
for i in np.arange(what,what+numberofgalaxies,1):
    object = SkyCoord(ra=list_ra[i]*u.degree, dec=list_dec[i]*u.degree, frame='icrs')
    #dist_obj = list_z[i]*c/H0 #Mpc
    #c_obj = SkyCoord(ra=list_ra[i]*u.degree, dec=list_dec[i]*u.degree, distance=dist_obj*u.mpc)
    distance = [99999 for i in range(n+1)]
    redshift = [99999 for i in range(n+1)]

    perc = 0
    for j in tqdm(range(len(ra))):
        surrounding = SkyCoord(ra=ra[j]*u.degree, dec=dec[j]*u.degree, frame='icrs')
        #dist_surr = z[j]*c/H0 #Mpc
        #c_surr = SkyCoord(ra=ra[j]*u.degree, dec=dec[j]*u.degree, distance=dist_surr*u.mpc)
        #if c_obj.separation_3d(c_surr).value<max(distance):distance[np.argmax(distance)]=c_obj.separation_3d(c_surr).value
        if list_ra[i]-ra[j]<5 and list_dec[i]-dec[j]<5 and abs(list_z[i]-z[j])<zcriterion and object.separation(surrounding)<max(distance)*u.arcsec:
            distance[np.argmax(distance)] = (object.separation(surrounding)*cosmo.kpc_proper_per_arcmin(list_z[i])).to(u.Mpc).value
            redshift[np.argmax(distance)] = z[j]*c
    
    plateifu.append(list_plate[i])
    redshift = [x for _,x in sorted(zip(distance,redshift))]
    distance.sort()
    # del redshift[-1]
    # del distance[-1]
    result_all.append(str(distance))
    redshift_all.append(str(redshift))
    # print(c*list_z[i])
    # print(distance, redshift)
    print('5th nearest galaxy of',list_plate[i],'has a 2d distance of',max(distance),'(mpc) and c*z difference',abs(redshift[np.argmax(distance)]-list_z[i]*c),'km/s')
    table = Table([plateifu,result_all,redshift_all], names=['plateifu','distance(mpc)','redshift'])
    ascii.write(table, 'D:\\yenting\\results\\countsurroundings\\5th_nearest_cylindrical_e1\\result1500.txt',overwrite=True)