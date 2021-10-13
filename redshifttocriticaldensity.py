# https://docs.astropy.org/en/stable/cosmology/
import numpy as np
from astropy.cosmology import WMAP7   # WMAP 7-year cosmology
from astropy.io import ascii
from astropy.table import Table, Column, MaskedColumn
from astropy import units as u
from astropy.coordinates import SkyCoord

h = 70/100 #dimensionless
mass_sun = 1.9891*10**33 #grams

file = ascii.read("D:\\yenting\\codes\\extracted_mpl11_my.txt")
plateifu = file['plateifu']
Mh_M = file['Mh_M']
z = file['z_1'] #galaxy redshift
objra = file['objra']
objdec = file['objdec']
gp_ra = file['gp_ra']
gp_dec = file['gp_dec']
gal_ra = file['gal_ra']
gal_dec = file['gal_dec']
print(plateifu[0],type(objra[0]),z[0])
result = []
objgp,galgp,objgal = [],[],[]
for i in range(len(plateifu)):
    obj = SkyCoord(ra=objra[i]*u.degree, dec=objdec[i]*u.degree, frame='icrs')
    gp = SkyCoord(ra=gp_ra[i]*u.degree, dec=gp_dec[i]*u.degree, frame='icrs')
    gal = SkyCoord(ra=gal_ra[i]*u.degree, dec=gal_dec[i]*u.degree, frame='icrs')
    objgp.append(obj.separation(gp).to(u.arcsec).value)
    galgp.append(gal.separation(gp).to(u.arcsec).value)
    objgal.append(obj.separation(gal).to(u.arcsec).value)
    
    rho_c = WMAP7.critical_density(z[i]).value
    mass = pow(10,float(Mh_M[i]))*(mass_sun/h)
    r_200 = (mass/((4*np.pi/3)*200*rho_c))**(1/3)*3.24077929*10**(-22) #kpc
    print(r_200)
    result.append(r_200)
print(len(plateifu),len(result),len(objgp),len(galgp),len(objgal))
table = Table([plateifu,result,objgp,galgp,objgal], names=['plateifu','r_200','objgp','galgp','objgal'])
ascii.write(table, 'D:\\yenting\\results\\countsurroundings\\r_200.txt',overwrite=True)
