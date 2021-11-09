import numpy as np
import math
import sys
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy.table import Table, Column, MaskedColumn
from astropy.io import ascii
from datahandling import filehandling
from plotbin.display_bins_generators import display_bins_generators,display_pixels
from scipy import interpolate,stats

def tran(filename):
    data = open(filename)
    wavelength = []
    transmission = []
    for i in data.readlines():
        wavelength.append(float(i.split()[0]))
        transmission.append(float(i.split()[1]))
    return wavelength,transmission
filename = f"/Users/michaelcheng/Downloads/CFHT_CFHT.U.dat"
wavelength,transmission = tran(filename)
plt.plot(wavelength,transmission,label='U band')
filename = f"/Users/michaelcheng/Downloads/CFHT_CFHT.cfh9402_g'.dat"
wavelength,transmission = tran(filename)
plt.plot(wavelength,transmission,label='g\' band')
filename = f"/Users/michaelcheng/Downloads/CFHT_CFHT.cfh9602_r'.dat"
wavelength,transmission = tran(filename)
plt.plot(wavelength,transmission,label='r\' band')
filename = f"/Users/michaelcheng/Downloads/CFHT_CFHT.cfh9703_i'.dat"
wavelength,transmission = tran(filename)
plt.plot(wavelength,transmission,label='i\' band')
plt.xlim(3000,11000)
plt.xlabel('Wavelength (A)')
plt.ylabel('Transmission')
plt.legend()
plt.tight_layout()
# hdu.info()
# asdf = hdu['DATA_DCBGC'].data[:,20,30]
# plt.plot(asdf)
# print(hdu['PRIMARY'].header.keys)
print('hi')

# data = filehandling("/Volumes/SDrive/yenting_pa_alignment/sauron_and_vla/ATLAS3D.csv")
# name = data.extract('plateifu')

# filename = f"/Volumes/SDrive/yenting_pa_alignment/MaNGA/Pipe3D/manga-10214-3703.Pipe3D.cube.fits.gz"
# hdu = fits.open(filename)
# hdu.info()
# print(hdu['SSP'].header.keys)
# stellar_mass = hdu['SSP'].data[19] #stellar mass density dust corrected in 'm_Sun/spaxels^2'
# plt.imshow(stellar_mass)
# plt.show()
# fits.open(f"/Volumes/SDrive/yenting_pa_alignment/MaNGA/Pipe3D/manga-10214-3703.Pipe3D.cube.fits.gz")['SSP'].data[19] #stellar mass density dust corrected in 'm_Sun/spaxels^2'

# filename = f"/Volumes/SDrive/yenting_pa_alignment/MaNGA/VOR10-MILESHC-MASTARSSP/manga-10512-3701-MAPS-VOR10-MILESHC-MASTARSSP.fits.gz"
# hdu = fits.open(filename)
# hdu.info()
# print(hdu['EMLINE_GFLUX'].header.keys)
# plt.imshow(hdu['SPX_ELLCOO'].data[1],origin='lower')
# X=hdu['SPX_SKYCOO'].data[0]
# Y=hdu['SPX_SKYCOO'].data[1]
# val=hdu['SPX_ELLCOO'].data[1]
# plt.subplot(211)
# display_pixels(X,Y,val,colorbar=True)
# val=hdu['SPX_ELLCOO'].data[0]
# plt.subplot(212)
# display_pixels(X,Y,val,colorbar=True)
# plt.imshow(hdu['SPX_MFLUX'].data*np.sqrt(hdu['SPX_MFLUX_IVAR'].data),origin='lower')
# plt.xticks(np.arange(0,44,2))
# plt.yticks(np.arange(0,44,2))
# plt.scatter(20,20)
# plt.colorbar()
# plt.show()
# spectrum = hdu[0].data
# table = hdu[2].data

# for p in name:
#     found = 0
#     for i in range(9):
#         try:
#             filename = f"/Volumes/SDrive/yenting_pa_alignment/sauron_and_atlas3d/Cappellari2011a_Atlas3D_Paper1_UnBinned_Cubes_v1.0/MS_{p}_r{i}_C2D.fits"
#             hdu = fits.open(filename)
#         except Exception:
#             found+=1
#             # print(filename)
#     # hdu.info()
#     # print(hdu[0].header.keys)
#     ra,dec = hdu[0].header['TCRVL6'],hdu[0].header['TCRVL7']
#     print(p,ra,dec)
# print(hdu[0].data)

# filename = f"/Volumes/SDrive/yenting_pa_alignment/files/dapall-v3_1_1-3.1.0.fits"
# hdu = fits.open(filename)
# hdu.info()
# a = hdu['PRIMARY'].header.keys
# # print(a)
# '''
# 35:
# ELG17   = 'OIII-5008'
# ELG24   = 'Ha-6564 '
# ELG12   = 'Hdel-4102'

# 46:
# SPI22   = 'HDeltaA '
# SPI24   = 'HDeltaF '
# '''
# b = hdu['VOR10-MILESHC-MASTARSSP'].data['SPECINDEX_1RE']
# # b = hdu['VOR10-MILESHC-MASTARSSP'].data.keys
# # spectrum = hdu[0].data
# # table = hdu[2].data
# print('hi')

# x = table["A"] # Coordinates of the original spaxels in arcsec (north is up)
# y = table["D"]
# flux = np.mean(spectrum, 1) # surface brightness

# filename = "/Volumes/SDrive/yenting_pa_alignment/vla_proposal/sauron_and_atlas3d/atlas3d_stellar_kinematics/cappellari2011a/PXF_bin_MS_NGC3665_r1_idl.fits.gz"
# hdu = fits.open(filename)
# table = hdu[1].data

# xgen = table["XS"] # Voronoi generators
# ygen = table["YS"]
# velbin = table["VPXF"] # Mean stellar velocity

# display_bins_generators(xgen, ygen, velbin, x, y)
# plt.tricontour(x, y, -2.5*np.log10(flux/np.max(flux)), levels=np.arange(20)) # 1 mag contours
# ax.tricontour(x, y, -2.5*np.log10(flux/np.max(flux).ravel()),levels=np.arange(20), colors='k')  # 1 mag contours
# plt.show()

# list = ascii.read("D:\\yenting\\codes\\r&absm.txt")
# r = list['z_1']
# absm = list['Mr01']

# post=fits.open('D:\\yenting\\post_catalog.dr72all0.fits')
# post.info()
# print(post[1].header.keys)
# post[1].columns.info()
# ra = post[1].data['ra']
# dec = post[1].data['dec']
# z1 = post[1].data['z']

# kcorrect=fits.open('D:\\yenting\\kcorrect.nearest.model.z0.10.fits')
# kcorrect.info()
# print(kcorrect[1].header.keys)
# kcorrect[1].columns.info()
# z2 = kcorrect[1].data['z']

# match = fits.open('D:\\yenting\\matchofpostandkcorrect.fits')
# match.info()
# print(match[1].header.keys)
# match[1].columns.info()
# z1 = match[1].data['z_1']
# z2 = match[1].data['z_2']
# absmag = [row[3] for row in match[1].data['absm']]
# #print(absmag)
# plt.scatter(z1,absmag,color='blue')
# plt.scatter(r,absm,color='red')
# plt.ylim(-5,-30)
# plt.plot([max(r),max(r)],[-5,-30])
# plt.plot([min(z1),max(z1)],[-20.5,-20.5])
# plt.show()

# x = [2.0115268, 1.5115345, 1.011542, 0.51154953, 0.011557105, -0.48843533, -0.98842776, 2.5115237, 2.0115304, 1.5115371, 1.0115438, 0.5115504, 0.011557105, -0.48843622]
# y = [-4.5015883, -4.501585, -4.501583, -4.5015817, -4.501581, -4.5015817, -4.501583, -4.0015917, -4.001588, -4.0015845, -4.0015826, -4.001581, -4.0015807, -4.001581]
# val = [-24.841093, -26.078999, -18.280184, 5.6414766, 22.37421, 21.009716, -1.6568977, -5.937025, 2.2560012, 7.6078825, 16.814026, 19.093033, 25.717766, 19.902508]

# display_pixels(x, y, val,colorbar=True)
# plt.show()

# xold = np.array([4,5,6,4,5,6,4,5,6])
# yold = np.array([1,1,1,2,2,2,3,3,3])
# xnew = np.array([1,2,3,1,2,3,1,2,3])
# ynew = np.array([1,1,1,2,2,2,3,3,3])
# vel = np.array([1,2,3,4,5,6,7,8,9])
# xyIn = np.column_stack([xold, yold])
# xyOut = np.column_stack([xnew, ynew])
# velOut = interpolate.griddata(xyIn, vel, xyOut, method='linear')
# # display_pixels(xold, yold, vel,colorbar=True)
# display_pixels(xnew, ynew, velOut,colorbar=True)
# plt.show()
