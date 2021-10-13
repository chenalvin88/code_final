import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.colorbar as cb
import math
from astropy.io import fits
from astropy.table import Table, Column, MaskedColumn
from astropy.io import ascii
import statistics as stats
from datahandling import filehandling,match,comparePA,newPAdiff,rewrite
from astropy import units as u
from astropy.coordinates import SkyCoord
from plotbin.display_pixels import display_pixels
from plotbin.display_bins import display_bins
from plotbin.display_bins_generators import display_bins_generators
from scipy import interpolate,stats
from scipy.spatial import distance
from astropy.cosmology import WMAP9 as cosmo

def interval(number,min,max):
    return (number>=min) & (number<=max)

def goodring(val,pa,x,y,perc,deg,debug=False):
    pamax = pa[np.argmax(val)]
    index = ((val>np.percentile(val,perc))&~(interval(pa,pamax-deg,pamax+deg))&~(interval(pa,pamax+360-deg,pamax+360+deg))&~(interval(pa,pamax+180-deg,pamax+180+deg))&~(interval(pa,pamax-180-deg,pamax-180+deg))&~(interval(pa,pamax-360-deg,pamax-360+deg))).nonzero()
    if debug:
        plt.scatter(x[index],y[index])
        display_pixels(x, y, val, cmap='viridis', pixelsize = np.min(distance.pdist(np.column_stack([x, y]))))
        plt.pause(0.01)
        # plt.show()
        plt.close()
    if len(index[0])==0:print('This is a good ring.')
    else:print('This is a bad ring.')
    return len(index[0])==0

def fold(val,pa,x,y,debug=False):
    val-=np.nanmin(val)
    maxind = np.argmax(val)
    pamax = np.radians(pa[np.argmax(val)]+90)
    if x[maxind]*np.cos(pamax)+y[maxind]*np.sin(pamax)>=0:
        half_move = (x*np.cos(pamax)+y*np.sin(pamax)<0).nonzero()
        half_stay = (x*np.cos(pamax)+y*np.sin(pamax)>=0).nonzero()
    if x[maxind]*np.cos(pamax)+y[maxind]*np.sin(pamax)<=0:
        half_move = (x*np.cos(pamax)+y*np.sin(pamax)>0).nonzero()
        half_stay = (x*np.cos(pamax)+y*np.sin(pamax)<=0).nonzero()
    temp = -2*(x[maxind]*x[half_move]+y[maxind]*y[half_move])/(x[maxind]**2+y[maxind]**2)
    newx = temp*x[maxind]+x[half_move]
    newy = temp*y[maxind]+y[half_move]
    xyIn = np.column_stack([newx, newy])
    xyOut = np.column_stack([x[half_stay], y[half_stay]])
    valInterp = interpolate.griddata(xyIn, val[half_move], xyOut, method='linear') 
    valPlot = np.where(np.isnan(valInterp),0,valInterp)
    valOut = valPlot + val[half_stay]
    padiff=abs(pa[half_stay] - pa[np.argmax(val)])
    padiff=np.where(padiff>180,360-padiff,padiff)
    err = np.sqrt(np.average((padiff)**2,weights = val[half_stay]))
    if debug:
        plt.subplot(221)
        display_pixels(x[half_move], y[half_move], val[half_move],colorbar=True,cmap='viridis',pixelsize = np.min(distance.pdist(np.column_stack([x[half_move], y[half_move]]))))
        # plt.scatter([x[half_move[0][5]]],[y[half_move[0][5]]])
        plt.title('to be folded')
        plt.subplot(222)
        display_pixels(x[half_stay], y[half_stay],valInterp,colorbar=True,cmap='viridis',vmin=0,vmax=np.nanmax(valInterp),pixelsize = np.min(distance.pdist(np.column_stack([x[half_stay], y[half_stay]]))))
        # plt.scatter(newx[5],newy[5])
        plt.title('folded')
        plt.subplot(223)
        display_pixels(x[half_stay], y[half_stay], val[half_stay],colorbar=True,cmap='viridis',pixelsize = np.min(distance.pdist(np.column_stack([x[half_stay], y[half_stay]]))))
        plt.scatter(x[maxind],y[maxind])
        plt.title('half with max')
        plt.subplot(224)
        display_pixels(x[half_stay], y[half_stay], valOut,colorbar=True,cmap='viridis',pixelsize = np.min(distance.pdist(np.column_stack([x[half_stay], y[half_stay]]))))
        plt.title('output')
        plt.tight_layout()
        plt.pause(0.1)
        plt.close()
    print('PA max =',pa[np.argmax(val)],'error =',err)
    return [pa[np.argmax(val)],err]


#top (100-perc) percent bright pixels in deg degrees >> a crucial number!
def jPA(plateifu, data, source, perc = 90, deg = 45 ,crop=0, redshift=None, top=0, givenpa=None, ax=None, measure=True): 
    arcsecpp={
        'FIRST':1.8,
        'VLASS':1.0,
        'NVSS':15
    }
    
    assert source in ['FIRST','VLASS','NVSS']
    if source=='FIRST':assert data[(plateifu,'first_data','float')]==1
    if source=='VLASS':assert data[(plateifu,'vlass_data','float')]==1
    print('starting ',plateifu)
    result,goodring_rad=[],[]
    
    # open source
    if source == 'FIRST' :hdu2 = fits.open('/Volumes/SDrive/yenting_pa_alignment/radio_and_optical_images/first/ifu=%s_ra=%4.3f_dec=%4.3f_first.fits' %(plateifu,data[(plateifu,'objra','float')],data[(plateifu,'objdec','float')]))
    if source == 'VLASS' :hdu2 = fits.open('/Volumes/SDrive/yenting_pa_alignment/radio_and_optical_images/vlass/ifu=%s_ra=%4.3f_dec=%4.3f_vlass.fits' %(plateifu,data[(plateifu,'objra','float')],data[(plateifu,'objdec','float')]))
    if source == 'NVSS' :hdu2 = fits.open('/Volumes/SDrive/yenting_pa_alignment/radio_and_optical_images/nvss/ifu=%s_ra=%4.3f_dec=%4.3f_nvss.fits' %(plateifu,data[(plateifu,'objra','float')],data[(plateifu,'objdec','float')]))

    # image post processing 
    if crop>0 : original_image=(hdu2[0].data)[crop:-crop+1,crop:-crop+1]
    else : original_image=hdu2[0].data
    assert original_image.shape[0] == original_image.shape[1]
    flat_image = original_image.ravel()-np.nanmin(original_image)+1e-5
    m=original_image.shape[0]
    side = np.arange(-m/2+0.5,m/2+0.5)
    X, Y = map(np.ravel, np.meshgrid(side, side)) 

    if source in ['FIRST','VLASS']:
        bottom = 1
        if top==0:top = int(40/arcsecpp[source])
        ringsize = 2.5
    if source in ['NVSS']:
        bottom = 1
        if top==0:top = int(160/arcsecpp[source])
        ringsize = 2.5

    # measure
    for rad in np.arange(bottom,top,ringsize/2):
        print(f'from {rad} to {rad+ringsize}')
        size = rad+ringsize
        side = np.linspace(-size,size,int(10*size))
        XX, YY = map(np.ravel, np.meshgrid(side, side)) 
        xyIn = np.column_stack([X, Y])
        xyOut = np.column_stack([XX, YY])
        interp_image = interpolate.griddata(xyIn, flat_image, xyOut, method='linear') 
        r = np.sqrt(XX**2+YY**2)
        index = ((r>=rad)&(r<rad+ringsize)).nonzero()
        PA = (np.degrees(np.arctan2(YY,XX))-90)
        PA = np.where(PA<-180,PA+360,PA)
        if goodring(interp_image[index],PA[index],XX[index],YY[index],perc,deg):
            result.append(fold(interp_image[index],PA[index],XX[index],YY[index]))
            goodring_rad.append(rad)

    # show results
    if ax==None:ax=plt.gca()
    print('result =',result)
    if redshift is not None:
        kpcparcsec = cosmo.kpc_proper_per_arcmin(redshift).to(u.kpc/u.arcsec).value
        ax.plot([0.6*m/2*arcsecpp[source]-20/kpcparcsec,0.6*m/2*arcsecpp[source]],[-0.8*m/2*arcsecpp[source],-0.8*m/2*arcsecpp[source]],color='k')
        ax.annotate('20 kpc', [0.5*m/2*arcsecpp[source]-20/kpcparcsec,-0.77*m/2*arcsecpp[source]],color='k',size=15)
    
    if givenpa is not None:
        x = np.linspace(-m,m,m+1)
        ax.plot(x, math.tan(math.radians(givenpa+90))*x,'k--', lw=1)
        im=display_pixels(X*arcsecpp[source], Y*arcsecpp[source], flat_image,colorbar=False,cmap='jet',norm=colors.LogNorm(vmin=np.nanmin(flat_image),vmax=np.nanmax(flat_image)))
        plt.colorbar(im)
        finalPA,finalerr,finalerr1=givenpa,0,0
        ax.set_title(f'jPA={givenpa:.2f}$^\circ$\n(given pa on {source})')

    elif len(result)==0:
        finalPA='nan'
        finalerr='nan'
        print('all bad rings')

        im=display_pixels(X*arcsecpp[source], Y*arcsecpp[source], flat_image,colorbar=False,cmap='jet',norm=colors.LogNorm(vmin=np.nanmin(flat_image),vmax=np.nanmax(flat_image)))
        plt.colorbar(im)
        ax.set_title('detected as not a jet')
            
    else:
        average = np.degrees(math.atan2(np.average(np.sin(np.radians([col[0]*2 for col in result])),weights = [col[1]**(-2) for col in result]),np.average(np.cos(np.radians([col[0]*2 for col in result])),weights = [col[1]**(-2) for col in result])))
        if average<0:average = (average+360)/2
        elif average>=0:average = average/2
        finalPA = np.round(average,decimals = 2)
        finalerr = np.round(np.sqrt(1/sum([col[1]**(-2) for col in result])),decimals = 2)
        padiff=abs(np.array([col[0] for col in result]) - finalPA)
        padiff=np.where(padiff>180,360-padiff,padiff)
        padiff=np.where(padiff>90,180-padiff,padiff)
        finalerr1 = np.sqrt(np.average(padiff**2))
        print('PA_max = ',finalPA,', propagated err = ',finalerr, ', overall err = ',finalerr1)

        x = np.linspace(-m*arcsecpp[source],m*arcsecpp[source],2)
        ax.plot(x, math.tan(math.radians(finalPA+90))*x,'k--', lw=1)
        if measure:
            for (e,p) in zip(goodring_rad,result):
                ax.add_artist(plt.Circle((0,0),e*arcsecpp[source],fill=False, lw=0.5))
                p[0] += 90
                ax.scatter(e*arcsecpp[source]*np.cos(np.radians(p[0])),e*arcsecpp[source]*np.sin(np.radians(p[0])),c='white',lw=0,alpha=0.5)
        im=display_pixels(X*arcsecpp[source], Y*arcsecpp[source], flat_image,colorbar=False,cmap='jet',norm=colors.LogNorm(vmin=np.nanmin(flat_image),vmax=np.nanmax(flat_image)),pixelsize = np.min(distance.pdist(np.column_stack([X*arcsecpp[source], Y*arcsecpp[source]]))))
        plt.colorbar(im)
        # ax.set_title(f'jPA={finalPA:.2f}$^\circ$({source})\npropagated err={finalerr:.2f}$^\circ$, overall err={finalerr1:.2f}$^\circ$')
        ax.set_title(f'jPA={finalPA:.2f}$^\circ\pm${finalerr:.2f}$^\circ$({source})')
        
    # plt.show()
    
    return [plateifu,finalPA,finalerr,finalerr1]

# data = filehandling("500_e6.csv")
# for plateifu in ['9093-12703']:
#     jPA(plateifu, data=data, source='FIRST', perc = 94, deg = 45 ,crop=0,redshift=0.02,givenpa=163,measure=True)
#     # plt.savefig('')
#     plt.show()
