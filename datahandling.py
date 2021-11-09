# from astropy.io import fits
# from astropy.table import Table, Column, MaskedColumn
# from astropy.io import ascii
import numpy as np
from astropy.io import fits
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import math
import csv
from scipy import interpolate,stats
# from measurejPA_class import jPA

class filehandling:
    def __init__(self, filename):
        self.filename = filename
        f = open(filename,newline='')
        self.rows = csv.DictReader(f)
        self.plateifus = []
        for i,row in enumerate(self.rows):
            if i==0:
                self.header = list(row.keys())
            self.plateifus.append(row['plateifu'])
        # print(self.header)
    def __getitem__(self, x):
        type = None
        if len(x)==2:plateifu,col= x
        if len(x)==3:plateifu,col,type = x
        assert plateifu in self.plateifus, 'invalid plateifu '+plateifu
        assert col in self.header, 'invalid column name '+col
        for row in csv.DictReader(open(self.filename,newline='')):
            if plateifu == row['plateifu']:
                if type=='float':
                    return float(row[col])
                if type=='int':
                    return int(row[col])
                else:
                    return row[col]

    def extract(self,col,selection='all',tofloat=False):
        assert selection in self.header or selection=='all', 'invalid selection'
        assert col in self.header, 'invalid column name'
        result = list()
        for row in csv.DictReader(open(self.filename,newline='')):
            if selection == 'all' and not tofloat:
                result.append(row[col])
            if selection == 'all' and tofloat:
                if row[col]=='':result.append(float('nan'))
                else: result.append(float(row[col]))
            if selection != 'all' and row[selection]=='1' and not tofloat:
                result.append(row[col])
            elif selection != 'all' and row[selection]=='1' and tofloat:
                if row[col]=='':result.append(float('nan'))
                else:result.append(float(row[col]))
        return result
        
# data = filehandling("/Volumes/SDrive/yenting/possible galaxies/500_e2.csv")
# print(data[('10214-3703','index')])    

class match:
    #match by elements
    def inxy(x,y):
        match=[]
        for i in x:
            if i in y:
                match.append(i)
        return match
    def inxnoty(x,y):
        match=[]
        for i in x:
            if i not in y:
                match.append(i)
        return match
    #match by criterion
    def xify(x,y):
        match=[]
        for i in range(len(x)):
            if y[i]==1:
                match.append(x[i])
        return match
    #output criterion
    def xify_criterion(x,y):
        criterion=[]
        for i in x:
            if i in y:criterion.append(1)
            else:criterion.append(0)
        return criterion

# def glue_images(data,index,nvss,first,vlass,crop_n, crop_f, crop_v,savefig=False,table=False):
#     list_first=data.extract('first_data')
#     list_vlass=data.extract('vlass_data')
#     list_ra=data.extract('objra')
#     list_dec=data.extract('objdec')
#     list_firstcatalog=data.extract('POSANG')
#     hdu = fits.open('D:\\yenting\\NVSS_FIRST_SDSS_MaNGACUBE_galaxies\\nvss\\ifu=%s_ra=%4.3f_dec=%4.3f_nvss.fits' %(list_plate,list_ra,list_dec))
#     nvss_image=hdu[0].data
#     if list_first==1:
#         hdu = fits.open('D:\\yenting\\NVSS_FIRST_SDSS_MaNGACUBE_galaxies\\first\\ifu=%s_ra=%4.3f_dec=%4.3f_first.fits' %(list_plate,list_ra,list_dec))
#         first_image=hdu[0].data
#     else:first_image=[[0]]
#     # if list_vlass==1:
#     #     hdu = fits.open('D:\\yenting\\NVSS_FIRST_SDSS_MaNGACUBE_galaxies\\vlass\\ifu=%s_ra=%4.3f_dec=%4.3f_vlass.fits' %(list_plate,list_ra,list_dec))
#     #     vlass_image=hdu[0].data
#     # else:vlass_image=[[0]]
#     opticalimg = mpimg.imread('D:\\yenting\\NVSS_FIRST_SDSS_MaNGACUBE_galaxies\\optical\\ifu=%s_ra=%4.3f_dec=%4.3f_optical_image.jpg'%(list_plate,list_ra,list_dec))
    
#     # if list_first==1:firstmeasurement = mpimg.imread('D:\\yenting\\results\\measurejPA\\all_first\\ifu=%s_ra=%3.2f_dec=%3.2f_jpa.png'%(list_plate,list_ra,list_dec))
#     # else:firstmeasurement=[[0]]
#     # if list_vlass==1:vlassmeasurement = mpimg.imread('D:\\yenting\\results\\measurejPA\\all_vlass\\ifu=%s_ra=%3.2f_dec=%3.2f_jpa.png'%(list_plate,list_ra,list_dec))
#     # else:vlassmeasurement=[[0]]
#     if list_first==1:firstandoptical = mpimg.imread('D:\\yenting\\results\\check_source\\first and optical\\ifu=%s_ra=%4.3f_dec=%4.3f_overlay.jpg'%(list_plate,list_ra,list_dec))
#     else:firstandoptical=[[0]]
#     if list_vlass==1:vlassandoptical = mpimg.imread('D:\\yenting\\results\\check_source\\vlass and optical\\ifu=%s_ra=%4.3f_dec=%4.3f_overlay.jpg'%(list_plate,list_ra,list_dec))
#     else:vlassandoptical=[[0]]
#     #plt.fig()
#     plt.figure(figsize=(20,9.5))
#     plt.subplot(2,4,1)
#     plt.imshow(nvss_image)
#     plt.ylim(0,nvss_image.shape[1])
#     plt.title('nvss (not to scale)')

#     plt.subplot(2,4,2)
#     nvss.append(jPA(index-2, data='nvss', ringsize = 2.5, plott = False, plotf = False, plotsave = False, tablee = False, handmiddle = [], perc = 90, deg = 45, crop=crop_n))

#     if list_first==1:
#         plt.subplot(2,4,3)
#         first.append(jPA(index-2, data='first', ringsize = 2.5, plott = False, plotf = False, plotsave = False, tablee = False, handmiddle = [], perc = 90, deg = 45,crop=crop_f))
        
#     if list_firstcatalog!='' and list_first==1:
#         # list_firstcatalog = 90
#         plt.subplot(2,4,4)
#         m=first_image.shape[0]
#         middle = [m/2-0.5,m/2-0.5]
#         x = np.linspace(-0.5,m-0.5,m)
#         plt.imshow(first_image)
#         plt.colorbar()
#         plt.ylim(-0.5,m-0.5) 
#         plt.plot(x, math.tan(math.radians(list_firstcatalog+90))*(x-middle[0])+middle[1], color='blue', lw=1)
#         plt.title('first catalog, jPA=%4.3f deg'%(list_firstcatalog))
#         # plt.title('measure by hand, jPA=%4.3f deg'%(list_firstcatalog))

#     if list_vlass==1:
#         plt.subplot(2,4,5)
#         vlass.append(jPA(index-2, data='vlass', ringsize = 2.5, plott = False, plotf = False, plotsave = False, tablee = False, handmiddle = [], perc = 90, deg = 45,crop=crop_v))

#     # if list_first==1:
#     #     plt.subplot(2,4,3)
#     #     plt.imshow(firstmeasurement)
#     #     plt.axis('off')
#     #     plt.title('first measurement')

#     # if list_vlass==1:
#     #     plt.subplot(2,4,7)
#     #     plt.imshow(vlassmeasurement)
#     #     plt.axis('off')
#     #     plt.title('vlass measurement')

#     plt.subplot(2,4,6)
#     plt.imshow(opticalimg)
#     plt.axis('off')
#     plt.title('optical')

#     if list_first==1:
#         plt.subplot(2,4,7)
#         plt.imshow(firstandoptical)
#         plt.axis('off')
#         plt.title('first and optical')

#     if list_vlass==1:
#         plt.subplot(2,4,8)
#         plt.imshow(vlassandoptical)
#         plt.axis('off')
#         plt.title('vlass and optical')

#     # plt.tight_layout()
#     if savefig==True:plt.savefig('D:\\yenting\\results\\check_source\\e10_all\\ifu=%s_ra=%4.3f_dec=%4.3f_glue.jpg' %(list_plate,list_ra,list_dec))
#     plt.show()
#     plt.clf()

#     print(nvss)
#     if table == True:
#         table = Table(np.array(nvss).T.tolist(), names=['plateifu','finalPA_avg','finalerr_avg','finalPA_max','finalerr_max'])
#         ascii.write(table, 'D:\\yenting\\results\\check_source\\e10_all\\nvss.txt',overwrite=True)
#         table = Table(np.array(first).T.tolist(), names=['plateifu','finalPA_avg','finalerr_avg','finalPA_max','finalerr_max'])
#         ascii.write(table, 'D:\\yenting\\results\\check_source\\e10_all\\first.txt',overwrite=True)
#         table = Table(np.array(vlass).T.tolist(), names=['plateifu','finalPA_avg','finalerr_avg','finalPA_max','finalerr_max'])
#         ascii.write(table, 'D:\\yenting\\results\\check_source\\e10_all\\vlass.txt',overwrite=True)
    
#     print('finish %s'%(list_plate))

#     return nvss,first,vlass

def comparePA(plateifu,pa1,paerr1,pa2,paerr2,err1,err2):
    plt.figure(figsize = (15,5))
    for k,_ in enumerate(err1):
        # x,y=[],[]
        # for i,_ in enumerate(plateifu):
        #     if paerr1[i]<err1[k] and paerr2[i]<err2[k]:
        #         x.append(pa1[i])
        #         y.append(pa2[i])
        assert len(pa1)==len(pa2)==len(plateifu)
        x = [pa1[i] for i,_ in enumerate(pa1) if paerr1[i]<err1[k] and paerr2[i]<err2[k]]
        y = [pa2[i] for i,_ in enumerate(pa2) if paerr1[i]<err1[k] and paerr2[i]<err2[k]]
        plate = [plateifu[i] for i,_ in enumerate(plateifu) if paerr1[i]<err1[k] and paerr2[i]<err2[k]]

        assert len(x)==len(y)==len(plate)
        anox = [x[i] for i,_ in enumerate(x) if abs(x[i]-y[i])>30 and abs(x[i]-y[i])<150]
        anoy = [y[i] for i,_ in enumerate(y) if abs(x[i]-y[i])>30 and abs(x[i]-y[i])<150]
        anoplate = [plate[i] for i,_ in enumerate(plate) if abs(x[i]-y[i])>40 and abs(x[i]-y[i])<140]

        # print(x,y)
        plt.subplot(1,len(err1),k+1)
        plt.title('errx<%d and erry<%d'%(err1[k],err2[k]))
        plt.plot([0,180],[0,180],color='red')
        plt.scatter(x,y)
        for i,ano in enumerate(anoplate):
            plt.annotate(ano, (anox[i],anoy[i]))

    plt.show()

def comparePAdiff(PAdiff1, PAdiff2, label1, label2, title, parameter, ax, binnum=30, setticks=None, combine=False):
    parameter1 = parameter[~np.isnan(PAdiff1)]
    parameter2 = parameter[~np.isnan(PAdiff2)]
    print(parameter1)
    print(parameter2)
    [dstat,pval] = stats.ks_2samp(parameter1, parameter2)
    if not combine:ax.text(0.05,0.1, r'D$_{KS}$=%.2f'%(dstat), fontsize=13, transform=ax.transAxes, horizontalalignment='left',verticalalignment='top')
    ma=np.nanmax(np.concatenate((parameter1,parameter2)))
    mi=np.nanmin(np.concatenate((parameter1,parameter2)))
    if not combine:
        ax.hist(parameter1,bins=np.linspace(mi,ma,binnum+1),alpha=0.5,label=label1+'\n'+r'$\mu$'+f'={np.nanmean(parameter1):.2f}, m={np.nanmedian(parameter1):.2f}')
        ax.hist(parameter2,bins=np.linspace(mi,ma,binnum+1),alpha=0.5,label=label2+'\n'+r'$\mu$'+f'={np.nanmean(parameter2):.2f}, m={np.nanmedian(parameter2):.2f}')
    else:
        ax.hist(np.concatenate((parameter1,parameter2)),bins=np.linspace(mi,ma,binnum+1),alpha=0.5)

    ax.set_title(title)
    ax.set_ylabel('number of galaxies')
    ax.legend()
    if setticks is not None:
        ax.set_xticks(setticks)

def comparePAdiff_scatter(PAdiff1, PAdiff2, label1, label2, title, parameter1, parameter2, xlabel, ylabel, ax):
    parameter1_label1=parameter1[~np.isnan(PAdiff1)]
    parameter2_label1=parameter2[~np.isnan(PAdiff1)]
    parameter1_label2=parameter1[~np.isnan(PAdiff2)]
    parameter2_label2=parameter2[~np.isnan(PAdiff2)]
    [dstat,pval] = stats.ks_2samp(parameter1_label1, parameter2_label1)
    ax.text(0.05,0.2, r'D$_{KS,x}$=%.2f'%(dstat), fontsize=13, transform=ax.transAxes, horizontalalignment='left',verticalalignment='top')
    [dstat,pval] = stats.ks_2samp(parameter1_label2, parameter2_label2)
    ax.text(0.05,0.1, r'D$_{KS,y}$=%.2f'%(dstat), fontsize=13, transform=ax.transAxes, horizontalalignment='left',verticalalignment='top')
    ax.scatter(parameter1_label1,parameter2_label1,c='r',label=label1)
    ax.scatter(parameter1_label2,parameter2_label2,c='b',label=label2)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend()

def findPAdiff(jPA,kPA):
    PAdiff = jPA.copy()
    for i,_ in enumerate(jPA):
        if isinstance(jPA[i], float) and not math.isnan(jPA[i]) and isinstance(kPA[i], float) and not math.isnan(kPA[i]):
            PAdiff[i] = abs(jPA[i]-(kPA[i]+90.0)%180.0)
            PAdiff[i] = min(PAdiff[i],180-PAdiff[i])
        else:
            PAdiff[i] = float('nan')
    return PAdiff

def newPAdiff(criterion,min,max,PAdiff):
    remove = 0
    for i,_ in enumerate(PAdiff):
        # print(isinstance(criterion[i],float),criterion[i])
        if criterion[i]<min or criterion[i]>max or np.isnan(criterion[i]):
            # print(criterion[i],PAdiff[i])
            if not np.isnan(PAdiff[i]):remove+=1
            PAdiff[i] = np.nan
    datanum = np.count_nonzero(~np.isnan(PAdiff))
    print('removed %d data and %d are left'%(remove, datanum))
    # print(new)
    return PAdiff

def rewrite(list,newelement,criterion,printt=False):
    new = list.copy()
    count = 0
    for i in range(len(new)):
        if str(criterion[i]) == newelement:
            new[i] = float(newelement)
            count+=1
    if printt:print('rewrite %d data'%(count))
    return new

    