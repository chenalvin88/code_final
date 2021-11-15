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

def comparePAdiff(PAdiff1, PAdiff2, label1, label2, title, parameter, ax, binnum=30, binsize=None, setticks=None, combine=False, std_dev=False):
    parameter1 = parameter[~np.isnan(PAdiff1)]
    parameter2 = parameter[~np.isnan(PAdiff2)]
    [dstat,pval] = stats.ks_2samp(parameter1, parameter2)
    # if not combine:ax.text(0.05,0.1, r'D$_{KS}$=%.2f,pval=%.2f'%(dstat,pval), fontsize=13, transform=ax.transAxes, horizontalalignment='left',verticalalignment='top')
    if not combine:ax.text(0.05,0.1, r'p=%.2f'%(pval), fontsize=13, transform=ax.transAxes, horizontalalignment='left',verticalalignment='top')
    ma=np.nanmax(np.concatenate((parameter1,parameter2)))
    mi=np.nanmin(np.concatenate((parameter1,parameter2)))
    if binsize is not None : binnum=int((ma-mi)/binsize)
    if not combine:
        if not std_dev:
            ax.hist(parameter1,bins=np.linspace(mi,ma,binnum+1),alpha=0.5,label=label1+'\n'+r'$\mu$'+f'={np.nanmean(parameter1):.2f}, m={np.nanmedian(parameter1):.2f}, n={np.count_nonzero(~np.isnan(parameter1)):d}')
            ax.hist(parameter2,bins=np.linspace(mi,ma,binnum+1),alpha=0.5,label=label2+'\n'+r'$\mu$'+f'={np.nanmean(parameter2):.2f}, m={np.nanmedian(parameter2):.2f}, n={np.count_nonzero(~np.isnan(parameter2)):d}')
        if std_dev:
            ax.hist(parameter1,bins=np.linspace(mi,ma,binnum+1),alpha=0.5,label=label1+'\n'+r'$\mu$'+f'={np.nanmean(parameter1):.2f}, m={np.nanmedian(parameter1):.2f}, n={np.count_nonzero(~np.isnan(parameter1)):d}'+r'$, \sigma=$'+f'{np.nanstd(parameter1):.2f}')
            ax.hist(parameter2,bins=np.linspace(mi,ma,binnum+1),alpha=0.5,label=label2+'\n'+r'$\mu$'+f'={np.nanmean(parameter2):.2f}, m={np.nanmedian(parameter2):.2f}, n={np.count_nonzero(~np.isnan(parameter2)):d}'+r'$, \sigma=$'+f'{np.nanstd(parameter2):.2f}')
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
    # [dstat,pval] = stats.ks_2samp(parameter1_label1, parameter2_label1)
    # ax.text(0.05,0.2, r'D$_{KS,x}$=%.2f'%(dstat), fontsize=13, transform=ax.transAxes, horizontalalignment='left',verticalalignment='top')
    # [dstat,pval] = stats.ks_2samp(parameter1_label2, parameter2_label2)
    # ax.text(0.05,0.1, r'D$_{KS,y}$=%.2f'%(dstat), fontsize=13, transform=ax.transAxes, horizontalalignment='left',verticalalignment='top')
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

    