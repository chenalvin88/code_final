import numpy as np
import math as math
from numpy import sin as sin
from numpy import cos as cos
from numpy import inner as dot
from math import acos as acos
from math import atan as atan
from mpl_toolkits import mplot3d 
import matplotlib.pyplot as plt
from random import uniform as uniform
from tqdm import tqdm

def discriminant(tp1,tp2):
    u1=np.array([sin(tp1[0])*cos(tp1[1]), sin(tp1[0])*sin(tp1[1]),cos(tp1[0])])
    u2=np.array([sin(tp2[0])*cos(tp2[1]), sin(tp2[0])*sin(tp2[1]),cos(tp2[0])])
    result=acos((dot(u1,u2))/(np.sqrt(dot(u1,u1))*np.sqrt(dot(u2,u2))))
    return result

# note that there are modelingsize^2 modeling points
def to3d(data_list, ax,modelingsize=100,err=2, quiet=True):
    degin3d=list()
    # fig, axes = plt.subplots(2,2, figsize=(15, 10))
    bins = np.linspace(0,np.pi/2,100)
    for data in tqdm([data for data in data_list if not np.isnan(data)]):
        assert data<=90
        for i in (np.arange(0,np.pi,np.pi/modelingsize)):
            for j in np.arange(0,np.pi,np.pi/modelingsize):
                temp = discriminant([i,0],[j,np.radians(data)])
                if temp >np.pi/2 :
                    temp = temp
                    # print('hi')
                    # temp = np.pi-temp
                degin3d.append(np.degrees(temp))
                # if temp>bins[8] and temp<bins[9]:
                #     axes[0][0].scatter(np.degrees(i),np.degrees(j))
                # if temp>bins[9] and temp<bins[10]:
                #     axes[0][1].scatter(np.degrees(i),np.degrees(j))
                # if temp>bins[10] and temp<bins[11]:
                #     axes[1][0].scatter(np.degrees(i),np.degrees(j))
    hist, edges, plot= ax.hist(degin3d, bins=np.degrees(bins), density=True)
    # hist = hist/sum(hist)/np.pi
    # xlist = [(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)]
    # plt.bar(xlist,hist,0.05)
    # print(sum(hist))
    if not quiet:
        # axes[1][1].hist(degin3d, bins=bins, density=True)
        plt.show()

# to3d([30],modelingsize=200, quiet=False)