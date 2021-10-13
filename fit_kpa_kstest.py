import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate,stats
from plotbin.display_pixels import display_pixels

def findbest(ind,steps):
    unitvecx = [np.cos(e/(steps-1)*2*np.pi) for e in ind]
    unitvecy = [np.sin(e/(steps-1)*2*np.pi) for e in ind]
    # print(int(np.degrees(np.arctan2(np.mean(unitvecy),np.mean(unitvecx)))%360))
    return int(np.degrees(np.arctan2(np.mean(unitvecy),np.mean(unitvecx)))%360)

# findbest([0,1,2,3,4,5,6,7,350,351,352,353,354,355,356,357,358,359,360],361)

def fit_kpa_ks(xold,yold,V_map,re,ax1=None,ax2=None,debug=False,interp_size=50,steps=361,binning=False):
    if xold.shape[0]<4:
        output=[np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
    else:
        xold, yold, vel = map(np.asarray, [xold, yold, V_map])
        
        size = max(max(abs(xold)),max(abs(yold)))
        side = np.linspace(-size, size, interp_size)
        xnew, ynew = map(np.ravel, np.meshgrid(side, side))  
        
        xyIn = np.column_stack([xold, yold])
        xyOut = np.column_stack([xnew, ynew])
        velOut = interpolate.griddata(xyIn, vel, xyOut, method='linear')

        stat = []
        for pa in np.radians(np.linspace(0.0001, 180, steps)):
            half1_list = [e for i,e in enumerate(velOut) if xnew[i]*np.sin(pa)+ynew[i]*np.cos(pa)>0 and not np.isnan(e)]
            half2_list = [e for i,e in enumerate(velOut) if xnew[i]*np.sin(pa)+ynew[i]*np.cos(pa)<0 and not np.isnan(e)]
            [dstat,pval] = stats.ks_2samp(half1_list, half2_list) #The null hypothesis is that the two distributions are identical
            bins = [i for i in np.linspace(np.nanmin(vel),np.nanmax(vel),50)]
            counthalf1,_ = np.histogram(half1_list, bins=bins)
            counthalf2,_ = np.histogram(half2_list, bins=bins)
            # [dstat,pval] = stats.ks_2samp(half1_list/np.nansum(half1_list), half2_list/np.nansum(half2_list)) #The null hypothesis is that the two distributions are identical
            # bins = [i for i in np.linspace(np.nanmin(half1_list/np.nansum(half1_list)),np.nanmax(half1_list/np.nansum(half1_list)),50)]
            # counthalf1,_ = np.histogram(half1_list/np.nansum(half1_list), bins=bins)
            # counthalf2,_ = np.histogram(half2_list/np.nansum(half2_list), bins=bins)
            stat.append([pa,dstat,pval,bins,counthalf1,counthalf2])

            if debug:
                plt.figure(figsize=(10,5))
                xplot = [xnew[i] for i,_ in enumerate(velOut) if xnew[i]*np.sin(pa)+ynew[i]*np.cos(pa)>0]
                yplot = [ynew[i] for i,_ in enumerate(velOut) if xnew[i]*np.sin(pa)+ynew[i]*np.cos(pa)>0]
                plt.subplot(121)
                plt.gca().set_title(f'{np.degrees(pa):.2f} deg\ninterpolated {len(xold)} pixels to {np.count_nonzero(~np.isnan(velOut))}')
                plt.gca().plot(side,side/np.tan(pa), 'limegreen', linewidth=3) # PA
                plt.gca().plot(side,-side*np.tan(pa), 'k--', linewidth=3) # Zero-velocity line
                # plt.gca().scatter(0,0,s=100)
                # plt.gca().annotate('center',[0,0])
                # display_pixels(xnew, ynew, velOut, vmin=np.nanmin(velOut),vmax=np.nanmax(velOut))
                display_pixels(xplot, yplot, half1_list, vmin=np.nanmin(velOut),vmax=np.nanmax(velOut))
                plt.gca().invert_xaxis()
                
                plt.subplot(122)
                plt.gca().plot([(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)],counthalf1/counthalf1.sum())
                plt.gca().plot([(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)],counthalf2/counthalf2.sum())
                plt.gca().set_title(f'D statistic:{dstat:.4f},p-value:{pval:.2E}')
                plt.pause(1e-10)
                plt.close()
        best = findbest(np.argwhere(np.array(stat).T[1]==np.amax(np.array(stat).T[1])).ravel(), steps)
        output = [np.degrees(stat[best][0]),stat[best][1],stat[best][2],len(xold),stat[best][3],stat[best][4],stat[best][5]]
        
        # plot
        if ax1 is not None:
            if not binning:ax1.set_title(f'interpolated {len(xold)} pixels to {np.count_nonzero(~np.isnan(velOut))}')
            if binning:ax1.set_title(f'interpolated {len(xold)} bins to {np.count_nonzero(~np.isnan(velOut))} pixels')
            ax1.plot(side,side/np.tan(stat[best][0]), 'limegreen', linewidth=3) # PA
            ax1.plot(side,-side*np.tan(stat[best][0]), 'k--', linewidth=3) # Zero-velocity line
            plt.sca(ax1)
            display_pixels(xnew, ynew, velOut, vmin=np.nanmin(velOut),vmax=np.nanmax(velOut),colorbar=True)
            ax1.invert_xaxis()
        if ax2 is not None:
            cv=15*np.sqrt((len(half1_list)+len(half2_list))/(len(half1_list)*len(half2_list)))
            ax2.axhline(cv)
            ax2.set_title(f'fitting process\ncv={cv:.2f}')
            ax2.set_xlabel('degrees')
            ax2.set_xticks(np.linspace(0,180,7))
            ax2.set_ylabel('D statistic')
            ax2.plot((180/np.pi)*np.array(stat).T[0],np.array(stat).T[1])

    return output
# output structure:
# PA,dstat,pval,countgoodpix,bins,counthalf1,counthalf2
        