from sympy import Symbol,symbols, diff, integrate, sin, cos, sqrt, acos, Function
from sympy.utilities.lambdify import lambdify, implemented_function
from sympy.abc import x,y,a,phi
from sympy.geometry import Curve
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib.widgets import Slider, Button, RadioButtons
from scipy.integrate import quad
import scipy.special as special
from tqdm import tqdm
import warnings

length=600

def x1(a,phi):
    return acos((2/(cos(phi)-1))*(cos(a)-(cos(phi)+1)/2))/sqrt(2)
    
def x2(a,phi):
    return (np.pi-a)/sqrt(2)
    
def g(x,a,phi):
    return acos((2/(cos(phi)+1))*(cos(a)-cos(sqrt(2)*x)*(cos(phi)-1)/2))/sqrt(2)

def f(x,a,phi):
    return diff(g(x,a,phi),a)
    
def p(newx1,newx2,newf,a,phi):
    warnings.filterwarnings("ignore")
    if a>=phi:
        result,err=quad(newf,0,newx2(a,phi),args=(a,phi))
    elif a<phi:
        result,err=quad(newf,newx1(a,phi),newx2(a,phi),args=(a,phi))
    warnings.filterwarnings("default")
    return result
    
def one_analyticmc(phi_in,debug=False,plot=False):
    xlist = np.linspace(0.0, np.pi, length)
    a_0 = np.radians(30)
    phi_0 = np.radians(phi_in)

    newx1 = lambdify((a,phi), x1(a,phi), modules=['numpy'])
    newx2 = lambdify((a,phi), x2(a,phi), modules=['numpy'])
    newg = lambdify((x,a,phi), g(x,a,phi), modules=['numpy'])
    newf = lambdify((x,a,phi), f(x,a,phi), modules=['numpy'])

    result = np.array([p(newx1,newx2,newf,a,phi_0) for a in xlist])
    result[0],result[-1]=0,0
    nans, w= np.array([a or b for a,b in zip(np.isnan(result),np.isinf(result))]), lambda z: z.nonzero()[0]
    result[nans]= np.interp(w(nans), w(~nans), result[~nans])
    result=result[:len(result)//2]#+np.flip(result[len(result)//2:],0)
    result = np.array([length*e/sum(result)/180 for e in result])
    # print(sum(result))

    if debug==True:
        # plotting
        fig = plt.figure()
        ax = fig.add_subplot(111)
        # Adjust the subplots region to leave some space for the sliders and buttons
        fig.subplots_adjust(bottom=0.25)


        # Draw the initial plot
        # The 'line' variable is used for modifying the line later
        [line1] = ax.plot(xlist, newg(xlist, a_0, phi_0), linewidth=2, color='red')
        [line2] = ax.plot(xlist, newf(xlist, a_0, phi_0), linewidth=2, color='blue')
        [line3] = ax.plot(xlist, np.full_like(xlist, p(newx1,newx2,newf,a_0,phi_0)), linewidth=2, color='blue')
        ax.set_xlim([0, np.pi])
        ax.set_ylim([0, np.pi])

        # Add two sliders for tweaking the parameters

        # Define an axes area and draw a slider in it
        a_slider_ax  = fig.add_axes([0.15, 0.15, 0.65, 0.03])
        a_slider = Slider(a_slider_ax, 'a', 0, np.pi, valinit=a_0)

        # Draw another slider
        phi_slider_ax = fig.add_axes([0.15, 0.1, 0.65, 0.03])
        phi_slider = Slider(phi_slider_ax, 'phi', 0, np.pi, valinit=phi_0)

        # Define an action for modifying the line when any slider's value changes
        def sliders_on_changed(val):
            line1.set_xdata(xlist)
            line1.set_ydata(newg(xlist, a_slider.val, phi_slider.val))
            line2.set_xdata(xlist)
            line2.set_ydata(newf(xlist, a_slider.val, phi_slider.val))
            line3.set_ydata(np.full_like(xlist, p(newx1,newx2,newf,a_slider.val, phi_slider.val)))
            fig.canvas.draw_idle()
        a_slider.on_changed(sliders_on_changed)
        phi_slider.on_changed(sliders_on_changed)

        # Add a button for resetting the parameters
        reset_button_ax = fig.add_axes([0.8, 0.025, 0.1, 0.04])
        reset_button = Button(reset_button_ax, 'Reset', hovercolor='0.975')
        def reset_button_on_clicked(mouse_event):
            phi_slider.reset()
            a_slider.reset()
        reset_button.on_clicked(reset_button_on_clicked)
        plt.show()
    if plot==True:
        plt.plot(np.linspace(0.0, np.pi/2, int(length/2)),result)
        # plt.show()
    return result

# result=one_analyticmc(30,debug=True,plot=True)
# print(result)
# plt.plot(np.degrees(np.linspace(0.0, np.pi, length)),result)
# plt.show()

def analyticmc(data,ax):
    alist = np.linspace(0.0, np.pi/2, int(length/2))
    result = np.zeros_like(alist)
    print('finding analytic result...')
    for i,element in enumerate(tqdm(data)):
        if not np.isnan(element):
            result+=one_analyticmc(element)
    result/=np.count_nonzero(~np.isnan(data))
    ax.plot(np.degrees(alist),result)
    ax.set_yticklabels(ax.get_yticks(),rotation=45)
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
    # plt.show()
    return result

# result=analyticmc([30])
# print(result)
# plt.plot(np.degrees(np.linspace(0.0, np.pi/2, int(length/2))),result)
# plt.show()