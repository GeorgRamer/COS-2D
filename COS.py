# -*- coding: utf-8 -*-
"""
Created on Thu May 22 16:34:45 2014

@author: 
    Georg Ramer
    georg@ramer.at
    
    
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import ListedColormap, BoundaryNorm



 






class COSbase(object):
    def __init__(self, rawdata, rawtimes = None, nu_1=None, nu_2=None):
        """
        data is a np.ndarray with one spectrum per column. 
        rawtimes is a onedimensional array containing the perturbation
        """
        self.rawdata = np.array(rawdata)
        if rawtimes is None:
            rawtimes = np.arange(rawdata.shape[1])
        if nu_1 is None:
            nu_1 = np.arange(rawdata.shape[0])
        if nu_2 is None:
            nu_2 = nu_1
        self.nu_1 = nu_1
        self.nu_2 = nu_2
        self.rawtimes = np.array(rawtimes)
    

    def means(self):
        """
        calculates row wise means of the data.
        """
        return np.mean(self.rawdata, axis=1)
    

    def meancentered(self):
        """
            row wise mean center of the data.
        """
        return col_plus(self.rawdata, -self.means())
    
    def nu_to_index(nu, nu_1=True):
        if nu_1:
            return np.argmin(abs(self.nu_1 - nu))
        else:
            return np.argmin(abs(self.nu_2 - nu))
        
        
        
class COS(COSbase):
    """
        calculates 2D COS. Implementation following 
        Noda, I., & Ozaki, Y. (2005). Two-dimensional correlation spectroscopy:
        applications in vibrational and optical spectroscopy. John Wiley & Sons.
        
    """
    def __init__(self, *args,**kwargs):
        super(COS, self).__init__(*args,**kwargs)
    

    def sync(self):
        return 1/(len(self.rawtimes)-1) * np.dot(self.meancentered(), self.meancentered().T)
    

    def async(self):
        disturb = np.hstack((2*self.rawtimes[0] - self.rawtimes[1], self.rawtimes, 2*self.rawtimes[-1]-self.rawtimes[-2]))
        diffs = disturb[2:] - disturb[:-2]
        hino = np.zeros((self.rawtimes.shape[0],self.rawtimes.shape[0]))
        for row in range(hino.shape[0]):
            for col in range(hino.shape[1]):
                if row != col:
                    hino[row,col] = (disturb[col+2] - disturb[col])/\
                     (2*np.pi*(disturb[col+1] - disturb[row+1]))

        zwave = np.dot(hino, self.meancentered().T)
        return 1/(self.rawtimes[-1] - self.rawtimes[0])*np.dot(row_multiply(self.meancentered(),diffs),zwave).T 


class CODS(COSbase):
    """
        applies CDOS to data.
    """
    
    def __init__(self,*args,**kwargs):
        """
        rawtimes    have to be equally spaced in an ndarray
        """
        super(CODS,self).__init__(*args,**kwargs)
     

    def T_matrix(self):
        meancentered= self.meancentered()
        correlation = np.atleast_2d(np.sum( meancentered * meancentered, axis=1))
        return np.dot(correlation.T,correlation)
    

    def characteristic_times(self):
        kbar = self.characteristic_time_index()
        m = self.rawtimes.shape[0]
        tbar = (self.rawtimes[-1]-self.rawtimes[0]) * (kbar-1)/(m-1) + self.rawtimes[0]
        return tbar
        

    def characteristic_time_index(self):
        Abar = self.means()
        Awave = self.meancentered()
        m = self.rawtimes.shape[0]
        k = np.arange(m)
        kbar = np.sum(row_multiply(Awave, k),axis = 1)/(m*Abar)
        return kbar


    def async(self, epsilon = 10**-20):
        """
        return asynchronous CDOS spectrum of the data.
        epsilon is the cut off below which the average of a signal is considered zero
        """
        Abar = self.means()
        Awave = self.meancentered()
        zeroed_lines = (np.abs(Abar) <= np.amax(np.abs(Abar))*epsilon)
        Awave[zeroed_lines] = np.nan
        Abar[zeroed_lines] = np.nan
        
        m = self.rawdata.shape[1]
        ks = np.arange(1,m+1)
        summands = np.sum(row_multiply(Awave, ks), axis=1) / Abar
        summands = np.repeat(np.atleast_2d(summands).T,len(summands), axis=1)
        delta = (summands - summands.T) * self.T_matrix() / (m * (m-1))
        delta[np.isnan(delta)] = 0
        return delta
    
    def get_peak_vals(nu_1, nu_2, epsilon=10**-20, absolut=False):
        """
        gives itnerpretation according to Noda's rules. `absolut'
        """
        ind_1 = self.nu_to_index(nu_1)
        ind_2 = self.nu_to_index(nu_2, nu_1=False)
        sync_res = self.sync()[ind_1, ind_2]
        async_res = self.async()[ind_1, ind_2]
        return {"sync": sync_res, "async": async_res, "mean": np.means()[ind_1]}
        
        
        

    def sync(self, epsilon=10**-20):
        """
        return synchronous CDOS spectrum
        """
        sync_square = self.T_matrix()**2 - self.async(epsilon)**2
        sync_square[sync_square<0] = 0
        return np.sqrt(sync_square)




# helper functions that make life easier

def COS_plot(cos, means=None,ax=None, wn=None, lines = 5, title=None):
    """
    Creates nice (and correct) plot of 2D correlation calculations.
    
    `cos'    rectangular ndarray containing the correlation spectrum
    `means'    the average spectra to be displayed on the axis for orientation
    `ax'    the axis in which to plot the data (None new axis)
    `wn'    wavenumbers for the cos matrix, none for 
    """
    if wn is None:
        wn =np.arange(cos.shape[0])
    X,Y = np.meshgrid(wn,wn)
    
 
    if ax is None:
       ax = plt.figure().add_subplot(111)
    
    # flip COS so it fits  NODA Rules
    cos = cos.T
    X = X.T
    Y = Y.T
    #calculate values for layout
    absmax = np.nanmax(np.abs(cos))
    conts = np.linspace(-absmax,absmax,2*lines )
    bounds = (conts[0], conts[lines-1], conts[-1])
    cmap = ListedColormap(("lightgrey","white"))
    #plot grey negative areas
    ax.contourf(X,Y, cos , cmap =cmap,\
    norm = BoundaryNorm(bounds, cmap.N))
    #plot height lines
    ax.contour(X,Y,cos, conts, colors = "black",linestyles="solid")
    ax.set_aspect(1.)
    ax.set_xlim((wn[0],wn[-1]))
    ax.set_ylim((wn[0],wn[-1]))
    ax.set_xlabel(r"$\nu_1$")
    ax.set_ylabel(r"$\nu_2$")
    if means is not None:
        divider = make_axes_locatable(ax)
        specX = divider.append_axes("top", size=1.2, pad=0.1, sharex=ax)
        specY = divider.append_axes("right", size=1.2, pad=0.1, sharey=ax)
        if not title is None:
            specX.set_title( title)
        specX.plot(wn, means)
        specY.plot(means,wn)
   
    

def cdos(data,times=None, epsilon=None):
    """
    data is a np.ndarray with one spectrum per column. 
    times
    """
    COS = CODS(data, times)
    if epsilon is None:
        sync = COS.sync()
        async = COS.async()
    else:
        sync = COS.sync(epsilon)
        async = COS.async(epsilon)
    means = COS.means()
    return namedtuple("CDOS", ("sync", "async", "means"))(sync, async,means)
        
def cos(data,times=None):
    """
    data is a np.ndarray with one spectrum per column. 
    """
    cos2d = COS(data, times)
    sync = cos2d.sync()
    async = cos2d.async()
    means = cos2d.means()
    return namedtuple("COS", ("sync", "async", "means"))(sync, async,means)        

       

def gaussian (A,m, w, x):
    return np.atleast_2d(A * np.exp(-(x-m)**2 / w**2)).T


def testdata (lenX, lenT):
    k1 = 0.2
    k2 = 0.8
    c1 = lambda t:np.exp(-k1 * t)
    c2 = lambda t:-(np.exp(-k2*t)-np.exp(-k1*t))
    c3 = lambda t:1-(k2*np.exp(-k1*t) - k1 * np.exp(-k2*t))/(k2 - k1)
    x = np.linspace(1400,1000,lenX)
    t = np.linspace(0,10, lenT)
    res1 = np.hstack(map(lambda t: gaussian(c1(t), 1320, 10,x)+gaussian(c1(t), 1080, 5,x),t))
    res2 = np.hstack(map(lambda t: gaussian(c2(t), 1280, 20,x) +  gaussian(c2(t), 1120, 10,x),t))
    res3 = np.hstack(map(lambda t: gaussian(c3(t), 1240,15,x)+ gaussian(c3(t), 1160, 10,x),t))

    return res1 + res2 + res3, x

def col_plus(array, summand):
    return array + np.repeat(np.atleast_2d(summand).T, array.shape[1], axis=1)

def col_multiply(array, multiplicator):
    return array * np.repeat(np.atleast_2d(multiplicator).T, array.shape[1], axis=1)

def row_multiply(array, multiplicator):
    return array * np.repeat(np.atleast_2d(multiplicator), array.shape[0], axis=0)
    





