# -*- coding: utf-8 -*-
"""
Created on Mon May 26 15:51:17 2014

@author: 
    Georg Ramer
    georg@ramer.at
    
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from COS import cos, testdata, COS_plot


# create testdata
data, wn =  testdata(200,400)

plt.figure("spectra")
plt.plot(wn, data)
plt.xlim((wn[0],wn[-1]))
plt.figure()
plt.contourf( range(400),wn,data)

#calculate 2D-CoS
result = cos(data =data)

sync = result.sync
async = result.async  



#display the 2D-CoS result
axsync = plt.figure("sync").add_subplot(111)
COS_plot(sync,result.means, ax=axsync, wn = wn, lines = 10)

axasync = plt.figure("async").add_subplot(111)
COS_plot(async,result.means, ax=axasync, wn = wn, lines=10)


plt.show()
