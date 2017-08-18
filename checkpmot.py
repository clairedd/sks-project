#!/usr/bin/env python

from obspy.core import read, UTCDateTime
import matplotlib.pyplot as plt
import numpy as np
import glob
import sys

dirs = glob.glob('/msd/*/2016/')
stas =[]
for sta in dirs:
    stas.append(sta.split('/')[2])


for idx, sta in enumerate(stas):
    print('On sta : ' + str(idx) + ' of ' + str(len(stas)))
    try:
        st = read('/msd/' + sta + '/2017/10*/00_LH*')

        for tr in st.select(channel='LH0'):
            st.remove(tr)
        for tr in st:
            tr.data = tr.data.astype(np.float64)    
        st.merge(fill_value=0.)
        st.detrend('constant')
        st.sort()


        fpairs = [[1./8.,1./4.],[1./22.,1./18.],[1./110., 1./90.],[0.,0.]]
        fig = plt.figure(1,figsize=(8,8))
        for idx, fpair in enumerate(fpairs):

            stTemp = st.copy()
            if fpair[0] != 0.:
                stTemp.filter('bandpass',freqmin=fpair[0],freqmax=fpair[1])
            stTemp.taper(0.05)
            stTemp.normalize(global_max=True)
            theta= np.arctan2(stTemp[0].data,stTemp[1].data)
            r = np.sqrt((stTemp[0].data)**2 + (stTemp[1].data)**2)
            r *= 1./np.max(r)

            his, thetas, rs = np.histogram2d(theta, r, bins=(360., 100.),range=[[-np.pi, np.pi],[0., 2*np.mean(r)]], normed=True)
            his = his.T
            THE, R = np.meshgrid(thetas, rs)
            ax = plt.subplot(2,2,(idx+1),projection='polar')
            if fpair[0] != 0.:
                plt.title('Band-Pass:' + str(int(1./fpair[1])) + ' to ' + str(int(1./fpair[0])))
            else:
                plt.title('Raw')
            plt.pcolormesh(THE,R,his)
            plt.ylim((0,2*np.mean(r)))
            plt.yticks([])
            
        plt.tight_layout()    

        plt.savefig(sta + '_2017_10.jpg',format='jpeg')
        plt.clf()
        plt.close()
    except:
        print('problem with ' + sta)



#print(r)
#print(theta)
#fig = plt.figure(1)
#ax = plt.subplot(111,projection='polar')
#plt.plot(theta,r)
#plt.show()
