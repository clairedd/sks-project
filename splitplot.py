#!/usr/bin/env python

from obspy.core import read, UTCDateTime
import matplotlib.pyplot as plt
import numpy as np
import glob
import sys


fig = plt.figure(1,figsize=(8,8))

f= open('ULNsplits','r')

theta =[]
r =[]
for line in f:
    line = ' '.join(line.split())
    line = line.split(' ')
    print(line)
    theta.append(np.deg2rad(float(line[5])))
    theta.append(np.deg2rad(float(line[5])+180.))
    r.append(float(line[8]))
    r.append(float(line[8]))
print(theta)
theta= np.asarray(theta)
r = np.asarray(r)
his, thetas, rs = np.histogram2d(theta, r, bins=(36., 10.), range=[[0.,2*np.pi],[0.,1.]],normed=True)
his = his.T
THE, R = np.meshgrid(thetas, rs)
ax = plt.subplot(1,1,1,projection='polar')
plt.pcolormesh(THE,R,his)
plt.yticks([])
plt.ylim((0.,1.))
plt.tight_layout()    
plt.show()

#print(r)
#print(theta)
#fig = plt.figure(1)
#ax = plt.subplot(111,projection='polar')
#plt.plot(theta,r)
#plt.show()
