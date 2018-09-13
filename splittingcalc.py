#!/usr/bin/env python

from obspy.core import read, UTCDateTime
from obspy.clients.fdsn import Client
import sys
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as lin
from obspy.signal.cross_correlation import xcorr
import matplotlib as mpl
from obspy.geodetics.base import gps2dist_azimuth
from obspy.core.event import Amplitude
from multiprocessing import Pool
from obspy.taup import TauPyModel

#Call model which we will use to compute delay times
model = TauPyModel(model="iasp91")

#Data parameters: Where to get the data and for what time period
client = Client("IRIS")
net = 'IU'
stat= 'ULN'
stime = UTCDateTime('2016-001T00:00:00.0')

#Event Parameters: Mw and distance of events you want to pull from
#Minimum Event Magnitude
minmag = 6.0
#Maximum Event Magnitude
maxmag = 8.0
#Minimum distance (in degrees) from center to pull events from
mindeg = 88.0
#Maximum distance (in degrees) from center to pull events from
maxdeg = 130.0
#What phase you are trying to use to determine your splits.
phasetype = 'SKS'

#Data Processing Parameters: Frequencies over which you want to filter your data.
freqmin = 0.06
freqmax = 0.09

#Window length. This is set for now, but ideally will implement the clustering method of SHEBA.
pickwin = int(6)

"""
From this point on, the code should not need to be manipulated by the user. All calculations
set forth should run based off of the parameters given above.
"""

#Initialize the vectors of phi and dt values that we will use to computed
#our energy matrix.
phis = np.deg2rad(np.arange(-90.,90.,1.))
dts = np.arange(0.0, 5.0, 0.1)
#Initialize vectors that will be used to store data later on.
stF = []
stW = []
spliteve = []
goodfast = []
badfast = []
gooddelay = []
baddelay = []
gc1 = []
gc2 = []
goodsplitf = []
goodsplitd = []


#Defining the covariance values that we will use later
def covmatrix(phis):
    M = np.zeros((2,2,len(phis)))
    M[0,0,:] = np.cos(phis)
    M[0,1,:] = -np.sin(phis)
    M[1,0,:] = np.sin(phis)
    M[1,1,:] = np.cos(phis)
    return M

#Function that will initialize all arrays we will use in the calculation.
def initial(phis,dts):
    #Energy matrix
    Emat = np.zeros((len(phis),len(dts)))
    #Fast lambda values
    lam1 = np.zeros((len(phis),len(dts)))
    #Slow lambda values
    lam2 = np.zeros((len(phis),len(dts)))
    M = covmatrix(phis)
    return Emat, lam1, lam2, M

#Does the calculations based off the minimum energy method of Silver and Chan, 1991.
def Ematrix(phis, dts):
    Emat, lam1, lam2, M = initial(phis,dts)
    #Loop through all possible values of phi
    for p in np.arange(1,len(phis)):
        Ms = M[:,:,p]
        dataset = [stF,stW]
        dataset = np.array([dataset])
        testfast = Ms * dataset.T
        #Compute the temporary value of the fast direction for that cell.
        #Need to run through splitlab code to figure out what pickwin is
        tmpF = testfast[0, pickwin]
        #Inverse of the covariance matrix
        Minv = lin.inv(Ms)
        #Now loop through all possible values of dt.
        for t in np.arange(1,len(dts), -0.1):
            shift = dts[t]
            tmpS = testfast[1, pickwin+shift]
            #Splitlab writes Minv as MM' in their code. Unsure if M should be multiplied
            #with Minv here or not.
            corr_FS = Minv * [tmpF, tmpS]
            Emat[p,t] = corr_FS[0,:] * lin.inv(corr_FS[1,:])
            comp1 = corr_FS[0,:]
            comp2 = corr_FS[1,:]
            comp1 += -np.mean(comp1)
            comp2 += -np.mean(comp2)
            # Make our coveraiance matrix
            c11 = sum(comp1*comp1)
            c22 = sum(comp2*comp2)
            c12 = sum(comp1*comp2)
            covar = np.matrix([[c11,c12],[c12,c22]])
            #Compute the values of lambda
            lamb = lin.eig(covar)
            lamb = lamb.real.astype('float')
            lam1[p,t] = lamb[1,1]
            lam2[p,t] = lamb[0,0]
    return Emat, tmpF, tmpS, lam1, lam2

#Determines the correct phi and dt given the data.
def splitparams(phis, dts):
    Emat, tmpF, tmpS, lam1, lam2 = Ematrix(phis, dts)
    #Find the minimum value of the energy matrix.
    [realphi, realdt] = np.amin(Emat)
    #Find the minimum value of the slow and fast lambdas.
    ind = np.amin(l1*l2)
    #Give the coordinates of the minimum lambda values.
    [indpEV, indtEV]= np.unravel_index(size(l2), ind)
    return realphi, realdt, indpEV, indtEV

#Function for converting the index of the phi and dt values into degrees and seconds.
def convert(phis, dts, stF, stW):
    realphi, realdt, indpEV, indtEV = splitparams(phis, dts)
    delay = dts(realdt)
    backazi = np.rad2deg(phis(realphi))
    if backazi>90:
        backazi = backazi - 180.0
    return delay, backazi

#Pull station information from the selected client and define the GPS coordinates and
#elevation of that particular station.
def stainfo(client,net, stat, stime):
    #Get data from the IRIS DMC
    inventory = client.get_stations(network=net, station=stat,
                                        channel = 'BH*', level="response",
                                        location="00",starttime=stime)
    for network in inventory:
        for station in network:
            for channel in station:
                stalat= station.latitude
                stalon = station.longitude
                staelev = station.elevation
    return stalat, stalon, staelev

#Get events from the DMC and determine which ones have SKS splits.
def getphaseevents(client, net, sta, stime):
    stalat, stalon, staelev = stainfo(client, net, stat, stime)
    #Get a catalogue of events for the station(s) of choice and in a given area.
    cat = client.get_events(starttime=stime, minmagnitude=minmag, latitude=stalat,
                            longitude=stalon, maxradius=maxdeg, minradius = mindeg)
    for eve in cat:
        (dis,azi, bazi) = gps2dist_azimuth(stalat, stalon, eve.origins[0].latitude,eve.origins[0].longitude)
        arrivals = model.get_travel_times(source_depth_in_km=eve.origins[0].depth/1000., distance_in_degree=dis)
        #Loop over all arrivals that ObsPy detects for a given event, and pick out the coordinates
        #that have SKS arrivals
        for arrival in arrivals:
            #print(arrival)
            if arrival.name == phasetype:
                #Create a time window in the data around the phase arrival
                phasetime = (eve.origins[0].time + arrival.time)- 1.
                etime = (eve.origins[0].time + arrival.time)+ 5.
                nstime = eve.origins[0].time - 300.
                netime = eve.origins[0].time - 294.
                print("Getting Waveform data")
                st = client.get_waveforms(net, stat, '00', 'BH*', phasetime, etime, attach_response = True)
                stN = client.get_waveforms(net, stat, '00', 'BH*', nstime, netime, attach_response = True)
                #Process the data
                st.remove_sensitivity()
                stN.remove_sensitivity()
                st.filter('bandpass', freqmin= freqmin, freqmax = freqmax)
                stN.filter('bandpass', freqmin= freqmin, freqmax = freqmax)
                st.taper(0.05)
                stN.taper(0.05)
                    #st.plot()
                    #stN.plot()
                #This is a hack for rotating the data into the radial and transverse directions.
                for tr in st.select(channel = 'BH1'):
                    tr.stats.channel = 'BHN'
                for tr in st.select(channel = 'BH2'):
                    tr.stats.channel = 'BHE'
                st.rotate('NE->RT', bazi)
                for tr in stN.select(channel = 'BH1'):
                    tr.stats.channel = 'BHN'
                for tr in stN.select(channel = 'BH2'):
                    tr.stats.channel = 'BHE'
                stN.rotate('NE->RT', bazi)
                stS = st.copy()
                #Calculate the signal-to-noise ratio based off of the standard deviations of
                #the SKS phase window and that of a random selection of "noise" for that given
                #station.
                S = stS[0].std()
                N = stN[0].std()
                #print("Calculating Signal-to-Noise Ratio")
                snr = S/N
                if snr >5:
                    stF = st.copy()
                    stW = st.copy()
                    delay, backazi = convert(phis, dts, stF, stW)
                    spliteve.append(eve)
                    goodfast.append(backazi)
                    gooddelay.append(delay)
                    print("Good Split")
                    print('Delay Time is:  '+ delay + '   '+'Phi is:  '+ backazi)
                else:
                    stF = st.copy()
                    stW = st.copy()
                    delay, backazi = convert(phis, dts, stF, stW)
                    spliteve.append(eve)
                    goodfast.append(backazi)
                    gooddelay.append(delay)
                    print("Bad Split")
    return eve, st, snr, spliteve, goodfast, gooddelay


#Using the functions defined above, we now compute "good" and "bad" measurements on the basis of
#signal-to-noise ratios. At some point, I would like particle motion to also be taken into account
#but that is a future project.
eve, st, snr, spliteve, goodfast, gooddelay = getphaseevents(client, net, stat, stime)
