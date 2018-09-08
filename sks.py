from obspy.core import read, UTCDateTime
from obspy.clients.fdsn import Client
import sys
import matplotlib.pyplot as plt
import numpy as np
from obspy.signal.cross_correlation import xcorr
import matplotlib as mpl
from obspy.geodetics.base import gps2dist_azimuth
from obspy.core.event import Amplitude
from multiprocessing import Pool

#mpl.rc('font',family='serif')
#mpl.rc('font',serif='Times')
#mpl.rc('text', usetex=True)
#mpl.rc('font',size=12)

from obspy.taup import TauPyModel
model = TauPyModel(model="iasp91")

print("Defining our Functions")
def getsksstuff(st, debug = True):
    # debug is a flag to help with debugging.  So if we switch it to false
    # then all the print statements turn off
    if debug:
        print(st)
    phis = np.deg2rad(np.arange(-90.,90.,1.))
    if debug:
        print(phis)
    # Start with a huge lambda and then try to update it
    minlam = 50000.
    if debug:
        print('Lambda Assigned')
    # Make a copy of st and window it at 6 s from the start
    stF = st.copy()
    if debug:
        print('stF has been initialized')
    stF = stF.trim(endtime = stF[0].stats.starttime + 5.)
    dt = 0.
    #phigood = 0.
    #dtgood = 0.
    # Slide through the slow axis with steps of .1
    if debug:
        print("Starting Calculation")
    for stW in st.slide(window_length=6., step=.1):
        dt += .1
        if debug:
            print('Rotating around all Phis')
        for phi in phis:
            # Rotate the data by the given phi
            comp1 = np.cos(phi)*stF[0].data - np.sin(phi)*stF[1].data
            comp2 = np.sin(phi)*stW[0].data + np.cos(phi)*stW[1].data
            # comp1 is our "fast" component and comp2 is slow
            if debug:
                print('Fast vector is done')
                print('Slow vector is done')
            comp1 += -np.mean(comp1)
            comp2 += -np.mean(comp2)
            if debug:
                print('Data is now demeaned')
            # Make our coveraiance matrix
            c11 = sum(comp1*comp1)
            c22 = sum(comp2*comp2)
            c12 = sum(comp1*comp2)
            cors = np.matrix([[c11,c12],[c12,c22]])
            if debug:
                print('Covariance Matrix computed. Now getting eigenvalues')
            # Compute some eigen values
            lamb, v = np.linalg.eig(cors)
            lamb = lamb.real.astype('float')
            if debug:
                print('Eigenvalue found. Determining splitting parameters')
            if lamb[0]*lamb[1] < minlam:
                dtgood = dt
                phigood = phi
                #print("Found Good Fit")
                if debug:
                    print('New minimum: ' + str(lamb[0]*lamb[1]))
                    #print('New phi: ' + str(phigood) + ' New dt: ' + str(dtgood))
                minlam = lamb[0]*lamb[1]
                    #if debug:
        #print('Returning: ' + str(phigood) + ' ' + str(dtgood))
    phigood = np.rad2deg(phigood)
    print('Here is our direction: ' + str(phigood) + ' here is our delay: ' + str(dtgood))
    if debug:
        print('Now saving our data')
    stW = st.copy()
    stW.trim(stW[0].stats.starttime + dtgood, stW[0].stats.starttime + 1. + dtgood)
    comp1 = np.cos(phigood)*stF[0].data - np.sin(phigood)*stF[1].data
    comp1 *= 1./np.max(np.abs(comp1))
    comp2 = np.sin(phigood)*stW[0].data + np.cos(phigood)*stW[1].data
    comp2 *= 1./np.max(np.abs(comp2))
    #print("Returning Splitting Variables")
    return phigood, dtgood, stF, stW, comp1, comp2

def plots(stW, stF):
    fig = plt.figure()
    t = np.arange(0., 40.1,0.1)
    plt.subplot(311)
    plt.plot(t,stF[0].data, color = 'c')
    plt.plot(t, stW[0].data, color = 'burlywood')
    plt.subplot(312)
    plt.plot(t, comp1, color = 'c')
    plt.plot(t, comp2, color = 'burlywood')
    plt.subplot(313)
    plt.plot(t, stW[0].data, color = 'c')
    plt.plot(t, phigood*stW[1].data, color = 'burlywood')
    #goods = plt.savefig(eve.origins[0].time + '.jpg', format= 'JPEG', dpi=200)
    return goods


print("Getting Station Info")
net = 'IU'
stat= 'ULN'

stime = UTCDateTime('2016-001T00:00:00.0')

#stalat = 40.128
#stalon= -107.51

client = Client("IRIS")
inventory = client.get_stations(network=net, station=stat,
                                    channel = 'BH*', level="response",
                                    location="00",starttime=stime)
                                    
print("getting station coordinates")
station_coordinates = []
for network in inventory:
    for station in network:
        for channel in station:
            station_coordinates.append(network.code)
            station_coordinates.append(station.code)
            station_coordinates.append(station.latitude)
            station_coordinates.append(station.longitude)
            station_coordinates.append(station.elevation)
            station_coordinates.append(channel.azimuth)


#print(station_coordinates)

for ele in enumerate(station_coordinates):
    stalat= station_coordinates[2]
    stalon = station_coordinates[3]
    staelev = station_coordinates[4]

print(stalat)
print(stalon)

spliteve = []
goodfast = []
badfast = []
gooddelay = []
baddelay = []
gc1 = []
gc2 = [] 
goodsplitf = []
goodsplitd = []


print("Getting Events")
# Map of events
cat = client.get_events(starttime=stime, minmagnitude=6., latitude=stalat,
                        longitude=stalon, maxradius=130., minradius = 88.)
                                           
#def splittingparam(cat):
for eve in cat:
    #print(eve)
    (dis,azi, bazi) = gps2dist_azimuth(stalat, stalon, eve.origins[0].latitude,eve.origins[0].longitude)
    arrivals = model.get_travel_times(source_depth_in_km=eve.origins[0].depth/1000., distance_in_degree=dis)
    #st = Stream()
    for arrival in arrivals:
        #print(arrival)
        if arrival.name == 'SKS':
            #print(eve)
            #print('Windowing SKS event time')
            skstime = (eve.origins[0].time + arrival.time)- 1.
            etime = (eve.origins[0].time + arrival.time)+ 5.
            nstime = eve.origins[0].time - 300.
            netime = eve.origins[0].time - 294.
            print("Getting Waveform data")
            st = client.get_waveforms(net, stat, '00', 'BH*', skstime, etime, attach_response = True)
            stN = client.get_waveforms(net, stat, '00', 'BH*', nstime, netime, attach_response = True)
            #print("Processing our data")
            st.remove_sensitivity()
            stN.remove_sensitivity()
            st.filter('bandpass', freqmin= 0.06, freqmax = 0.09)
            stN.filter('bandpass', freqmin= 0.06, freqmax = 0.09)
            st.taper(0.05)
            stN.taper(0.05)
                #st.plot()
                #stN.plot()
            #This is a hack for rotation
            #print("Rotating Data")
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
            stN.std()
            S = stS[0].std()
            N = stN[0].std()
            #print("Calculating Signal-to-Noise Ratio")
            snr = S/N
            print(snr)
                
            if snr>5:
                fast, delay, stF, stW, comp1, comp2 = getsksstuff(st)
                spliteve.append(eve)
                goodfast.append(fast)
                gooddelay.append(delay)
                gc1.append(comp1)
                gc2.append(comp2)
                print("Good Split")
            else:
                fast, delay, stF, stW, comp1, comp2 = getsksstuff(st)
                badfast.append(fast)
                baddelay.append(delay)
                print("Bad Split")
                

        
plt.scatter(gooddelay, goodfast, color = 'c')
plt.scatter(baddelay, badfast, color = 'r')
plt.xlabel("Delay Time (Seconds)")
plt.ylabel("Fast Direction (Degrees)")
plt.xticks()
plt.yticks()
plt.show()
#plt.savefig('SNR_more_than_5.jpg', format = 'JPEG')

#for eveidx, eve in enumerate(spliteve):
#    if eveidx == 0:
#        goodcompf = gc1[eveidx]
#        goodcompw = gc2[eveidx]
#        goodsplitf.append(goodfast[eveidx])
#        goodsplitd.append(gooddelay[eveidx])
#        print("Reference Added")
#    else:
#        compf = gc1[eveidx]
#        compw = gc2[eveidx]
#        fcompr = np.abs(compf-goodcompf/goodcompf)
#        wcompr = np.abs(compw-goodcompw/goodcompw)
#        if fcompr[eveidx] < 0.1:
#            if wcompr[eveidx]<0.1:
#                goodsplitf.append(goodfast[eveidx])
#                goodsplitd.append(gooddelay[eveidx])
#            print("Good Comparison, another event added")
#        else:
#            print("Bad Comparison, data thrown out")
#
#plt.savefig('Data_Within_10p.jpg', format = 'JPEG')
#print("Getting Events")
# Map of events
#cat = []

#pool = Pool(10)
#pool.map(splittingparam,cat)
#pool.close()
                
        ##trim data for S-N-R calculation
        
        ##Add the results to a file with magnitude, time of event, distance, backazimuth
        
        ##For each one, plot of the SKS splitting after it's been rotated. 
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        #arrival = arrivals[0]
        #print(arrival)
        #estime = eve.origins[0].time
        #sttemp = st.copy()
        #sttemp += stG.copy()
        #sttemp.trim(estime-20.,estime + 4020.)
        #sttemp.detrend('linear')
        #for tr in sttemp.select(channel='LGZ'):
            #tr.data *= 10.**-18
        #for tr in sttemp.select(channel='LHZ'):
            #tr.simulate(paz_remove=paz)
            ##tr.data *= 1./(.98*20.*1500.*(2**24)/40.)
            ##print(tr.data)
        #sttemp.filter('bandpass',freqmin=0.02, freqmax = 0.15)
        #print(sttemp)
        #lag, corr = xcorr(sttemp[0],sttemp[1], 500)
        #print('Here is the lag: '  + str(lag))
        #print('Here is the corr: ' + str(corr))
        #sttemp[1].stats.starttime += lag
        #sttemp.trim(estime + arrival.time - 60.,estime + arrival.time + 120.)
        #rel = sttemp[0].std()/sttemp[1].std()
        #fig = plt.figure(1)
        #t = np.arange(0, sttemp[0].stats.npts)
        #plt.plot(t,sttemp[0].data*10**9, label=sttemp[0].id)
        #plt.plot(t,sttemp[1].data*10**9, label=sttemp[1].id)
        #mag = eve.magnitudes[0].mag
        #magstr = eve.magnitudes[0].magnitude_type
        
        #plt.title(str(estime.year) + ' ' + str(estime.julday).zfill(3) + ' ' + str(estime.hour).zfill(2) + ':' + str(estime.minute).zfill(2) 
                    #+ ' ' + magstr + '=' + str(mag) + ' lag=' + str(lag) + ' s Rel. Gain: '+ str(round(rel,2)))
        #plt.xlim((min(t),max(t)))
        #plt.xlabel('Time (s)')
        #plt.ylabel('Acceleration $(nm/s^2)$')
        #plt.legend()
        ##plt.show()
        #plt.savefig(str(estime.year) + '_' + str(estime.julday).zfill(3) + '.jpg', format='JPEG')
        #plt.clf()
        #plt.close()
    #except:
        #print('Problem')
