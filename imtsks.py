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


# Parameters for making a nice plot
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=12)

# Making a global ray tracing model
from obspy.taup import TauPyModel
model = TauPyModel(model="iasp91")

# Client will be a global variable
client = Client("IRIS")


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
    goods = plt.savefig(eve.origins[0].time + '.jpg', format= 'JPEG', dpi=200)
    return goods

def stat():
    stat = 'ULN'
    return stat

def grabsta():
    # Function to get station information
    print("Getting Station Info")
    net = 'IU'
    #stat= stat()

    stime = UTCDateTime('2016-001T00:00:00.0')

    
    #inventory = client.get_stations(network=net, station=stat,
                                    #channel = 'BH*', level="response",
                                    #location="00",starttime=stime)
                                    
    #print("getting station coordinates")
    #station_coordinates = []
    #for network in inventory:
        #for station in network:
            #for channel in station:
                #station_coordinates.append((network.code, station.code, 
                                            #station.latitude, station.longitude, 
                                            #station.elevation,channel.azimuth))
    

    #print(station_coordinates)

    stalat= 47.8651
    stalon = 107.0532
    staelev = 1610.0
    return stime, stalat, stalon #station_coordinates

def definelist():
    # Makes a bunch of variables to add
    spliteve = []
    goodfast = []
    badfast = []
    gooddelay = []
    baddelay = []
    gc1 = []
    gc2 = [] 
    goodsplitf = []
    goodsplitd = []
    return spliteve, goodfast, badfast, gooddelay, baddelay, gc1, gc2, goodsplitf, goodsplitd

def grabSKSeves(stime, stalat, stalon, debug=True):
    # Grab all events for the station
    if debug:
        print("Getting Events")
    # Map of events
    cat = client.get_events(starttime=stime, minmagnitude=6., latitude=stalat, 
                            longitude=stalon, maxradius=130., minradius = 88.)
                            
    SKStimes = []
    Origintimes = []
    diss = []
    mags = []
    bazis = []
    if debug:
        print('Trying to find SKS events')
    # Now only grab good phases
    for eve in cat:
        try:
            (dis,azi, bazi) = gps2dist_azimuth(stalat, stalon, eve.origins[0].latitude,eve.origins[0].longitude)
            if debug:
                print('Here is the distance' + str(dis))
            arrivals = model.get_travel_times(source_depth_in_km=eve.origins[0].depth/1000., distance_in_degree=dis)
            for arrival in arrivals:
                if (arrival.name == 'SKS') or (arrival.name == 'SKKS'):
                    if debug:
                        print('Got an SKS event')
                    Origintimes.append(eve.origins[0].time)
                    SKStimes.append(arrival.time)
                    mags.append(eve.magnitudes[0].mag)
                    diss.append(dis)
                    if debug:
                        print(dis)
                    bazis.append(bazi)
        except:
            print('Got a problem')
    
    return diss, bazis, SKStimes, Origintimes, mags
                                           

    
def SKSpick(Origintime, arrivaltime):        
    skstime = (Origintime + arrivaltime)-20.
    etime = (Origintime + arrivaltime)+ 30.
    nstime = Origintime - 300.
    netime =  Origintime - 250.
    return skstime, etime, nstime, netime

def waveformdat(skstime, etime, nstime, netime, net, stat):            
    print("Getting Waveform data")
    try:
        st = client.get_waveforms(net, stat, '00', 'BH*', skstime, etime, attach_response = True)
        stN = client.get_waveforms(net, stat, '00', 'BH*', nstime, netime, attach_response = True)
    except:
        pass
    return st, stN


def process(st, stN):
    st.remove_sensitivity()
    stN.remove_sensitivity()
    st.filter('bandpass', freqmin= 0.06, freqmax = 0.09)
    stN.filter('bandpass', freqmin= 0.06, freqmax = 0.09)
    st.taper(0.05)
    stN.taper(0.05)
    return st, stN

def rotate(st, stN, bazi):
    #This is a hack for rotation
    #print("Rotating Data")
    st, stN = process(st, stN)
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
    return st, stN
    
def snr(st, stN):    
    stS = st.copy()
    stN.std()
    S = stS[0].std()
    N = stN[0].std()
    snr = S/N
    return snr


print("Defining our Functions")
def getsksstuff(st, debug = False):
    # debug is a flag to help with debugging.  So if we switch it to false
    # then all the print statements turn off
    if debug:
        print(st)
    phis = np.deg2rad(np.arange(-90.,90.,1.))
    if debug:
        print(phis)
    # Start with a huge lambda and then try to update it
    minlam = 50000.
    # Make a copy of st and window it at 20 s from the start
    stF = st.copy()
    stF = stF.trim(endtime = stF[0].stats.starttime + 20.)
    dt = 0.
    # Slide through the slow axis with steps of .1
    #print("Starting Calculation")
    for stW in st.slide(window_length=20., step=.1):
        dt += .1
        for phi in phis:
            # Rotate the data by the given phi
            comp1 = np.cos(phi)*stF[0].data - np.sin(phi)*stF[1].data
            comp2 = np.sin(phi)*stW[0].data + np.cos(phi)*stW[1].data
            # comp1 is our "fast" component and comp2 is slow
            comp1 += -np.mean(comp1)
            comp2 += -np.mean(comp2)
            # Make our coveraiance matrix
            c11 = sum(comp1*comp1)
            c22 = sum(comp2*comp2)
            c12 = sum(comp1*comp2)
            cors = np.matrix([[c11,c12],[c12,c22]])
            # Compute some eigen values
            lamb, v = np.linalg.eig(cors)
            lamb = lamb.real.astype('float')
            if lamb[0]*lamb[1] < minlam:
                dtgood = dt
                phigood = phi
                #print("Found Good Fit")
                if debug:
                    print('New minimum: ' + str(lamb[0]*lamb[1]))
                    print('New phi: ' + str(phigood) + ' New dt: ' + str(dtgood))
                minlam = lamb[0]*lamb[1]
    if debug:
        print('Returning: ' + str(phigood) + ' ' + str(dtgood))
    phigood = np.rad2deg(phigood)
    print('Here is our direction: ' + str(phigood) + ' here is our delay: ' + str(dtgood))
    stW = st.copy()
    stW.trim(stW[0].stats.starttime + dtgood, stW[0].stats.starttime + 20. + dtgood)
    comp1 = np.cos(phigood)*stF[0].data - np.sin(phigood)*stF[1].data
    comp1 *= 1./np.max(np.abs(comp1))
    comp2 = np.sin(phigood)*stW[0].data + np.cos(phigood)*stW[1].data
    comp2 *= 1./np.max(np.abs(comp2))
    #print("Returning Splitting Variables")
    return phigood, dtgood, stF, stW, comp1, comp2
    
    
def goodsnr(snr):                
    snr = snr(st, stN)
    spliteve, goodfast, badfast, gooddelay, baddelay, gc1, gc2, goodsplitf, goodsplitd = definelist()
    try:
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
    except:
        pass
    return spliteve, goodfast, gooddelay, gc1, gc2, badfast, baddelay
#plt.scatter(gooddelay, goodfast, color = 'c')
#plt.scatter(baddelay, badfast, color = 'r')        
#plt.savefig('SNR_more_than_5.jpg', format = 'JPEG')


def constrain(spliteve, goodfast, gooddelay, gc1, gc2, badfast, baddelay):
    spliteve, goodfast, gooddelay, gc1, gc2, badfast, baddelay = goodsnr(snr)
    try:
        for eveidx, eve in enumerate(spliteve):
            if eveidx == 0:
                goodcompf = gc1[eveidx]
                goodcompw = gc2[eveidx]
                goodsplitf.append(goodfast[eveidx])
                goodsplitd.append(gooddelay[eveidx])
                print("Reference Added")
            else:
                compf = gc1[eveidx]
                compw = gc2[eveidx]
                fcompr = np.abs(compf-goodcompf/goodcompf)
                wcompr = np.abs(compw-goodcompw/goodcompw)
                if fcompr[eveidx] < 0.1: 
                    if wcompr[eveidx]<0.1:
                        goodsplitf.append(goodfast[eveidx])
                        goodsplitd.append(gooddelay[eveidx])
                    print("Good Comparison, another event added")
                else:
                    print("Bad Comparison, data thrown out")
    except:
        pass
    return goodsplitf, goodsplitd

def runprogram(stat):
    net ='IU'
    # Make a bunch of empty variables
    spliteve, goodfast, badfast, gooddelay, baddelay, gc1, gc2, goodsplitf, goodsplitd = definelist()
    
    # Grab station coordinates
    stime, stalat, stalon = grabsta()
    
    # Get parameters for SKS stuff for the station
    diss, bazis, SKStimes, Origintimes, mags = grabSKSeves(stime, stalat, stalon)

    # For each event we now have the necessary info
    for fiple in zip(diss, bazis, SKStimes, Origintimes, mags):
        skstime, etime, nstime, netime = SKSpick(fiple[3],fiple[2])
        
        st, stN = waveformdat(skstime, etime, nstime, netime, net, stat)
        st, stN = process(st, stN)
        st, stN = rotate(st, stN, fiple[1])
        snrVAL = snr(st, stN)
        spliteve, goodfast, gooddelay, gc1, gc2, badfast, baddelay = goodsnr(snrVAL)
        goodsplitf, goodsplitd = constrain(spliteve, goodfast, gooddelay, gc1, gc2, badfast, baddelay)
    return goodsplitf, goodsplitd



# Here is the start of your program

if __name__ == "__main__":

# This is a big of an odd way to do this (you are mixing functional programming with procedural)
    stat = stat()

    runprogram(stat)



#pool = Pool(10)
#pool.map(runprogram, stat)
#pool.close()
                
#print(goodsplitf, goodsplitd)                
        
        ##Add the results to a file with magnitude, time of event, distance, backazimuth
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
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
