#!/usr/bin/env python 

from General_Fetch import Fetch
from obspy import UTCDateTime


def main():

	network='TA'
	station='E25K' #None
	starttime = "2016-08-23"
	endtime = "2016-09-28"
	centercoords = [58, -145]
	minradius = 30
	maxradius = 120
	minmag = 6.0

	#Set up test instance. Ensure that the stations you request (can leave this blank to request all stations)
	#are within the boundary box given

	test = Fetch(network=network,station=station,starttime=UTCDateTime(starttime),endtime=UTCDateTime(endtime),\
		minlatitude=55,maxlatitude=70,minlongitude=-160,maxlongitude=-140)

	#Fetch the details of the events that are within the request region 
	test.fetchEvents(centercoords=centercoords,minradius=minradius,maxradius=maxradius,minmag=minmag)

	#Write the event details to a file
	test.writeEvents()

	#Get all the station information
	test.fetchInventory()

	#Write the stations to a file
	test.writeStations()

	#Get the data and store in a file called waveforms_dir. You can give it any path
	#Note that the data are downloaded as mseed files, one per component
	print("Getting data")
	test.GetData(req_type='event',datadirpath='waveforms_dir')

	#Remove instrument response. Default is to displacement
	test.CorrectResponse()

if __name__ == '__main__': 

	main()