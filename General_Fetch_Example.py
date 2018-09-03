#!/usr/bin/env python 

from General_Fetch import Fetch
from obspy import UTCDateTime


def main():

	network='TA'
	station=None
	starttime = "2008-08-01"
	endtime = "2009-08-01"
	centercoords = [39, -111]
	minradius = 30
	maxradius = 120
	minmag = 6.0

	#Set up test instance. Ensure that the stations you request (can leave this blank to request all stations)
	#are within the boundary box given

	test = Fetch(network=network,station=station,starttime=UTCDateTime(starttime),endtime=UTCDateTime(endtime),\
		minlatitude=38,maxlatitude=40,minlongitude=-110,maxlongitude=-112)

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
	test.GetData(req_type='event',datadirpath='sheba_test')

	#Remove instrument response. Default is to displacement
	#test.CorrectResponse()

if __name__ == '__main__': 

	main()