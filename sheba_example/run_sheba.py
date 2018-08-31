#!/usr/bin/env python

#RMS 2018
#Simple example that shows how to build an input file for sheba in python, run it and obtain a result

import os 
import glob
import obspy as op

def main():

	#NOTE, need to set sac headers user1 and user3 for this to work.
	#These sac headers represent the start and end time of the user-picked 
	#window. We will want to come up with a good way of defining them automatically

	sacfiles = glob.glob('CAN.BH*')
	for fname in sacfiles:
		st = op.read(fname,format='SAC')
		tr = st[0]
		print(tr.stats.sac)

		#These SAC headers set the start and end points of the windows for 
		#cluster analysis. They will need to be determined and added to the headers each time
		tr.stats.sac.user0 = 605
		tr.stats.sac.user1 = 606
		tr.stats.sac.user2 = 622
		tr.stats.sac.user3 = 623
		tr.write(fname,format='SAC')

	#create sheba input file
	ofile = open('sheba.in','w')
	ofile.write('#sheba.in\nCAN\nBHE\nBHN\nBHZ\n1\n10 10\n4.0\n0\n0')
	ofile.close()

	#run sheba
	os.system('sheba_exec sheba.in')

	#the result is placed in CAN_sheba.result

if __name__ == '__main__':

	main()