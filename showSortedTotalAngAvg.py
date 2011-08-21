#!/usr/bin/env python

# Usage:
# In this directory, type:
#    ./viewRun.py -rxxxx
# For details, type 
#	 ./viewRun.py --help
# where rxxxx is the run number of hits and nonhits found using the hitfinder executable. 
# By default, this script looks into the h5 files that are in the appropriate rxxxx directory
#

import numpy as N
import h5py as H
import matplotlib
import matplotlib.pyplot as P
import sys
import os
import re
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-r", "--run", action="store", type="string", dest="runNumber", 
					help="run number you wish to view", metavar="xxxx", default="")
parser.add_option("-f", "--inputFile", action="store", type="string", dest="fileList", 
				  help="text file containing h5 ang_avg files", metavar="file.txt", default="")
parser.add_option("-o", "--outputdir", action="store", type="string", dest="outputDir", 
				  help="store output in this directory (default: output_rxxxx)", metavar="myDir", default="output")		
(options, args) = parser.parse_args()

########################################################
# Edit this variable accordingly
# Files are read for source_dir/runtag and
# written to write_dir/runtag.
# Be careful of the trailing "/"; 
# ensure you have the necessary read/write permissions.
########################################################
source_dir = "/reg/d/psdm/cxi/cxi25410/scratch/cleaned_hdf5/"
ang_avg_dir = "/reg/d/psdm/cxi/cxi25410/scratch/cleaned_hdf5/ice_runs/"
runtag = "r%s"%(options.runNumber)
write_dir = options.outputDir + '_' + runtag + '/'
h5files= []

if(options.fileList is ""):
	searchDir = ang_avg_dir + runtag
	print "Now examining H5 files in %s/ ..."%(searchDir)	
	searchstring="[a-zA-Z0-9\_]+"+runtag+"[a-z0-9\_]+-angavg.h5"
	h5pattern = re.compile(searchstring)
	h5files = [x for x in os.listdir(searchDir) if h5pattern.findall(x)]
else:
	f = open(options.fileList, 'r')
	h5files = f.read().split()
	f.close()

numFiles = len(h5files)
print "Found %d angavg files"%(numFiles)
integratedIntens = N.zeros(numFiles)
for i in range(numFiles):
	fullFilePath = ang_avg_dir+runtag+"/"+h5files[i]
	f = H.File(fullFilePath, 'r')
	integratedIntens[i] = N.array(f['/data/data']).sum()
	f.close()

ordering = integratedIntens.argsort()
P.plot(integratedIntens[ordering])
P.xlabel("sorted frame number")
P.ylabel("integrated radial intensities of each run")
P.show()
orderedFiles = N.array(h5files)[ordering]

if not os.path.exists(write_dir):
	os.mkdir(write_dir)
recordTag = write_dir + runtag + "_orderedFiles.txt"
orderedFiles.tofile(recordTag, sep="\n") 
print "wrote to file "+ recordTag
