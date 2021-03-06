#!/usr/bin/env python

import numpy as N
import h5py as H
import glob as G
import matplotlib
import matplotlib.pyplot as P
import sys, os, re, shutil, subprocess
from myModules import extractDetectorDist as eDD
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-r", "--run", action="store", type="string", dest="runNumber", help="run number you wish to view", metavar="xxxx", default="")
parser.add_option("-W", "--weakHitsTreatment", action="store", type="int", dest="weakHitsTreatment", help="(default)0, stores ang_avg.h5 file names as rxxxx_weakAvgFiles.txt;\n1, also shows averages both ang_avg and 2D patterns (slow);", metavar="0 or 1", default=0)
parser.add_option("-S", "--strongHitsTreatment", action="store", type="int", dest="strongHitsTreatment", help="(default)0, stores ang_avg.h5 filenames as rxxxx_strongAvgFiles.txt in output dir;\n1, also shows averages both ang_avg and 2D patterns (slow);", metavar="0 or 1", default=0)
parser.add_option("-c", "--copyFiles", action="store_true", dest="store_files", help="copy *.angavg.h5 files into output directory",default=False)
parser.add_option("-o", "--outputdir", action="store", type="string", dest="outputDir", help="output directory (default: output_rxxxx)", metavar="myOutputDir", default="output")  
parser.add_option("-v", "--verbose", action="store_true", dest="verbose", help="prints out the frame number as it is processed", default=False)
(options, args) = parser.parse_args()

########################################################
# Edit this variable accordingly
# Files are read for source_dir/runtag and
# written to write_dir/runtag.
# Be careful of the trailing "/"; 
# ensure you have the necessary read/write permissions.
########################################################
source_dir = "/reg/d/psdm/cxi/cxi25410/scratch/cleaned_hdf5/"
ang_avg_dir = "/reg/d/psdm/cxi/cxi25410/scratch/cleaned_hdf5/"

runtag = "r%s"%(options.runNumber)
write_dir = options.outputDir + '_' + runtag + '/'
canonicalOrderedHitsFN = write_dir + runtag + "_orderedFiles.txt"
canonicalOrderedIntensInHitsFN = write_dir + runtag + "_orderedIntens.txt"

h5files= []
integratedIntens = []
orderedH5Files = []
foundFiles = []
#Check if write_dir already exists, and which angavg files have already
#been copied over.
if not os.path.exists(write_dir):
	os.mkdir(write_dir)
else:
	write_anomaly_dir = write_dir 
	originaldir=os.getcwd()
	foundTypes=[write_dir]
	anomalousTypes=G.glob(write_anomaly_dir+"type[1-9]")
	if (len(anomalousTypes)>0):
		foundTypes+=anomalousTypes
	numTypes=len(foundTypes)
	for cDir in foundTypes:
		os.chdir(cDir)
		foundFiles+=G.glob("LCLS*angavg.h5")
		os.chdir(originaldir)
	print "Found %d types (including type0) with %d files" % (numTypes, len(foundFiles))	

if(options.verbose):
	print "detector at distance: %lf" % eDD.get_detector_dist_in_meters(runtag)

#Always check the angavg files in searchDir.
searchDir = ang_avg_dir + runtag
print "Now examining H5 files in %s/ ..."%(searchDir)   
searchstring="[a-zA-Z0-9\_]+"+runtag+"[a-z0-9\_]+-angavg.h5"
h5pattern = re.compile(searchstring)
h5files = [x for x in os.listdir(searchDir) if h5pattern.findall(x)]
numFiles = len(h5files)

#Check if file sorting has already been done.
if(os.path.exists(canonicalOrderedHitsFN) and os.path.exists(canonicalOrderedIntensInHitsFN)):
	print "Found sorted filenames and intens lists."
	f = open(canonicalOrderedHitsFN, 'r')
	orderedH5Files = N.array(f.read().split())
	f.close()
	f = open(canonicalOrderedIntensInHitsFN, 'r')
	integratedIntens = N.array([float(i) for i in f.read().split()])
	f.close()

#Check if presorted files are identical to those in searchDir:
#Re-sort only if different.	
if(len(set(orderedH5Files) - set(h5files)) == 0):
	print "Number of pre-sorted filenames matches those in %s" % (searchDir)
	P.plot(integratedIntens)
	P.title("Note the strong/weak hits cutoff")
	P.xlabel("sorted frame number")
	P.ylabel("integrated radial intensities of each run")
	P.show()
else:
	print "Found %d angavg files, sorting them now.."%(numFiles)
	integratedIntens = N.zeros(numFiles)
	for i in range(numFiles):
		fullFilePath = ang_avg_dir+runtag+"/"+h5files[i]
		f = H.File(fullFilePath, 'r')
		integratedIntens[i] = N.abs(N.array(f['/data/data'])).sum()
		f.close()
	ordering = integratedIntens.argsort()
	orderedH5Files = N.array(h5files)[ordering]
	P.plot(integratedIntens[ordering])
	P.title("Note the strong/weak hits cutoff")
	P.xlabel("sorted frame number")
	P.ylabel("integrated radial intensities of each run")
	P.show()
	orderedH5Files.tofile(canonicalOrderedHitsFN, sep="\n")
	(integratedIntens[ordering]).tofile(canonicalOrderedIntensInHitsFN, sep="\n")
	print "wrote to files %s and %s"%(canonicalOrderedIntensInHitsFN, canonicalOrderedHitsFN) 

colmax=1000
colmin=0
########################################################
# Imaging class modified from Ingrid Ofte's pyana_misc code
########################################################
class img_class (object):
	def __init__(self, inarr, inangavg , filename, meanWaveLengthInAngs=eDD.nominalWavelengthInAngs):
		self.inarr = inarr*(inarr>0)
		self.filename = filename
		self.inangavg = inangavg
		self.HIceQ ={}
		for i,j in eDD.iceHInvAngQ.iteritems():
			self.HIceQ[i] = eDD.get_pix_from_invAngsQ(runtag,j, meanWaveLengthInAngs)
		global colmax
		global colmin
	
	def on_keypress(self,event):
		global colmax
		global colmin
		if event.key == 'p':
			if not os.path.exists(write_dir):
				os.mkdir(write_dir)
			pngtag = write_dir +"%s.png" % (self.filename)	
			print "saving image as " + pngtag 
			P.savefig(pngtag)
		if event.key == 'r':
			colmin = self.inarr.min()
			colmax = self.inarr.max()
			P.clim(colmin, colmax)
			P.draw()
	
	def on_click(self, event):
		global colmax
		global colmin
		if event.inaxes:
			lims = self.axes.get_clim()
			colmin = lims[0]
			colmax = lims[1]
			range = colmax - colmin
			value = colmin + event.ydata * range
			if event.button is 1 :
				if value > colmin and value < colmax :
					colmin = value
			elif event.button is 2 :
				colmin = self.inarr.min()
				colmax = ((self.inarr)).max()
			elif event.button is 3 :
				if value > colmin and value < colmax:
					colmax = value
			P.clim(colmin, colmax)
			P.draw()
	
	def draw_img(self):
		global colmax
		global colmin
		fig = P.figure(num=None, figsize=(13.5, 5), dpi=100, facecolor='w', edgecolor='k')
		cid1 = fig.canvas.mpl_connect('key_press_event', self.on_keypress)
		cid2 = fig.canvas.mpl_connect('button_press_event', self.on_click)
		canvas = fig.add_subplot(121)
		canvas.set_title(self.filename)
		self.axes = P.imshow(self.inarr, origin='lower', vmax = colmax, vmin = colmin)
		self.colbar = P.colorbar(self.axes, pad=0.01)
		self.orglims = self.axes.get_clim()
		canvas = fig.add_subplot(122)
		canvas.set_title("angular average")
		P.plot(self.inangavg)
		maxAngAvg = (self.inangavg).max()
		numQLabels = len(self.HIceQ.keys())+1
		labelPosition = maxAngAvg/numQLabels
		for i,j in self.HIceQ.iteritems():
			P.axvline(j,0,colmax,color='r')
			P.text(j,labelPosition,str(i), rotation="45")
			labelPosition += maxAngAvg/numQLabels
		P.show() 

instructions="Right-click on colorbar to set maximum scale.\nLeft-click on colorbar to set minimum scale.\nCenter-click on colorbar (or press 'r') to reset color scale.\nInteractive controls for zooming at the bottom of figure screen (zooming..etc).\nPress 'p' to save PNG of image (with the current colorscales) in the appropriately named folder.\nHit Ctl-\ or close all windows (Alt-F4) to terminate viewing program."

#######################################################
#Records strong and weak file names into a textfile in output dir
######################################################
cutoff = int(input("weak/strong hits cutoff? "))

weakFiles = orderedH5Files[:cutoff]
recordTag = write_dir + runtag + "_weakAvgFiles.txt" 
N.array(weakFiles).tofile(recordTag, sep='\n')

strongFiles = orderedH5Files[cutoff:] 
recordTag = write_dir + runtag + "_strongAvgFiles.txt"  
N.array(strongFiles).tofile(recordTag, sep='\n')

#Copies only new files that were 
if(options.store_files):
	sStrong = set(strongFiles)
	sWeak = set(weakFiles)
	sFound = set(foundFiles)
	diffFiles = sStrong - sFound # != sFound-sStrong
	
	#Check we found angavg files in outputDir that weren't in searchDir
	if(len(sFound & set(h5files)) != len(sFound)):
		print "There are angavg files in the %s that are not in %s. Confused. Aborting." % (write_dir, searchDir)
		sys.exit(1)
	#Check if we found more(or fewer) angavg files than already in outputDir
	elif(len(diffFiles)==0):
		diffFiles = sFound - sStrong
		if(len(diffFiles) > 0):
			recordTag = write_dir + runtag + "_excessFiles.txt" 
			N.array(list(diffFiles)).tofile(recordTag, sep='\n')
			print "Excess angavg files in %s that did not make the cutoff stored in %s" % (write_dir, recordTag)
		else:
			print "All strong hits already in %s. Nothing copied." % (write_dir)
	else:
		print "%d strong hits copied earlier. Will copy only %d of the remainder" % (len(foundFiles), len(diffFiles))
		for fname in diffFiles: 
			diffractionName = source_dir+runtag+"/"+re.sub("-angavg",'',fname)
			angAvgName = ang_avg_dir + runtag + '/' + fname
			if(os.path.exists(diffractionName)):
				shutil.copyfile(angAvgName, write_dir+fname)

########################################################
# Weakly scattering files
#########################################################
waveLengths=[]
if (options.weakHitsTreatment == 1):
	print "averaging weak hits..."
	arr  = []
	avg = N.zeros(1233)
	fcounter = 0
	for fname in weakFiles:
		diffractionName = source_dir+runtag+"/"+re.sub("-angavg",'',fname)
		if(options.verbose and (fcounter%10==0)):
			print str(fcounter) + " of " + str(cutoff) + " weak files"
		fcounter += 1
		if(os.path.exists(diffractionName)):
			f = H.File(diffractionName, 'r')
			d = N.array(f['/data/data'])
			waveLengths.append(f['LCLS']['photon_wavelength_A'][0])
			if arr == []:
				arr = d.astype('float64')
			arr += d
			f.close()
			angAvgName = ang_avg_dir + runtag + '/' + fname
			f = H.File(angAvgName, 'r')
			davg = N.array(f['data']['data'][0])
			avg += davg
			f.close()
	
	arr /= len(weakFiles)
	[colmax, colmin] = [arr.max(), 0.]
	currImg = img_class(arr, avg/len(weakFiles), runtag+"_AvgPattWeakHits", meanWaveLengthInAngs=N.mean(waveLengths))
	print instructions
	currImg.draw_img()
	avg.tofile(write_dir + runtag + "_weak_avg.txt", sep = "\n", format="%lf")

########################################################
# Strongly scattering files
########################################################
waveLengths=[]
if(options.strongHitsTreatment == 1):
	print "averaging strong hits"
	arr  = []
	avg = N.zeros(1233)
	fcounter = 0
	numStrongFiles = len(strongFiles)
	numPresentStrongFiles = 0
	for fname in strongFiles:
		diffractionName = source_dir+runtag+"/"+re.sub("-angavg",'',fname)
		if(options.verbose and (fcounter%10==0)):
			print str(fcounter) + " of " + str(numStrongFiles) + " strong files"
		fcounter += 1
		if(os.path.exists(diffractionName)):
			f = H.File(diffractionName, 'r')
			d = N.array(f['/data/data'])
			waveLengths.append(f['LCLS']['photon_wavelength_A'][0])
			if arr == []: 
				arr = d.astype('float64')
			arr += d
			f.close()
			angAvgName = ang_avg_dir + runtag + '/' + fname
			f = H.File(angAvgName, 'r')
			davg = N.array(f['data']['data'][0])
			avg += davg
			f.close()
			numPresentStrongFiles += 1.
		else:
			print diffractionName + " not found!"
	
	arr /= numPresentStrongFiles
	[colmax, colmin] = [arr.max(), 0.] 
	currImg = img_class(arr, avg/len(strongFiles), runtag+"_AvgPattStrongHits", meanWaveLengthInAngs=N.mean(waveLengths))
	print instructions
	currImg.draw_img()
	
	avg.tofile(write_dir + runtag + "_strong_avg.txt", sep = "\n", format="%lf")
