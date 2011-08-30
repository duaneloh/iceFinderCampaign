#!/usr/bin/env python
import numpy as N
from numpy import linalg as LA
import h5py as H
import glob as G
import os 
import re
import pylab as P
from myModules import extractDetectorDist as eDD
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-r", "--run", action="store", type="string", dest="runNumber", help="run number you wish to view", metavar="rxxxx")
parser.add_option("-i", "--inspectOnly", action="store_true", dest="inspectOnly", help="inspect output directory", default=False)
parser.add_option("-o", "--outputDir", action="store", type="string", dest="outputDir", help="output directory (also inspection directory) will be appended by run number (default: output_rxxxx); separate types will be stored in output_rxxxx/anomaly/type[1-3]", default="output")
parser.add_option("-M", "--maxIntens", action="store", type="int", dest="maxIntens", help="doesn't plot intensities above this value", default=10000)
parser.add_option("-v", "--verbose", action="store_true", dest="verbose", help="prints out the frame number as it is processed", default=False)
(options, args) = parser.parse_args()

source_dir = "/reg/d/psdm/cxi/cxi25410/scratch/cleaned_hdf5/"
ang_avg_dir = "/reg/d/psdm/cxi/cxi25410/scratch/cleaned_hdf5/"
runtag = "r%s"%(options.runNumber)
write_dir = options.outputDir + '_' + runtag + '/' 
write_anomaly_dir = write_dir 
originaldir=os.getcwd()

foundTypes=G.glob(write_anomaly_dir+"type[1-9]")
numTypes=len(foundTypes)
print "Found %d types in addition to type0" % numTypes
write_anomaly_dir_types = [write_dir]
for i in foundTypes:
	write_anomaly_dir_types.append(i+"/")

viewTypes = N.zeros(numTypes+1)
files = []
storeFlag = 0
for cDir in write_anomaly_dir_types:
	viewTypes[storeFlag] = int(input("view "+cDir+" (1 for yes, 0 for no)? "))
	os.chdir(cDir)
	files.append(G.glob("LCLS*angavg.h5"))
	curr_dir = os.getcwd()
	write_anomaly_dir_types[storeFlag] = curr_dir
	os.chdir(originaldir)
	storeFlag += 1

colmax = options.maxIntens
colmin = 0
storeFlag = 0
########################################################
# Imaging class copied from Ingrid Ofte's pyana_misc code
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
		global storeFlag

	def on_keypress(self,event):
		global colmax
		global colmin
		global storeFlag
		if event.key == 'p':
			if(options.inspectOnly):
				print "Inspection only mode."
			else:
				pngtag = write_anomaly_dir_types[storeFlag] + "/%s.png" % (self.filename)
				print "saving image as " + pngtag
				P.savefig(pngtag)
		if event.key == 'r':
			colmin = self.inarr.min()
			colmax = ((self.inarr<options.maxIntens)*self.inarr).max()
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
				colmax = self.inarr.max()
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
		maxAngAvg = (self.inangavg).max()
		numQLabels = len(self.HIceQ.keys())+1
		labelPosition = maxAngAvg/numQLabels
		for i,j in self.HIceQ.iteritems():
			P.axvline(j,0,colmax,color='r')
			P.text(j,labelPosition,str(i), rotation="45")
			labelPosition += maxAngAvg/numQLabels

		P.plot(self.inangavg)
		P.show()


avgArr = N.zeros((numTypes+1,1760,1760))
avgRadAvg = N.zeros((numTypes+1,1233))
typeOccurences = N.zeros(numTypes+1)
waveLengths={}
for i in range(numTypes+1):
	waveLengths[i] = []

storeFlag=0
for dirName in write_anomaly_dir_types:
	fcounter = 0
	numFilesInDir = len(files[storeFlag])
	if (viewTypes[storeFlag]==1):
		print "now examining "+ dirName
		os.chdir(dirName)
		for fname in files[storeFlag]:
			if (options.verbose and (fcounter%10)==0):
				print "%d of %d files" % (fcounter, numFilesInDir)
			diffractionName = source_dir+runtag+"/"+re.sub("-angavg",'',fname)
			f = H.File(diffractionName, 'r')
			d = N.array(f['/data/data'])
			currWavelengthInAngs=f['LCLS']['photon_wavelength_A'][0]
			f.close()
			angAvgName = fname
			f = H.File(angAvgName, 'r')
			davg = N.array(f['data']['data'][0])
			f.close()
			waveLengths[storeFlag].append(currWavelengthInAngs)
			avgArr[storeFlag] += d
			avgRadAvg[storeFlag] += davg
			typeOccurences[storeFlag] += 1
			fcounter += 1
		os.chdir(originaldir)
	storeFlag += 1

########################################################
# Loop to display all H5 files with ice anomalies. 
########################################################
print "Right-click on colorbar to set maximum scale."
print "Left-click on colorbar to set minimum scale."
print "Center-click on colorbar (or press 'r') to reset color scale."
print "Interactive controls for zooming at the bottom of figure screen (zooming..etc)."
print "Press 'p' to save PNG in correct type-folder."
print "Hit Ctl-\ or close all windows (Alt-F4) to terminate viewing program."

storeFlag=0
for dirName in write_anomaly_dir_types:
	if (typeOccurences[storeFlag] > 0.):
		avgArr[storeFlag] /= typeOccurences[storeFlag]
		avgRadAvg[storeFlag] /= typeOccurences[storeFlag]		
		if(storeFlag>0):
			typeTag = runtag+'_'+(dirName.split('/')[-1])
		else:
			typeTag = runtag+'_type0'
		currImg = img_class(avgArr[storeFlag], avgRadAvg[storeFlag], typeTag, meanWaveLengthInAngs=N.mean(waveLengths[storeFlag]))
		currImg.draw_img()
		if(options.inspectOnly):
			print "Inspection only mode."
		else:
			f = H.File(dirName +'/'+ typeTag + ".h5", "w")
			entry_1 = f.create_group("/data")
			entry_1.create_dataset("diffraction", data=avgArr[storeFlag])
			entry_1.create_dataset("angavg", data=avgRadAvg[storeFlag])	
			f.close()
	storeFlag += 1
