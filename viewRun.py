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
parser.add_option("-s", "--step", action="store_true", dest="stepThroughImages", help="step through individual images in types (except type 0)", default=False)
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
viewTypes = N.zeros(numTypes+1, dtype='int64')
files = [""]*(numTypes+1)
write_anomaly_dir_types = [""]*(numTypes+1)
storeFlag = 0
write_anomaly_dir_types[0] = write_dir
os.chdir(write_dir)
files[0] = G.glob("LCLS*angavg.h5")
curr_dir = os.getcwd()
write_anomaly_dir_types[0] = curr_dir+'/'
os.chdir(originaldir)
viewTypes[0] = int(input("view "+write_dir+" (1 for yes, 0 for no)? "))

#write_anomaly_dir_types_jumbled = [write_dir]
for cDir in foundTypes:
	subString = re.sub(r""+write_dir+"type(\d+)",r"\1", cDir)
	storeFlag = int(subString)
	viewTypes[storeFlag] = int(input("view "+cDir+" (1 for yes, 0 for no)? "))
	os.chdir(cDir)
	files[storeFlag] = G.glob("LCLS*angavg.h5")
	curr_dir = os.getcwd()
	write_anomaly_dir_types[storeFlag] = curr_dir+'/'
	os.chdir(originaldir)

#	write_anomaly_dir_types_jumbled.append(i+"/")
#	write_anomaly_dir_types.append(i+"/")


#for cDir in write_anomaly_dir_types_jumbled:
#	if (cDir == write_dir):
#		storeFlag = 0
#	else:
#		subString = re.sub(r""+write_dir+"type(\d+)/",r"\1", cDir)
#		storeFlag = int(subString)
#	viewTypes[storeFlag] = int(input("view "+cDir+" (1 for yes, 0 for no)? "))
#	os.chdir(cDir)
#	files[storeFlag] = G.glob("LCLS*angavg.h5")
#	curr_dir = os.getcwd()
#	write_anomaly_dir_types[storeFlag] = curr_dir+'/'
#	os.chdir(originaldir)

colmax = options.maxIntens
colmin = 0
storeFlag = 0
########################################################
# Imaging class copied from Ingrid Ofte's pyana_misc code
########################################################
class img_class (object):
	def __init__(self, inarr, inangavg , filename, currTag, meanWaveLengthInAngs=eDD.nominalWavelengthInAngs, detectorDistance=eDD.get_detector_dist_in_meters(runtag)):
		self.origarr = inarr.copy()
		self.inarr = inarr*(inarr>0)
		self.filename = filename
		self.inangavg = inangavg
		self.wavelength = meanWaveLengthInAngs
		self.detectorDistance = detectorDistance
		self.HIceQ ={}
		global colmax
		global colmin
		global storeFlag
		self.tag = currTag
	
	def on_keypress_for_retagging(self,event):
		global colmax
		global colmin
		global storeFlag
		if event.key in [str(i) for i in range(numTypes+1)]:
			storeFlag = int(event.key)
			
			if(options.inspectOnly):
				print "Inspection only mode."
			else:
				if(not os.path.exists(write_anomaly_dir_types[storeFlag])):
					os.mkdir(write_anomaly_dir_types[storeFlag])
				pngtag = write_anomaly_dir_types[storeFlag] + "%s.png" % (self.filename)
				if(self.tag != 0):
					#delete previous assignment
					pngtag = write_anomaly_dir_types[self.tag] + "%s.png" % (self.filename)
					if os.path.isfile(pngtag):
						os.remove(pngtag)
						print "%s removed!" % (pngtag)
					else:
						print "No action taken."
					#Save new assignment if it's store flag not type 0
					if (storeFlag !=0):
						pngtag = write_anomaly_dir_types[storeFlag] + "%s.png" % (self.filename)
						P.savefig(pngtag)
						print "%s saved." % (pngtag)
						self.tag = storeFlag
				else:
					P.savefig(pngtag)
					print "%s saved." % (pngtag)
					self.tag = storeFlag
		if event.key == 'r':
			colmin = self.inarr.min()
			colmax = ((self.inarr<options.maxIntens)*self.inarr).max()
			P.clim(colmin, colmax)
			P.draw()
	
	def on_keypress_for_viewing(self,event):
		global colmax
		global colmin
		global storeFlag
		if event.key == 'p':
			pngtag = write_anomaly_dir_types[storeFlag] + "%s.png" % (self.filename)
			if(options.inspectOnly):
				print "Inspection only mode."
			else:
				P.savefig(pngtag)
				print "%s saved." % (pngtag)
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
		
	def draw_img_for_viewing(self):
		if(not options.inspectOnly):
			print "Press 'p' to save PNG."
		global colmax
		global colmin
		for i,j in eDD.iceHInvAngQ.iteritems():
			self.HIceQ[i] = eDD.get_pix_from_invAngsQ_and_detectorDist(runtag,j,self.detectorDistance, wavelengthInAngs=self.wavelength)

		fig = P.figure(num=None, figsize=(13.5, 6), dpi=100, facecolor='w', edgecolor='k')
		cid1 = fig.canvas.mpl_connect('key_press_event', self.on_keypress_for_viewing)
		cid2 = fig.canvas.mpl_connect('button_press_event', self.on_click)
		canvas = fig.add_axes([0.05,0.05,0.5,0.9])
		canvas.set_title(self.filename)
		self.axes = P.imshow(self.inarr, origin='lower', vmax = colmax, vmin = colmin)
		self.colbar = P.colorbar(self.axes, pad=0.01)
		self.orglims = self.axes.get_clim()


		(rTemp,cTemp) = self.inarr.shape
		#Approximate, center
		for i,j in self.HIceQ.iteritems():
			circ = P.Circle((rTemp/2, cTemp/2), radius=j)
			circ.set_fill(False)
			circ.set_edgecolor('k')
			canvas.add_patch(circ)

		canvas = fig.add_axes([0.6,0.05,0.35,0.4])
		canvas.set_title("angular average")
		maxAngAvg = (self.inangavg).max()

		numQLabels = len(self.HIceQ.keys())+1
		labelPosition = maxAngAvg/numQLabels
		for i,j in self.HIceQ.iteritems():
			P.axvline(j,0,colmax,color='r')
			P.text(j,labelPosition,str(i), rotation="45")
			labelPosition += maxAngAvg/numQLabels

		P.plot(self.inangavg)

		canvas = fig.add_axes([0.6,0.55,0.35,0.4], xlabel="pixel value", ylabel="count")
		canvas.set_title("histogram (omiting values<1) ")
		nonzeroarr = N.compress(self.origarr.flat>1, self.origarr.flat)
		canvas.hist(nonzeroarr, bins=N.arange(0,2000,10),log=True)

		P.show()
	
	def draw_img_for_retagging(self):
		if(not options.inspectOnly):
			print "Press 1-"+ str(numTypes)+ " to save png (overwrites old PNGs); Press 0 to undo (deletes png if wrongly saved)."
		global colmax
		global colmin
		for i,j in eDD.iceHInvAngQ.iteritems():
			self.HIceQ[i] = eDD.get_pix_from_invAngsQ_and_detectorDist(runtag,j,self.detectorDistance, wavelengthInAngs=self.wavelength)
		
		fig = P.figure(num=None, figsize=(13.5, 6), dpi=100, facecolor='w', edgecolor='k')
		cid1 = fig.canvas.mpl_connect('key_press_event', self.on_keypress_for_retagging)
		cid2 = fig.canvas.mpl_connect('button_press_event', self.on_click)
		canvas = fig.add_axes([0.05,0.05,0.5,0.9])
		canvas.set_title(self.filename)
		self.axes = P.imshow(self.inarr, origin='lower', vmax = colmax, vmin = colmin)
		self.colbar = P.colorbar(self.axes, pad=0.01)
		self.orglims = self.axes.get_clim()
		
		
		(rTemp,cTemp) = self.inarr.shape
		#Approximate, center
		for i,j in self.HIceQ.iteritems():
			circ = P.Circle((rTemp/2, cTemp/2), radius=j)
			circ.set_fill(False)
			circ.set_edgecolor('k')
			canvas.add_patch(circ)
		
		canvas = fig.add_axes([0.6,0.05,0.35,0.4])
		canvas.set_title("angular average")
		maxAngAvg = (self.inangavg).max()
		
		numQLabels = len(self.HIceQ.keys())+1
		labelPosition = maxAngAvg/numQLabels
		for i,j in self.HIceQ.iteritems():
			P.axvline(j,0,colmax,color='r')
			P.text(j,labelPosition,str(i), rotation="45")
			labelPosition += maxAngAvg/numQLabels
		
		P.plot(self.inangavg)
		
		canvas = fig.add_axes([0.6,0.55,0.35,0.4], xlabel="pixel value", ylabel="count")
		canvas.set_title("histogram (omiting values<1) ")
		nonzeroarr = N.compress(self.origarr.flat>1, self.origarr.flat)
		canvas.hist(nonzeroarr, bins=N.arange(0,2000,10),log=True)
		
		P.show()
	
	def draw_spectrum(self):
		print "Press 'p' to save PNG to correct type folder."
		global colmax
		global colmin
		localColMax=self.inarr.max()
		localColMin=self.inarr.min()
		aspectratio = 1.5*(self.inarr.shape[1])/(float(self.inarr.shape[0]))
		fig = P.figure(num=None, figsize=(13, 10), dpi=100, facecolor='w', edgecolor='k')
		cid1 = fig.canvas.mpl_connect('key_press_event', self.on_keypress_for_viewing) 
		cid2 = fig.canvas.mpl_connect('button_press_event', self.on_click)
		canvas = fig.add_axes([0.05,0.05,0.6,0.9], xlabel="q", ylabel="normalized angular average (will prompt to examine data larger than cutoff)")
		canvas.set_title(self.filename)
		self.axes = P.imshow(self.inarr, origin='lower', aspect=aspectratio, vmax = localColMax, vmin = localColMin)
		self.colbar = P.colorbar(self.axes, pad=0.01)
		self.orglims = self.axes.get_clim() 
		canvas2 = fig.add_axes([0.7,0.05,0.25,0.9], xlabel="log(sorting score)", ylabel="data")
		canvas2.set_ylim([0,numData])
		canvas2.plot(N.log(N.array(scoreKeeper)[ordering]),range(numData))
		P.show() 


avgArr = N.zeros((numTypes+1,1760,1760))
avgRadAvg = N.zeros((numTypes+1,1233))
typeOccurences = N.zeros(numTypes+1)
waveLengths={}
for i in range(numTypes+1):
	waveLengths[i] = []

#global storeFlag gets modified from retagging
#but currentlyExamining flag steps sequentially through types 
for currentlyExamining in range(numTypes+1):
	dirName = write_anomaly_dir_types[currentlyExamining]
	fcounter = 0
	numFilesInDir = len(files[currentlyExamining])
	if (viewTypes[currentlyExamining]==1):
		print "now examining "+ dirName
		os.chdir(dirName)
		for fname in files[currentlyExamining]:
			storeFlag = currentlyExamining
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
			if(options.stepThroughImages and (currentlyExamining !=0) and (options.inspectOnly == False)):
				currImg = img_class(d, davg, fname, currentlyExamining, meanWaveLengthInAngs=currWavelengthInAngs)
				currImg.draw_img_for_retagging()

				if(storeFlag != currentlyExamining):
					diffractionName = write_anomaly_dir_types[currentlyExamining]+"/"+re.sub("-angavg",'',fname)
					if(storeFlag != 0):
						print "moving angavg and pattern from " + str(currentlyExamining) + " to " + str(storeFlag)
						os.system("mv " + angAvgName + " " + write_anomaly_dir_types[storeFlag])
						os.system("mv " + diffractionName + " " + write_anomaly_dir_types[storeFlag])
					elif(storeFlag == 0):
						print "moving angavg from " + str(currentlyExamining) + " to " + str(storeFlag) + " , deleting pattern"
						os.system("mv " + angAvgName + " " + write_anomaly_dir_types[storeFlag])
						os.system("rm " + diffractionName)
					
			waveLengths[storeFlag].append(currWavelengthInAngs)
			avgArr[storeFlag] += d
			avgRadAvg[storeFlag] += davg
			typeOccurences[storeFlag] += 1
			fcounter += 1
		os.chdir(originaldir)

########################################################
# Loop to display all H5 files with ice anomalies. 
########################################################
print "Right-click on colorbar to set maximum scale."
print "Left-click on colorbar to set minimum scale."
print "Center-click on colorbar (or press 'r') to reset color scale."
print "Interactive controls for zooming at the bottom of figure screen (zooming..etc)."
print "Hit Ctl-\ or close all windows (Alt-F4) to terminate viewing program."

storeFlag=0
for dirName in write_anomaly_dir_types:
	if (typeOccurences[storeFlag] > 0.):
		avgArr[storeFlag] /= typeOccurences[storeFlag]
		avgRadAvg[storeFlag] /= typeOccurences[storeFlag]		
		if(storeFlag>0):
			typeTag = runtag+'_'+(dirName.split('/')[-2])
		else:
			typeTag = runtag+'_type0'
		currImg = img_class(avgArr[storeFlag], avgRadAvg[storeFlag], typeTag, storeFlag, meanWaveLengthInAngs=N.mean(waveLengths[storeFlag]))
		currImg.draw_img_for_viewing()
		if(options.inspectOnly):
			print "Inspection only mode."
		else:
			f = H.File(dirName +'/'+ typeTag + ".h5", "w")
			entry_1 = f.create_group("/data")
			entry_1.create_dataset("diffraction", data=avgArr[storeFlag])
			entry_1.create_dataset("angavg", data=avgRadAvg[storeFlag])	
			f.close()
	storeFlag += 1
