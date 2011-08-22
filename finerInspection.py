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
parser.add_option("-i", "--inspectDir", action="store", type="string", dest="inspectDirectory", help="data directory with the *angavg.h5 files you wish to inspect (defaults to output_rxxxx)", metavar="inputDir", default='')
parser.add_option("-o", "--outputDir", action="store", type="string", dest="outputDir", help="output directory will be appended by run number (default: output_rxxxx); separate types will be stored in output_rxxxx/anomaly/type[1-3]", default="output")
parser.add_option("-M", "--maxIntens", action="store", type="int", dest="maxIntens", help="doesn't plot intensities above this value", default=10000)
(options, args) = parser.parse_args()

source_dir = "/reg/d/psdm/cxi/cxi25410/scratch/cleaned_hdf5/"
ang_avg_dir = "/reg/d/psdm/cxi/cxi25410/scratch/cleaned_hdf5/"
runtag = "r%s"%(options.runNumber)
write_dir = options.outputDir + '_' + runtag + '/' 
write_anomaly_dir = write_dir + 'anomaly/'
inspectDirectory=""
if (options.inspectDirectory==''):
	inspectDirectory=write_dir
else:
	inspectDirectory=options.inspectDirectory
originaldir=os.getcwd()
os.chdir(inspectDirectory)
files = G.glob("*angavg.h5")

arr = []
print "reading ang_avgs.."
for i in files:
	f = H.File(i)
	arr.append(N.array(f['data']['data'][0]))
	f.close()

os.chdir(originaldir)

if(not os.path.exists(write_anomaly_dir)):
        os.mkdir(write_anomaly_dir)

masterArr = N.array(arr)
arr = masterArr

numData = len(arr)
angAvgLen = len(arr[0])

#Round 1: Sort by total intensities
scoreKeeper = N.zeros(numData)
scoreKeeper = [N.sum(i) for i in masterArr]
ordering = N.argsort(scoreKeeper)
sorted_arr = N.zeros((numData, angAvgLen))
files_round1 = []
filterLen = 5
medianFiltered_arr = N.zeros((numData, angAvgLen-filterLen))
print "sorting by total intensities.."
for i in range(numData):
	temp = masterArr[ ordering[i] ]
	max_temp = N.max(temp[530:560])
	#max_temp = N.max(temp[50:1153])
	min_temp = N.min(temp[50:1153]) #N.min(temp[700:800])
	sorted_arr[i] = (temp - min_temp) / (max_temp - min_temp) 
	for j in range(len(sorted_arr[i])-filterLen):
		medianFiltered_arr[i][j] = N.median(sorted_arr[i][j:j+filterLen])
	files_round1.append(files[ordering[i]])

round1FileNames = N.array([files[i] for i in ordering])

#Sort by maximum slope of median-filtered radial average
print "sorting by maximum of median filtered ang_avgs.."
sorted_arr2 = N.zeros((numData, angAvgLen))
scoreKeeper = [N.max(N.abs(i[201:1001]-i[200:1000])) for i in medianFiltered_arr]
ordering = N.argsort(scoreKeeper)
sorted_arr2 = sorted_arr[ordering]
round2FileNames = N.array([round1FileNames[i] for i in ordering])

colmax = options.maxIntens
colmin = 0
storeFlag = 0
numTypes = 5
write_anomaly_dir_types = [write_dir]
for i in range(numTypes):
	write_anomaly_dir_types.append(write_anomaly_dir+"type"+str(i+1)+"/") 
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
		if event.key in [str(i) for i in range(1,numTypes+1)]:
			storeFlag = int(event.key)
			recordtag = write_anomaly_dir_types[storeFlag] + runtag + "_" + event.key + ".txt"
			print "recording filename in " + recordtag
			storeFlag = int(event.key)
			if(not os.path.exists(write_anomaly_dir_types[storeFlag])):
				os.mkdir(write_anomaly_dir_types[storeFlag])
			f = open(recordtag, 'a+')
			#Could check here if it already has been recorded. Ouch!
			f.write(self.filename+"\n")
			f.close()
			pngtag = write_anomaly_dir_types[storeFlag] + "%s.png" % (self.filename)
			print "saving image as " + pngtag
			P.savefig(pngtag)
		if event.key == 'p':
			pngtag = write_anomaly_dir_types[storeFlag] + "%s.png" % (self.filename)
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
		self.axes = P.imshow(self.inarr, vmax = colmax, vmin = colmin)
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

	def draw_spectrum(self):
		global colmax
		global colmin
		localColMax=self.inarr.max()
		localColMin=self.inarr.min()
		aspectratio = (self.inarr.shape[1])/(self.inarr.shape[0])
		fig = P.figure(num=None, figsize=(8.5, 5), dpi=100, facecolor='w', edgecolor='k')
		cid1 = fig.canvas.mpl_connect('key_press_event', self.on_keypress) 
		cid2 = fig.canvas.mpl_connect('button_press_event', self.on_click)
		canvas = fig.add_subplot(111, xlabel="q", ylabel="normalized angular average")
		canvas.set_title(self.filename)
		self.axes = P.imshow(self.inarr, aspect=aspectratio, vmax = localColMax, vmin = localColMin)
		self.colbar = P.colorbar(self.axes, pad=0.01)
		self.orglims = self.axes.get_clim() 
		P.show()  	


print "Right-click on colorbar to set maximum scale."
print "Left-click on colorbar to set minimum scale."
print "Center-click on colorbar (or press 'r') to reset color scale."
print "Interactive controls for zooming at the bottom of figure screen (zooming..etc)."
print "Hit Ctl-\ or close all windows (Alt-F4) to terminate viewing program."

currImg = img_class(sorted_arr2,None, runtag+"_spectrum")
currImg.draw_spectrum()
cutoff = int(input("ice/water cutoff? "))

avgArr = N.zeros((numTypes+1,1760,1760))
avgRadAvg = N.zeros((numTypes+1,1233))
typeOccurences = N.zeros(numTypes+1)
waveLengths={}
for i in range(numTypes):
	waveLengths[i] = []
#Make average of nonanomalies
print "averaging water-only types.."
for fname in round2FileNames[:cutoff]:
	diffractionName = source_dir+runtag+"/"+re.sub("-angavg",'',fname)
	f = H.File(diffractionName, 'r')
	d = N.array(f['/data/data'])
	waveLengths[0].append(f['LCLS']['photon_wavelength_A'][0])
	f.close()
	avgArr[0] += d
	angAvgName = inspectDirectory + '/' + fname
	f = H.File(angAvgName, 'r')
	davg = N.array(f['data']['data'][0])
	f.close()
	avgRadAvg[0] += davg
	typeOccurences[0] += 1.

avgArr[0] /= typeOccurences[0]
avgRadAvg[0] /= typeOccurences[0]

storeFlag = 0
currImg = img_class(avgArr[0], avgRadAvg[0], runtag+"_AvgPattStrongWater",meanWaveLengthInAngs=N.mean(waveLengths[0]))
currImg.draw_img()

anomalies = round2FileNames[cutoff:]

########################################################
# Loop to display all H5 files with ice anomalies. 
########################################################
print "Right-click on colorbar to set maximum scale."
print "Left-click on colorbar to set minimum scale."
print "Center-click on colorbar (or press 'r') to reset color scale."
print "Interactive controls for zooming at the bottom of figure screen (zooming..etc)."
print "Press any single digit from '1-3' to save the H5 filename of current image to the appropriate file (e.g. r0079/r0079_1.txt), and also the taggedPNG ."
print "Hit Ctl-\ or close all windows (Alt-F4) to terminate viewing program."
waveLengths={}
for i in range(numTypes):
	waveLengths[i] = []
	
for fname in anomalies:
	storeFlag=0
	diffractionName = source_dir+runtag+"/"+re.sub("-angavg",'',fname)
	f = H.File(diffractionName, 'r')
	d = N.array(f['/data/data'])
	currWavelengthInAngs=f['LCLS']['photon_wavelength_A'][0]
	f.close()
	angAvgName = inspectDirectory + '/' + fname
	f = H.File(angAvgName, 'r')
	davg = N.array(f['data']['data'][0])
	f.close()
	currImg = img_class(d, davg, fname, meanWaveLengthInAngs=currWavelengthInAngs)
	currImg.draw_img()
	if(storeFlag == 1 or storeFlag == 2 or  storeFlag == 3):
		waveLengths[storeFlag].append(currWavelengthInAngs)
		avgArr[storeFlag] += d
		avgRadAvg[storeFlag] += davg
		typeOccurences[storeFlag] += 1
		if(not os.path.exists(write_anomaly_dir_types[storeFlag])):
			os.mkdir(write_anomaly_dir_types[storeFlag])

		print "mv " + angAvgName + " " + write_anomaly_dir_types[storeFlag]
		os.system("mv " + angAvgName + " " + write_anomaly_dir_types[storeFlag])
		os.system("cp " + diffractionName + " " + write_anomaly_dir_types[storeFlag])

for i in range(numTypes):
	if (typeOccurences[i] > 0.):
		storeFlag=i
		avgArr[i] /= typeOccurences[i]
		avgRadAvg[i] /= typeOccurences[i]
		typeTag = 'type'+str(i)+'_for_'+runtag
		currImg = img_class(avgArr[i], avgRadAvg[i], runtag+"_"+typeTag, meanWaveLengthInAngs=N.mean(waveLengths[i]))
		currImg.draw_img()
		(sx,sy) = avgArr[i].shape
		f = H.File(write_anomaly_dir_types[i] + typeTag + ".h5", "w")
		entry_1 = f.create_group("/data")
		entry_1.create_dataset("diffraction", data=avgArr[i])
		entry_1.create_dataset("angavg", data=avgRadAvg[i])	
		f.close()

"""
files_round2 = []
sorted_arr2 = N.zeros((numData-cutoff, angAvgLen))
scoreKeeper = [N.max(i[200:800]) for i in sorted_arr]
ordering = N.argsort(scoreKeeper)
for i in range(numData-cutoff):
	sorted_arr2[i] = sorted_arr[ ordering[i] ]
	files_round2.append(files_round1[ordering[i]])


psexportCmd = 'duaneloh@psexport.slac.stanford.edu:'
datadir='/reg/d/psdm/cxi/cxi25410/scratch/cleaned_hdf5/'
for i in files_round2:
	runnumSearch = re.search("(r[0-9]+)", i)
	runnum = runnumSearch.groups(1)[0]
	fileTag = re.sub("-angavg",r'',i)
	cmd = "scp " + psexportCmd + datadir + runnum + "/" + fileTag + " " + runnum + "/" + "h5files"
	os.system(cmd)
"""
