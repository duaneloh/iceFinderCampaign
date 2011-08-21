#!/usr/bin/env python

import numpy as N
import h5py as H
import matplotlib
import matplotlib.pyplot as P
import sys
import os
import re
from optparse import OptionParser

csvFileName="CXI25410.csv"
pixSize=110.E-6
nominalWavelengthInAngs=1.318
f = open(csvFileName, 'r')
lines = f.readlines()
f.close()
splitValues = [i.split(',') for i in lines]	
[tags, values] = [splitValues[0], splitValues[1:]]
expDescriptors = {}
runDescriptor="RUN"
detectorDistDescriptor="SAMPLE-DETECTOR Z (mm)"

for i in values:
	temp_dict = {}
	for j,k in zip(tags, i):
		temp_dict[j] = k
	curr_run = "r%04d"%(int(temp_dict[runDescriptor]))
	expDescriptors[curr_run] = temp_dict

def get_detector_dist_in_meters(run_tag):
	return (1.E-3)*float(expDescriptors[run_tag].get(detectorDistDescriptor))

def get_invAngsQ_from_pix(run_tag, pixPos):
	detDist=get_detector_dist_in_meters(run_tag)
	ang = N.arctan(pixPos*pixSize/detDist)
	return (2*N.pi/nominalWavelengthInAngs)*2.*N.sin(0.5*ang)

def get_pix_from_invAngsQ(run_tag, invAngsQ):
	detDist=get_detector_dist_in_meters(run_tag)
	temp = 2*N.arcsin(0.5*invAngsQ*nominalWavelengthInAngs/(2*N.pi))
	return detDist*N.tan(temp)/pixSize
	