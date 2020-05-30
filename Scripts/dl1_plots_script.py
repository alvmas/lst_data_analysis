#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 11:28:58 2020

@author: alvarom
"""

import pandas as pd
import matplotlib
from lstchain.io.io import get_dataset_keys
from lstchain.io.io import dl1_params_lstcam_key
from lstchain.io.io import dl2_params_lstcam_key
import tables
import matplotlib.pyplot as plt
import numpy as np
#from fpdf import FPDF
import argparse
from ctapipe.instrument import CameraGeometry
from ctapipe.visualization import CameraDisplay

parser = argparse.ArgumentParser()
parser.add_argument('--file', '-f', action='store', type=str,
                    dest='filename',
                    default=None
                    )
args = parser.parse_args()

df= pd.read_hdf(args.filename,key=dl1_params_lstcam_key)



#We do some cuts on the data
mask=(df['width']/df['length']>0.1) & (df['leakage']<0.2) & (df['log_intensity']>2) & (df['log_intensity']<4)
bin_num=int(np.sqrt(len(df[mask]['log_intensity'])))

#We represent the dl1 data features
f1=plt.figure()
plt.hist(df[mask]['log_intensity'],bins=bin_num)
plt.savefig('log_intens_plot.png')

f2=plt.figure()
plt.hist(df[mask]['width'],bins=bin_num)
plt.savefig('width_plot.png')

f3=plt.figure()
plt.hist(df[mask]['length'],bins=bin_num)
plt.savefig('length_plot.png')

f4=plt.figure()
plt.hist(df[mask]['dragon_time'],bins=bin_num)
plt.savefig('time_plot.png')

f5=plt.figure()
plt.hist(df[mask]['psi'],bins=bin_num)
plt.savefig('psi_plot.png')

get_dataset_keys(args.filename)
group=tables.open_file(args.filename).root.dl1.event.telescope.image.LST_LSTCam
images=[x['image'] for x in group.iterrows()]


#Signal of each dl1 event on the camera
plt.figure()
for i in range(0,len(images)):
		if df['log_intensity'][i]>3: 
			geom = CameraGeometry.from_name("LSTCam-002")
			display=CameraDisplay(geom)
			display.image=images[i]
			display.add_colorbar()
			#display.set_limits_minmax()
			plt.show()

