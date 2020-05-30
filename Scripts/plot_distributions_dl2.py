import pandas as pd
#from lstchain.tests.test_lstchain import dl2_params_lstcam_key
#from lstchain.io.io import dl1_params_lstcam_key
import argparse
import matplotlib.pyplot as plt
import csv
from astropy.time import Time
import numpy as np
from scipy.optimize import curve_fit
import matplotlib
import scipy.integrate as integrate
import os
from astropy.time import Time





def plot_hists(filelist,filelist_off,dl2_params_lstcam_key):
	fig,ax=plt.subplots(1,2)
	df=pd.read_hdf(filelist[0],key=dl2_params_lstcam_key)
	df_off=pd.read_hdf(filelist_off[0],key=dl2_params_lstcam_key)
	ax[0].hist(df['gammaness'],histtype='step')
	ax[0].set_title('ON events')
	ax[0].set_xlabel('gammaness')
	ax[0].set_ylabel('a.u')
	ax[1].hist(df_off['gammaness'],histtype='step')	
	ax[1].set_title('OFF events')
	ax[1].set_xlabel('gammaness')
	ax[1].set_xlabel('a.u')

	fig,ax=plt.subplots(1,2)
	ax[0].hist(df['length'],histtype='step')
	ax[0].set_title('ON events')
	ax[0].set_xlabel('length')
	ax[0].set_ylabel('a.u')
	ax[1].hist(df_off['length'],histtype='step') 
	ax[1].set_title('OFF events')
	ax[1].set_xlabel('length')
	ax[1].set_ylabel('a.u')



	fig,ax=plt.subplots(1,2)
	ax[0].hist(df['width'],histtype='step')
	ax[0].set_title('ON events')
	ax[0].set_xlabel('width')
	ax[0].set_xlabel('a.u')
	ax[1].hist(df_off['width'],histtype='step')
	ax[1].set_title('OFF events')
	ax[1].set_xlabel('width')
	ax[1].set_ylabel('a.u')


	fig,ax=plt.subplots(1,2)
	ax[0].hist(df['leakage'],histtype='step')
	ax[0].set_title('ON events')
	ax[0].set_xlabel('leakage')
	ax[0].set_ylabel('a.u')
	ax[1].hist(df_off['leakage'],histtype='step')
	ax[1].set_title('OFF events')
	ax[1].set_xlabel('leakage')
	ax[1].set_ylabel('a.u')


	fig,ax=plt.subplots(1,2)
	ax[0].hist(df['intensity'],histtype='step')
	ax[0].set_title('ON events')
	ax[0].set_xlabel('intensity')
	ax[0].set_ylabel('a.u')
	ax[0].set_xlim(0,10000)
	#ax[0].set_xscale('log')
	ax[1].hist(df_off['intensity'],histtype='step')
	ax[1].set_title('OFF events')
	ax[1].set_xlabel('intensity')
	ax[1].set_ylabel('a.u')
	ax[1].set_xlim(0,10000)
	#ax[1].set_xscale('log')
	plt.show()

def main():

	dl1_params_lstcam_key = 'dl1/event/telescope/parameters/LST_LSTCam'
	dl1_images_lstcam_key = 'dl1/event/telescope/image/LST_LSTCam'
	dl2_params_lstcam_key = 'dl2/event/telescope/parameters/LST_LSTCam'

	parser = argparse.ArgumentParser()
	parser.add_argument('--dir', '-d', action='store',type=str,dest='directory',default=None)

	args = parser.parse_args()
	filelist=[]
	filelist_off=[]
	for x in os.listdir(args.directory):
		rel_dir = os.path.relpath(args.directory)
		rel_file = os.path.join(rel_dir, x)
		if 'ON' in x:
			filelist.append(rel_file)
		elif 'OFF' in x:
			filelist_off.append(rel_file)
	plot_hists(filelist,filelist_off,dl2_params_lstcam_key)
        
if __name__ == "__main__":
    main()
                      
