'''
This code makes use of the scripts and part of the Jupyter notebooks from the cta-lstchain repository: https://github.com/cta-observatory/cta-lstchain
Copyright (c) 2018, the cta-lstchain developers All rights reserved.
'''
#Based on notebook from Rubén López-Coto


import pandas as pd
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


def mask_events(df,cuts):
	df_mask=[(df['leakage']<cuts['leakage']) & (df['intensity']>cuts['intensity'])& (df['intensity']<cuts['intensity_max']) & (df['wl']>cuts['wl']) & (df['gammaness']>cuts['gammaness']) & (df['reco_energy']<cuts['energy_up'])& (df['reco_energy']>cuts['energy_low'])]


	return(df_mask)


def calculate_theta2_par(df,df_mask):
	reco_src_x=df['reco_src_x'][df_mask[0]]
	reco_src_y=df['reco_src_y'][df_mask[0]]
	reco_alt=df.reco_alt[df_mask[0]]*180/3.14 #Convert to radian
	reco_az=np.arcsin(np.sin(df.reco_az[df_mask[0]]))*180/3.14
	theta2=np.power(reco_src_x*2,2)+np.power(reco_src_y*2,2)
	theta2=np.array(theta2)

	theta2_from_alt_az=np.power(reco_alt,2)+np.power(reco_az,2)
	theta2_from_alt_az=np.array(theta2_from_alt_az)
	return(theta2,theta2_from_alt_az)

def write_time_csv(filelist,cuts,dl2_params_lstcam_key):
	column_names = ["MJD Time", "Size", "Theta2",'Gammaness','Estimated energy (GeV)']
	df_final = pd.DataFrame(columns=column_names)
	for i in range(0,len(filelist)):
		df_i=pd.read_hdf(filelist[i],key=dl2_params_lstcam_key)
		df_i_mask=mask_events(df_i,cuts)
		theta2_i,theta2_from_altz_az_i=calculate_theta2_par(df_i,df_i_mask)
		fig,ax=plt.subplots(1,2)
		#fig.suptitle(str(os.listdir(args.directory)[i]))
		ax[0].hist(df_i['gammaness'].values,histtype='step')
		ax[0].set_xlabel('gammaness')
		ax[1].hist(theta2_i,histtype='step',bins=np.linspace(0,10))
		ax[1].set_xlabel('theta2_i')
		print(filelist[i])
		plt.figure()
		plt.hist(df_i.dragon_time[df_i_mask[0]]-df_i.ucts_time[df_i_mask[0]])
		plt.show()
		times=df_i.dragon_time[df_i_mask[0]].values
		t = Time(times,format='unix', scale='utc')
		data_i={'MJD Time': t.mjd,'Size':df_i.intensity[df_i_mask[0]].values,'Theta2':theta2_i,'Gammaness':df_i.gammaness[df_i_mask[0]].values,'Estimated energy (GeV)':1000*df_i.reco_energy[df_i_mask[0]].values}
		df_times=pd.DataFrame(data_i)
		df_times=df_times[df_times['Theta2']<0.4]
		df_final=pd.concat([df_final,df_times],names=column_names)

	df_final.to_csv('./../csv_files/Prueba_time_data.csv',index=None, header=True)
	print('DONE')



def obtain_info_theta2(filelist,cuts,dl2_params_lstcam_key):
	theta2_list=[]
	N_norm=0
	N=0
	norm_range_th2_min = 0.5
	norm_range_th2_max = 2
	theta2_cut=0.2
	for i in range(0,len(filelist)):
		df_i=pd.read_hdf(filelist[i],key=dl2_params_lstcam_key)
		df_i_mask=mask_events(df_i,cuts)
		theta2_i,theta2_from_altz_az_i=calculate_theta2_par(df_i,df_i_mask)
		N_norm += np.sum((theta2_i > norm_range_th2_min) & (theta2_i < norm_range_th2_max))
		N += np.sum(theta2_i < theta2_cut)
		theta2_list.append(theta2_i)
	theta2=np.concatenate(theta2_list)

	return(df_i,N_norm,N,theta2)

def calculate_theta2_sign(filelist,off_filelist,cuts,dl2_params_lstcam_key):

	df_i,Non_norm,Non,theta2=obtain_info_theta2(filelist,cuts,dl2_params_lstcam_key)
	df_i_off,Noff_norm,Noff,theta2_off=obtain_info_theta2(off_filelist,cuts,dl2_params_lstcam_key)	

	#Lets normalize
	Norm_theta2 = Non_norm / Noff_norm

	theta2_cut = 0.2
	Nex = Non - Noff * Norm_theta2
	S=np.sqrt(2)*(Non*np.log((1+Norm_theta2)/(Norm_theta2)*(Non/(Non+Noff)))+Noff*np.log((1+Norm_theta2)*(Noff/(Non+Noff))))**(1/2) #Formula 17 de Li and Ma (1983)

	#HISTOGRAM OF THETHA2
	nbins = 30
	range_max = 1
	obstime=158.52 #minutes
	fig, ax = plt.subplots()
	h_on = ax.hist(theta2,  label = 'ON data', bins=nbins, alpha=0.2, color = 'C1', range=[0,range_max])
	h_off = ax.hist(theta2_off, weights = Norm_theta2 * np.ones(len(theta2_off)),color='C2', range=[0,range_max], histtype='step', label = 'OFF data', bins=nbins, alpha=0.5)
	
	#Errors
	errorh=np.sqrt(h_on[0])
	errorh_off=np.sqrt(h_off[0])
	bins_centre=(h_on[1][1:]+h_on[1][:-1])/2
	bins_centre_off=(h_off[1][1:]+h_off[1][:-1])/2
	plt.errorbar(bins_centre, h_on[0], yerr=errorh, xerr=None,color='C1',fmt='.')
	plt.errorbar(bins_centre_off,h_off[0],yerr=errorh_off,xerr=None,color='C2',fmt='.')
	
	#Annotate and legends
	ax.annotate(s=f'Time={obstime:.0f} min \n Significance = {S:.2f} $\sigma$ \n Rate={Nex/obstime:.1f} $\gamma$/min',xy=(np.max(h_on[1]/3),8/9*max(h_on[0])),bbox=dict(boxstyle='Square',fc='w'),size=10,color='k')
	ax.set_xlabel(r'$\theta^2$ [deg$^2$]')
	ax.set_ylabel(r'Number of events')
	ax.set_ylim(2/5*min(h_off[0]),20/19*max(h_on[0]))
	ax.legend()
	plt.show()
	
def calculate_alpha_par(df,df_mask):
	cog_x = df.x[df_mask[0]]
	cog_y = df.y[df_mask[0]]
	psi = df.psi[df_mask[0]]
	hip=np.sqrt(np.power(cog_x,2)+np.power(cog_y,2))
	alpha = np.rad2deg(np.arccos((cog_x * np.cos(psi) + cog_y * np.sin(psi))/hip))
	alpha = alpha * (alpha < 90) + (180- alpha) * (alpha >= 90)

	return(alpha)


def obtain_info_alpha(filelist,cuts,dl2_params_lstcam_key):
	alpha_list=[]
	norm_range_alpha_min = 20.
	norm_range_alpha_max = 80.
	alpha_cut=8.
	N_norm=0
	N=0
	for i in range(0,len(filelist)):
		df_i=pd.read_hdf(filelist[i],key=dl2_params_lstcam_key)
		df_i_mask=mask_events(df_i,cuts)
		alpha_i=calculate_alpha_par(df_i,df_i_mask)
		N_norm += np.sum((alpha_i > norm_range_alpha_min) & (alpha_i < norm_range_alpha_max))
		N += np.sum(alpha_i < alpha_cut)
		alpha_list.append(alpha_i)
		
	alpha=np.concatenate(alpha_list)

	return(df_i,N_norm,N,alpha)

def obtain_energies(filelist,cuts,dl2_params_lstcam_key):
	energy_list=[]
	for i in range(0,len(filelist)):
                df_i=pd.read_hdf(filelist[i],key=dl2_params_lstcam_key)
                df_i_mask=mask_events(df_i,cuts)
                energy_i=df_i.reco_energy[df_i_mask[0]]
                energy_list.append(energy_i)
	energy=pd.concat(energy_list)
	return(energy)

def plot_energies(filelist,cuts,dl2_params_lstcam_key):
	energy=obtain_energies(filelist,cuts,dl2_params_lstcam_key)       
	
	#Distribution of energy
	fig,ax=plt.subplots()
	array=plt.hist(energy,bins=np.linspace(min(energy),800,30),label='Energy distribution',color='C5',alpha=0.1,density=True)
	plt.hist(energy,bins=np.linspace(min(energy),800,30),color='C5',alpha=0.5,histtype='step',density=True)
	h_e,bins=np.histogram(energy,bins=np.linspace(min(energy),800,30))

	#Errores
	h=array[0]
	edges=array[1]
	bin_centres=(edges[:-1]+edges[1:])/2
	error=np.sqrt(h_e)*h/h_e
	ax.errorbar(bin_centres,h,yerr=error,xerr=None,color='C5',fmt='.')

	#Define the pdf to fit to the distribution
	
	def landau(x,mu,c):
		return (1/(np.pi*c)*integrate.quad(lambda y:(np.exp(-y)*np.cos(y*((x-mu)/c)+2*y/np.pi*np.log(y/c))),0,np.inf)[0])
	landau=np.vectorize(landau,excluded=[2,3])

	#Do the fit
	fmax=9 #Maximum number to fit
	ppot,pcov=curve_fit(landau,bin_centres[0:fmax],h[0:fmax],p0=[160,40],maxfev=5000)
	plt.plot(np.arange(bin_centres[0],bin_centres[fmax]),landau(np.arange(bin_centres[0],bin_centres[fmax]),*ppot),'r-',linewidth=1,label='Landau fit',linestyle='dashed')
	ax.set_xlabel('Reconstructed energy(GeV)')
	ax.set_ylabel('a.u')

	#Identify the peak 
	valores=list(landau(np.arange(bin_centres[0],bin_centres[fmax]),*ppot))
	index=valores.index(max(valores))
	eth=np.arange(bin_centres[0],bin_centres[fmax])[index]

	#Annotation and legend
	plt.axvline(x=eth,linewidth=0.5,color='k',linestyle='dashed',label='Peak location')
	ax.annotate(s=f'Eth = ({eth:.0f} $\pm$ {np.sqrt(pcov[0,0]):.0f}) GeV',xy=(230,0.0033),size=10,color='k')
	plt.legend()
	plt.show()



def calculate_alpha_sign(filelist,off_filelist,cuts,dl2_params_lstcam_key):
	df_i,Non_norm,Non,alpha=obtain_info_alpha(filelist,cuts,dl2_params_lstcam_key)
	df_i_off,Noff_norm,Noff,alpha_off=obtain_info_alpha(off_filelist,cuts,dl2_params_lstcam_key)
		
	#Calculate the significance of the detection	
	obstime=158.52 #in minutes
	Norm_alpha = Non_norm / Noff_norm	
	Nex=Non-Noff*Norm_alpha
	S=np.sqrt(2)*(Non*np.log((1+Norm_alpha)/(Norm_alpha)*(Non/(Non+Noff)))+Noff*np.log((1+Norm_alpha)*(Noff/(Non+Noff))))**(1/2) #Formula 17 de Li and Ma (1983)
	nbins = 50

	#Plot the alpha histogram
	fig, ax = plt.subplots()
	h= ax.hist(alpha, label = 'ON data', bins=nbins, alpha=0.2, color = 'C1')
	hoff = ax.hist(alpha_off, weights = Norm_alpha * np.ones(len(alpha_off)),histtype='step', label = 'OFF data', bins=nbins, alpha=0.5, color = 'C2')
	#Errores		
	errorh=np.sqrt(h[0])
	errorh_off=np.sqrt(hoff[0])
	bins_centre=(h[1][1:]+h[1][:-1])/2
	bins_centre_off=(hoff[1][1:]+hoff[1][:-1])/2
	plt.errorbar(bins_centre, h[0], yerr=errorh, xerr=None,color='C1',fmt='.')
	plt.errorbar(bins_centre_off,hoff[0],yerr=errorh_off,xerr=None,color='C2',fmt='.')
	
	#Annotation and legend	
	ax.annotate(s=f'Time={obstime:.0f} min \n Significance = {S:.2f} $\sigma$ \n Rate={Nex/obstime:.1f} $\gamma$/min',xy=(np.max(h[1]/3),11/12*max(h[0])),bbox=dict(boxstyle='Square',fc='w'),size=10,color='k')
	ax.set_xlabel(r'$\alpha$ [deg]')
	ax.set_ylabel(r'Number of events')
	ax.set_ylim(4/5*min(hoff[0]),20/19*max(h[0]))
	ax.legend()
	plt.show()
	
	
cuts =  {
        "intensity": 500,
	"intensity_max":10000,
        "leakage": 0.2,
        "wl": 0.1,
        "gammaness":0.55,
        "energy_up":100000,
        "energy_low":0
}


def main():

	dl1_params_lstcam_key = 'dl1/event/telescope/parameters/LST_LSTCam'
	dl1_images_lstcam_key = 'dl1/event/telescope/image/LST_LSTCam'
	dl2_params_lstcam_key = 'dl2/event/telescope/parameters/LST_LSTCam'

	parser = argparse.ArgumentParser()
	parser.add_argument('--dir', '-d', action='store',type=str,
                    dest='directory',
                    default=None
                    )

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

	#calculate_theta2_sign(filelist,filelist_off,cuts,dl2_params_lstcam_key)
	#write_time_csv(filelist,cuts,dl2_params_lstcam_key)
	#plot_energies(filelist,cuts,dl2_params_lstcam_key)
	calculate_alpha_sign(filelist,filelist_off,cuts,dl2_params_lstcam_key)


if __name__ == "__main__":
    main()
