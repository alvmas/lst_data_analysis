#Import the different modules
from ctapipe.io import event_source, EventSeeker
from ctapipe.utils import get_dataset_path
from ctapipe.visualization import CameraDisplay
from ctapipe.calib import CameraCalibrator
from ctapipe.instrument import CameraGeometry
from matplotlib import pyplot as plt
from astropy import units as u
from lstchain.calib.camera.calibrator import LSTCameraCalibrator
import numpy as np
from lstchain.calib.camera.r0 import LSTR0Corrections
import argparse
from traitlets.config.loader import Config
from ctapipe.image.extractor import LocalPeakWindowSum
from ctapipe.image.reducer import NullDataVolumeReducer
from lstchain.calib import lst_calibration
from ctapipe.image import hillas_parameters, tailcuts_clean
import lstchain.reco.utils as utils
from lstchain.reco.utils import get_event_pos_in_camera,disp,source_dx_dy
from lstchain.reco import dl0_to_dl1
from scipy.optimize import curve_fit
import matplotlib

#matplotlib.use('Qt5Agg')
config=Config({
	'LSTR0COrrections':{
		'pedestal_path':'/remote/gamma2/users/alvaromas/Pruebas/pedestal_drs4_1830_3.fits',
		'tel_id': 1,
		"r1_sample_start":None,
		"r1_sample_end":None
	}
})

parser = argparse.ArgumentParser()
parser.add_argument('--file', '-f', action='store', type=str,
                    dest='filename',
                    default=None
                    )
parser.add_argument('--pedestal', '-p', action='store', type=str,
                    dest='pedestal',
                    default=None
                    )
parser.add_argument('--cal', '-c', action='store', type=str,
                    dest='cal',
                    default=None
                    )

args = parser.parse_args()

def pedestal_correction(event,tel_id,lst_r0):
	lst_r0.subtract_pedestal(event,tel_id)
	event.r1.tel[tel_id].waveform=event.r1.tel[tel_id].waveform-lst_r0.offset


def integrate_pulse(event,tel_id):
        integrator=LocalPeakWindowSum()
        integration,peakpos=integrator(event.r1.tel[tel_id].waveform)
        return(integration,peakpos)

def calibrate_R0data(event,pe_cal,tel_id,lst_r0):
        #Firts, we make R0->R1
        pedestal_correction(event,tel_id,lst_r0)

        #We eliminate the edges where the signal seems saturated
        event.r1.tel[1].waveform[:,:,:2]=0
        event.r1.tel[1].waveform[:,:,38:]=0

        #We do the calibration
        integration,peakpos=integrate_pulse(event,1)

        #Convert to phe
        signals=integration.astype(float)
        signals=signals*pe_cal
        return(signals,peakpos)

def clean_data(image,cleaning_method,cleaning_parameters,geom):
        combined_clean=image.copy()
        signal_pixels=cleaning_method(geom,combined_clean,**cleaning_parameters)
        cleaned_image=combined_clean
        cleaned_image[~signal_pixels]=0
        return(cleaned_image)

def disp_param(src_pos,hillas):
        d=source_dx_dy(src_pos[0],src_pos[1],hillas.x,hillas.y)
        dis=disp.disp(hillas.x,hillas.y,src_pos[0],src_pos[1])
        uu=dis[4]*dis[2]*np.cos(hillas.psi)
        vv=dis[4]*dis[2]*np.sin(hillas.psi)
        rec_pos=disp.disp_to_pos(uu,vv,hillas.x,hillas.y)
        return(d,dis,uu,vv,rec_pos)


def main():
	#PEDESTAL
	#Let's visualize the 2D pixel distribution in a pedestal file for one event
	pedestal_source=event_source(args.pedestal,max_events=1000)
	seeker=EventSeeker(pedestal_source)

	#Calculate the pedestal RMS distribution
	#Define the windows to integrate the distribution in one pixel for each event
	window_start=12
	window_end=19

	integrated_ped=[]
	lst_r0=LSTR0Corrections(config=config)
	
	for event in pedestal_source:
		pedestal_correction(event,1,lst_r0)
		integrated_ped.append(np.sum(event.r1.tel[1].waveform[:,:,window_start:window_end],axis=2))

	integrated_ped=np.array(integrated_ped)
	pedestal_var=np.var(integrated_ped,axis=0)
	pedestal_mean=np.mean(integrated_ped,axis=0)

	#signals,peakpos=calibrate_R0data(event,pe_cal,1,lst_r0)
	
	def gaus(x,a,x0,sigma):
		return a*np.exp(-(x-x0)**2/(2*sigma**2))
	'''
	plt.figure()	
	plt.subplot(1,2,1)
	geom = CameraGeometry.from_name("LSTCam-002")
	display=CameraDisplay(geom)
	display.image=signals[0]
	display.add_colorbar()
	plt.tight_layout()
	display.axes.text(2.0, 0, 'ADC Counts', rotation=90)
	plt.title('Pedestal signal (HG)')

	plt.subplot(1,2,2)
	display=CameraDisplay(geom)
	display.image=signals[1]
	display.add_colorbar()
	display.axes.text(2.0, 0, 'ADC Counts', rotation=90)
	plt.title('Pedestal signal (LG)')
	plt.tight_layout()
	plt.show()
	plt.savefig('ped_image.png')
	'''
	#Plots, one for each channel (HG and LG) with the RMS dsitribution
	plt.figure()
	array=plt.hist(integrated_ped[:,0,0],bins=50,alpha=0.4,label='HG Channel RMS distribution')
	#h,edges=np.histogram(integrated_ped[:,0,0],bins=50)
	h=array[0]
	edges=array[1]
	bin_centres=(edges[:-1]+edges[1:])/2

	mean=sum(h*bin_centres)/sum(h)
	sigma=sum(h*(bin_centres-mean)**2)/sum(h)	
	ppot,pcov=curve_fit(gaus,bin_centres,h,p0=[1,2800,200],maxfev=5000)
	plt.plot(np.arange(min(edges),max(edges)),gaus(np.arange(min(edges),max(edges)),*ppot),'r-',label='Gaussian fit')
	#plt.hist(integrated_ped[:,1,0],bins=50,alpha=0.4,label='LG Channel')
	plt.xlabel('ADC counts')
	plt.legend()
	plt.xlim(min(edges),2950)
	plt.annotate(s=f'Gauusian parameters: \n Mean = {ppot[1]:.2f} \n Standard deviation:{ppot[2]:.2f}', xy=[2870, 70], color = 'k')	
	plt.show()
	plt.savefig('dist_ped.png')

	#FLAT FIELD EVENTS
	cal_source=event_source(args.cal,max_events=1000)

	#We make the calibration
	charge_1=[]
	charge_2=[]
	num_events=0
	threshold=1000

	for event in cal_source:
		if event.count == 0:  # We don't use the first event.
        		continue
		#Firts, we make R0->R1
		pedestal_correction(event,1,lst_r0)

		#We eliminate the edges where the signal seems saturated
		event.r1.tel[1].waveform[:,:,:3]=0
		event.r1.tel[1].waveform[:,:,37:]=0
	
		integration,peakpos=integrate_pulse(event,1)
	
		if np.mean(integration[0])>threshold:
			charge_1.append(integration[0])
			charge_2.append(integration[1])
			num_events+=1

	channel1_charge=np.zeros((num_events,1855))
	channel2_charge=np.zeros((num_events,1855))

	channel1_charge=np.array(charge_1)
	channel2_charge=np.array(charge_2)


	#Let's see how the integration of the peak is done for the last  event
	t=np.arange(2,38,1)
	pixel=12
	fig8,ax=plt.subplots(1,2)
	ax[0].step(t,event.r1.tel[1].waveform[0,pixel,2:38])
	ax[1].step(t,event.r1.tel[1].waveform[1,pixel,2:38])
	ax[0].axvline(peakpos[0][pixel],linestyle='--',color='green')
	ax[1].axvline(peakpos[1][pixel],linestyle='--',color='green')
	ax[0].title.set_text('HG FF signal and peak for one pixel')
	ax[1].title.set_text('LG FF signal and peak for one pixel')
	ax[0].set_xlabel('t(ns)')
	ax[1].set_xlabel('t(ns)')


	#Distribution of charge in a pixel
	fig9,ax=plt.subplots(1,2)
	ax[0].hist(channel1_charge[:,pixel],bins=50)
	ax[1].hist(channel2_charge[:,pixel],bins=50)
	ax[0].title.set_text('Distribution of charge in one pixel (HG)')
	ax[1].title.set_text('Distribution of charge in one pixel (LG)')
	ax[0].set_xlabel('FADC')
	ax[1].set_xlabel('FADC')
	ax[0].set_ylabel('Counts')
	ax[1].set_ylabel('Counts')
	
	final_charge=np.zeros((2,num_events,1855))
	final_charge[0]=channel1_charge
	final_charge[1]=channel2_charge
	plt.show()
	
	#Let's estimate the number o phe using excess noise factor method
	mean_charge=np.mean(final_charge,axis=1)
	var_charge=np.var(final_charge,axis=1)
	F_square=1.2 
	phe=mean_charge*mean_charge/(var_charge-pedestal_var)*F_square
	pe_cal=phe/mean_charge
	
	
	fig10,ax=plt.subplots(1,2)
	ax[0].title.set_text('Distribution of mean pe (HG)')
	ax[1].title.set_text('Distribution of mean pe (LG)')
	
	ax[0].hist(phe[0])
	ax[1].hist(phe[1],bins=np.arange(0,120,5))
	ax[0].set_xlabel('pe')
	ax[1].set_xlabel('pe')
	ax[0].set_ylabel('Counts')
	ax[1].set_ylabel('Counts')
		

	fig11,ax=plt.subplots(1,2)
	
	plt.subplot(1,2,1)
	geom = CameraGeometry.from_name("LSTCam-002")
	display=CameraDisplay(geom)
	display.image=phe[0]
	display.add_colorbar()
	plt.tight_layout()	
	display.axes.text(2.0, 0, 'pe', rotation=90)
	plt.title('FF Calibrated signal (HG)')
	
	plt.subplot(1,2,2)
	display=CameraDisplay(geom)
	display.image=phe[1]
	display.add_colorbar()
	display.axes.text(2.0, 0, 'pe', rotation=90)
	plt.title('FF Calibrated signal (LG)')
	plt.tight_layout()
	plt.show()

	plt.savefig('ff.png')
	
	#REPRESENT PEDESTAL EVENTS
	integrated_ped=[]
	for event in pedestal_source:
		pedestal_correction(event,1,lst_r0)
		integrated_ped.append(np.sum(event.r1.tel[1].waveform[:,:,window_start:window_end],axis=2))

	integrated_ped=np.array(integrated_ped)
	pedestal_var=np.var(integrated_ped,axis=0)
	pedestal_mean=np.mean(integrated_ped,axis=0)

	signals,peakpos=calibrate_R0data(event,pe_cal,1,lst_r0)
	plt.figure()
	plt.subplot(1,2,1)
	geom = CameraGeometry.from_name("LSTCam-002")
	display=CameraDisplay(geom)
	display.image=signals[0]
	display.add_colorbar()
	plt.tight_layout()
	display.axes.text(2.0, 0, 'pe', rotation=90)
	plt.title('Pedestal signal (HG)')

	plt.subplot(1,2,2)
	display=CameraDisplay(geom)
	display.image=signals[1]
	display.add_colorbar()
	display.axes.text(2.0, 0, 'pe', rotation=90)
	plt.title('Pedestal signal (LG)')
	plt.tight_layout()
	plt.show()
	plt.savefig('ped_image.png')

	#RAW DATA->DL1
	source=event_source((args.filename),max_events=1000) #We read the raw data
	seeker=EventSeeker(source) #This is used to look for an event
	prob_event=seeker[49] 
	prob_event_data=prob_event.r0.tel[1] 
		
	#We choose one pixel to see its signal and see the pulse
	signal=prob_event_data.waveform[0][308] #Pixel 700
	fig3=plt.figure()
	plt.step(np.linspace(2,len(signal),len(signal)-2),signal[2:],label='Pixel with intense pulse')
	plt.xlabel('Sample number')
	
	signal=prob_event_data.waveform[0][100] #Pixel 257
	plt.step(np.linspace(2,len(signal),len(signal)-2),signal[2:],label='Pixel with no intense pulse')
	plt.xlabel('Sample number')
	plt.title('Signal in one pixel')
	#plt.show()
	plt.legend()
	
	#We do the correction for raw data
	signals,peakpos=calibrate_R0data(prob_event,pe_cal,1,lst_r0)
	#fig8,ax=plt.subplots(1,2)
	plt.figure()
	t=np.arange(2,38,1)
	plt.step(t,prob_event.r1.tel[1].waveform[0,308,2:38],label='Signal with subtracted baseline')
	plt.fill_between(t[int(peakpos[0][308])-5:int(peakpos[0][308])+2],prob_event.r1.tel[1].waveform[0,308,int(peakpos[0][308])-3:int(peakpos[0][308])+4],step='pre',alpha=0.4,color='orange',label='Window of integration')
	#ax[1].step(t,prob_event.r1.tel[1].waveform[1,685,2:38])
	plt.axvline(peakpos[0][308],linestyle='--',color='red',label='Peak position')
	#ax[1].axvline(peakpos[1][685],linestyle='--',color='red')
	plt.title('Signal and peak for one pixel')
	#ax[1].title.set_text('LG FF signal and peak for one pixel')
	plt.xlabel('t(ns)')
	plt.legend()
	plt.tight_layout()
	#ax[1].set_xlabel('t(ns)')

	
	#Choose the correct caibration and correct possible saturations and combine
	combined=signals[0].copy() #We will only change to the low gain channel if the pixel is saturated
	
	for pixel in range(0,combined.size):
		if np.any(prob_event_data.waveform[0][pixel]>4094):
			combined[pixel]=signals[1][pixel]
			
	
	#Cleaning
	prob_pixels=[1533,1535,1536,1534]
	for i in prob_pixels:
		combined[i]=0

	cleaning_parameters={'boundary_thresh':3,
	                'picture_thresh':6,
	                'keep_isolated_pixels':True,
	                'min_number_picture_neighbors':2
          }
	cleaned_image=clean_data(combined,tailcuts_clean,cleaning_parameters,geom)
	
		
	cleaning_parameters_less={'boundary_thresh':2,
       	       		'picture_thresh':4,
        	        'keep_isolated_pixels':False,
                	'min_number_picture_neighbors':1
                	}
	cleaned_image_less=clean_data(combined,tailcuts_clean,cleaning_parameters_less,geom)

	cleaning_parameters_less_less={'boundary_thresh':1,
			'picture_thresh':3,
			'keep_isolated_pixels':False,
			'min_number_picture_neighbors':1
			}
	cleaned_image_less_less=clean_data(combined,tailcuts_clean,cleaning_parameters_less_less,geom)	
	
	#We make two pots, one for HG and one for LG to see the shape of the signal once we have apploed the R0->R1 calibration.
	#fig11,ax=plt.subplots(1,2)
	plt.figure()
	plt.step(t,prob_event.r0.tel[1].waveform[0,685,2:38],label='R0')
	plt.step(t,prob_event.r1.tel[1].waveform[0,685,2:38],label='R1')#Pixel number 0
	#ax[1].step(t,prob_event.r0.tel[1].waveform[1,685,2:38],label='R0')
	#ax[1].step(t,prob_event.r1.tel[1].waveform[1,685,2:38],label='R1')
	plt.title('HG Integrated pedestal')
	#ax[1].title.set_text('LG Integrated pedestal')
	plt.legend()
	#ax[1].legend()
	plt.xlabel('t(ns)')
	#ax[1].set_xlabel('t(ns)')
	plt.ylabel('ADC counts')
	#ax[1].set_ylabel('ADC counts')
	

	#Visualize R0 data for HG and LG
	fig,ax=plt.subplots(1,2)
	
	plt.subplot(1,2,1)
	geom = CameraGeometry.from_name("LSTCam-002")
	display=CameraDisplay(geom)
	display.image=prob_event.r0.tel[1].waveform[0][:,0]
	display.add_colorbar()
	display.axes.text(2.0, 0, 'FADC', rotation=90)
	plt.title('R0 signal (HG)')
		
	plt.subplot(1,2,2)
	geom = CameraGeometry.from_name("LSTCam-002")
	display=CameraDisplay(geom)
	display.image=prob_event.r0.tel[1].waveform[1][:,0]
	display.add_colorbar()
	display.axes.text(2.0, 0, 'FADC', rotation=90)
	plt.title('R0 signal (LG)')
	plt.show()
	
	#Visualize calibrated  data for LG and HG
	plt.figure()
	geom = CameraGeometry.from_name("LSTCam-002")
	display=CameraDisplay(geom)
	display.image=signals[0]
	display.add_colorbar()
	display.axes.text(2.0, 0, 'pe', rotation=90)
	plt.title('Calibrated data (HG)')
	
	plt.figure()
	geom = CameraGeometry.from_name("LSTCam-002")
	display=CameraDisplay(geom)
	display.image=signals[1]
	display.add_colorbar()
	plt.title('Calibrated data (LG)')
	display.axes.text(2.0, 0, 'pe', rotation=90)

	#Visualize calibrated combined data
	plt.figure()
	display=CameraDisplay(geom)
	display.image=combined
	display.add_colorbar()
	plt.title('Calibrated combined data')
	display.axes.text(2.0, 0, 'pe', rotation=90)
	plt.show()
	
	
	#Visualize cleaned data
	
	plt.figure()
	display=CameraDisplay(geom)
	display.edgecolor='black'
	display.cmap=plt.cm.magma_r
	display.image=cleaned_image
	display.add_colorbar()
	display.axes.text(2.0, 0, 'pe', rotation=90)
	plt.title('Cleaned image')
	plt.savefig('cleaned_fig_iso')
	j=0
	k=0
	for i in range(0,len(cleaned_image)):
		if cleaned_image[i]>k:
			k=cleaned_image[i]
			j=i
			print(j)
	plt.figure()
	display=CameraDisplay(geom)
	display.cmap=plt.cm.magma_r
	display.image=cleaned_image_less
	display.add_colorbar()
	display.axes.text(2.0, 0, 'pe', rotation=90)
	plt.title('Cleaned image')

	plt.savefig('cleaned_less')

	plt.figure()
	display=CameraDisplay(geom)
	display.cmap=plt.cm.magma_r
	display.image=cleaned_image_less_less
	display.add_colorbar()
	display.axes.text(2.0, 0, 'pe', rotation=90)
	plt.title('Cleaned image')
	plt.savefig('cleaned_less_less')
	plt.show()

	#Calculate hillas parameters
	h=hillas_parameters(geom,cleaned_image)
	#src_pos=get_event_pos_in_camera(prob_event,tel)
	src_pos=[0,0]
	d,dis,uu,vv,rec_pos=disp_param(src_pos,h)
	
	plt.figure()
	display=CameraDisplay(geom)
	#display.cmap=plt.cm.twilight
	display.image=cleaned_image
	display.add_colorbar()
	display.axes.text(2.0, 0, 'pe', rotation=90)
	plt.title('Cleaned image')
	display.overlay_moments(h,color='white',linewidth=2)
	plt.scatter(src_pos[0],src_pos[1],color='black',label='actual source position')
	plt.quiver(h.x,h.y,d[0].value,d[1].value,units='xy',scale=1,label='actual disp')
	plt.quiver(h.x,h.y,uu.value,vv.value,units='xy',color='red',scale=1,label='reconstructed disp')
	plt.scatter(rec_pos[0].value,rec_pos[1].value,color='red',label='reconstructed source position')
	plt.legend()
	plt.show()
	plt.savefig('disp.png')

	fig11,ax=plt.subplots(1,3)
	plt.subplot(1,3,2)
	geom = CameraGeometry.from_name("LSTCam-002")
	display=CameraDisplay(geom)
	display.image=combined
	plt.title('Calibrated signal')

	plt.subplot(1,3,1)
	display=CameraDisplay(geom)
	display.image=prob_event.r0.tel[1].waveform[0][:,0]
	plt.title('Raw Event (R0)')
	

	plt.subplot(1,3,3)
	display=CameraDisplay(geom)
	display.cmap=plt.cm.magma_r
	display.image=cleaned_image
	plt.title('Cleaned image')
	plt.tight_layout()
	plt.savefig('conjunto.png')
	plt.show()

if __name__== "__main__":
	main()
