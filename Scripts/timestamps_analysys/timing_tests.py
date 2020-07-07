#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 10:02:25 2019

@author: alvarom
"""


import time
import numpy as np
import os
import matplotlib.pyplot as plt
from sympy import*
import math
from pylab import *
from scipy.stats import linregress
from matplotlib.backends.backend_pdf import PdfPages
from fpdf import FPDF

def func(x, a, b,d):
    return a*np.exp(b*x)+d

def read_file(filename):
    event_dic={}
    os.chdir("/home/pablo/Escritorio/AlvaroM_Pruebas") #PATH OF THE DIRECTORY OF THE FILE!!!!!!!!!
    with open(str(filename), 'r') as inputfileok:
        file_inok = inputfileok.readlines()
        bch_counter = np.zeros(len(file_inok)-2)
        event_counter = np.zeros(len(file_inok)-2) 
        event_busy_counter = np.zeros(len(file_inok)-2) 
        timestamp = np.zeros(len(file_inok)-2)
        pps = np.zeros(len(file_inok)-2) 
        system_time= np.zeros(len(file_inok)-2) 
        for i in range(2,len(file_inok)):
            linesstrip = file_inok[i].strip()
            linessplit = linesstrip.split('\t')   
            bch_counter[i-2]=linessplit[0]
            event_counter[i-2]=linessplit[1]
            event_busy_counter[i-2]=linessplit[2]
            timestamp[i-2]=linessplit[3]
            pps[i-2]=linessplit[4]
            system_time[i-2]=linessplit[-1]
    inputfileok.close()
    
    event_dic['bch_counter']=bch_counter
    event_dic['event_counter']=event_counter 
    event_dic['event_busy_counter']=event_busy_counter 
    event_dic['timestamp']= timestamp
    event_dic['pps']=pps
    event_dic['system_time']=system_time
    
    return event_dic


def Infovstime(event_dic):
    relative_time=(event_dic['timestamp']-event_dic['timestamp'][0])*10**(-9)
    
    plt.figure()
    fig=plt.plot(event_dic['bch_counter'],relative_time,'.')
    ax = plt.gca()
    plt.ylabel('T-T$_{0}$ (s)')
    plt.xlabel('Bunch counter')
    plt.title(' TimeStamp vs Bunch Counter')
    regress_par=linregress(event_dic['bch_counter'],relative_time)
    fit_fn = np.poly1d(regress_par[0:2])
    plt.plot(event_dic['bch_counter'],fit_fn(event_dic['bch_counter']),'r')
    text(0.7, 0.1, 'T$_{0}$='+str(event_dic['timestamp'][0]), horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
    text(0.25, 0.65, 'Linear Regression Parameters:'+'\n'+'r$^{2}$='+str("{0:.3f}".format(regress_par[2]))+'\n'
         +'Slope:'+str("{0:.3f}".format(regress_par[0]))+'s'+'\n'+'Intercept:'+str("{0:.2f}".format(regress_par[1]))+'s', horizontalalignment='center',
         verticalalignment='center', transform=ax.transAxes,bbox=dict(facecolor='white', alpha=0.5))
    
    
    message=0
    s_points=[]
    s_times=[]
    for i in range(0,len(event_dic['bch_counter'])-1):
        if (event_dic['bch_counter'][i]>event_dic['bch_counter'][i+1]) or (event_dic['bch_counter'][i+1]-event_dic['bch_counter'][i])>1:
            s_times.append(relative_time[i+1])
            s_points.append(event_dic['bch_counter'][i+1])
            message='Suspicious points detected'
    
    plt.plot(s_points,s_times,'k.')  
    
    if message==0:
        message='No suspicious points detected'
        plt.legend(['Bunch-Time points', 'Linear Fit'])
    else:
         plt.legend(['Bunch-Time points', 'Linear Fit','Suspicious points'])
        
         
    text(0.7, 0.3, str(message), horizontalalignment='center',
         verticalalignment='center', transform=ax.transAxes,bbox=dict(facecolor='white', alpha=0.5))
   
    plt.savefig('plot1.png')

    
    plt.figure()
    ax = plt.gca()
    fig=plt.plot(event_dic['event_counter'],relative_time,'.')
    plt.ylabel('T-T$_{0}$ (s)')
    plt.xlabel('Event counter')
    plt.title(' TimeStamp vs Event Counter')
    regress_par=linregress(event_dic['event_counter'],relative_time)
    fit_fn = np.poly1d(regress_par[0:2])
    plt.plot(event_dic['event_counter'],fit_fn(event_dic['event_counter']),'r')
    text(0.7, 0.1, 'T$_{0}$='+str(event_dic['timestamp'][0]), horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
    text(0.25, 0.65, 'Linear Regression Parameters:'+'\n'+'r$^{2}$='+str("{0:.3f}".format(regress_par[2]))+'\n'
         +'Slope:'+str("{0:.4f}".format(regress_par[0])+'s')+'\n'+'Intercept:'+str("{0:.2f}".format(regress_par[1]))+'s', horizontalalignment='center',
         verticalalignment='center', transform=ax.transAxes,bbox=dict(facecolor='white', alpha=0.5))
    
    text(0.7, 0.3, 'Estimated rate='+str("{0:.2f}".format(1/regress_par[0]))+'Hz', horizontalalignment='center',verticalalignment='center', bbox=dict(facecolor='white', alpha=0.7),transform=ax.transAxes)
    
    message=0
    s_points=[]
    s_times=[]
    for i in range(0,len(event_dic['event_counter'])-1):
        if (event_dic['event_counter'][i]>event_dic['event_counter'][i+1]) or (event_dic['event_counter'][i+1]-event_dic['event_counter'][i])>1:
            #print(event_dic['bch_counter'][i+1]-event_dic['bch_counter'][i])
            s_times.append(relative_time[i+1])
            s_points.append(event_dic['event_counter'][i+1])
            message='Suspicious points detected'
    
    plt.plot(s_points,s_times,'k.')
    
    if message==0:
        message='No suspicious points detected'
        plt.legend(['Event-Time points', 'Linear Fit'])
    else:
        plt.legend(['Event-Time points', 'Linear Fit','Suspicious points'])
        
    text(0.75, 0.4, str(message), horizontalalignment='center',
         verticalalignment='center', transform=ax.transAxes,bbox=dict(facecolor='white', alpha=0.5))
    
    plt.savefig('plot2.png')

    
    plt.figure()
    ax = plt.gca()
    fig=plt.plot(event_dic['event_busy_counter'],relative_time,'.')
    plt.title('TimeStamp vs Busy Event Counter')
    plt.ylabel('T-T$_{0}$ (s)')
    plt.xlabel('Busy Event counter')
    regress_par=linregress(event_dic['event_busy_counter'],relative_time)
    fit_fn = np.poly1d(regress_par[0:2])
    plt.plot(event_dic['event_busy_counter'],fit_fn(event_dic['event_busy_counter']),'r')
    text(0.75, 0.1, 'T$_{0}$='+str(event_dic['timestamp'][0]), horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
    text(0.25, 0.7, 'Linear Regression Parameters:'+'\n'+'r$^{2}$='+str("{0:.3f}".format(regress_par[2]))+'\n'
         +'Slope:'+str("{0:.2f}".format(regress_par[0]))+'s'+'\n'+'Intercept:'+str("{0:.2f}".format(regress_par[1]))+'s', horizontalalignment='center',
         verticalalignment='center', transform=ax.transAxes,bbox=dict(facecolor='white', alpha=0.5))
   
    
    
    message=0
    s_points=[]
    s_times=[]
    for i in range(0,len(event_dic['event_busy_counter'])-1):
        if (event_dic['event_busy_counter'][i]>event_dic['event_busy_counter'][i+1]) or (event_dic['event_busy_counter'][i+1]-event_dic['event_busy_counter'][i])>1:
            #print(event_dic['bch_counter'][i+1]-event_dic['bch_counter'][i])
            s_times.append(relative_time[i+1])
            s_points.append(event_dic['event_busy_counter'][i+1])
            message='Suspicious points detected'
            
    plt.plot(s_points,s_times,'k.')
    
    if message==0:
        message='No suspicious points detected'
        plt.legend(['BusyEvent-Time points', 'Linear Fit'])
    else:
        plt.legend(['BusyEvent-Time points', 'Linear Fit','Suspicious points'])
        
    text(0.75, 0.3, str(message), horizontalalignment='center',
         verticalalignment='center', transform=ax.transAxes,bbox=dict(facecolor='white', alpha=0.5))
    
    plt.savefig('plot3.png')

    
    plt.figure()
    ax = plt.gca()
    plt.plot(event_dic['pps'],relative_time,'.')
    plt.ylabel('T-T$_{0}$ (s)')
    plt.xlabel('PPS counter')
    plt.title('TimeStamp vs PPS Counter')
    regress_par=linregress(event_dic['pps'],relative_time)
    fit_fn = np.poly1d(regress_par[0:2])
    plt.plot(event_dic['pps'],fit_fn(event_dic['pps']),'r')
    text(0.7, 0.1, 'T$_{0}$='+str(event_dic['timestamp'][0]), horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
    text(0.25, 0.7, 'Linear Regression Parameters:'+'\n'+'r$^{2}$='+str("{0:.3f}".format(regress_par[2]))+'\n'
         +'Slope:'+str("{0:.2f}".format(regress_par[0]))+'s'+'\n'+'Intercept:'+str("{0:.2f}".format(regress_par[1]))+'s', horizontalalignment='center',
         verticalalignment='center', transform=ax.transAxes,bbox=dict(facecolor='white', alpha=0.5))

    
    message=0
    s_points=[]
    s_times=[]
    for i in range(0,len(event_dic['pps'])-1):
        if (event_dic['pps'][i]>event_dic['pps'][i+1]):
            s_times.append(relative_time[i+1])
            s_points.append(event_dic['pps'][i+1])
            message='Suspicious points detected'
    
    plt.plot(s_points,s_times,'k.')
    
    if message==0:
        message='No suspicious points detected'
        plt.legend(['PPS-Time points', 'Linear Fit'])
    else:
        plt.legend(['PPS-Time points', 'Linear Fit','Suspicious points'])
        
    text(0.7, 0.3, str(message), horizontalalignment='center',
         verticalalignment='center', transform=ax.transAxes,bbox=dict(facecolor='white', alpha=0.5))
        
    plt.savefig('plot4.png')
    
    plt.figure()
    ax = plt.gca()
    plt.plot(event_dic['event_counter'],event_dic['pps'],'.')
    plt.xlabel('Event counter')
    plt.ylabel('PPS counter')
    plt.title('PPS Counter vs Event Counter')
    regress_par=linregress(event_dic['event_counter'],event_dic['pps'])
    fit_fn = np.poly1d(regress_par[0:2])
    plt.plot(event_dic['event_counter'],fit_fn(event_dic['event_counter']),'r')
    text(0.25, 0.7, 'Linear Regression Parameters:'+'\n'+'r$^{2}$='+str("{0:.3f}".format(regress_par[2]))+'\n'
         +'Slope:'+str("{0:.2f}".format(regress_par[0]))+'s'+'\n'+'Intercept:'+str("{0:.2f}".format(regress_par[1]))+'s', horizontalalignment='center',
         verticalalignment='center', transform=ax.transAxes,bbox=dict(facecolor='white', alpha=0.5))

    
    message=0
    s_points=[]
    s_times=[]
    for i in range(0,len(event_dic['pps'])-1):
        if (event_dic['pps'][i]>event_dic['pps'][i+1]):
            s_times.append(event_dic['pps'][i+1])
            s_points.append(event_dic['event_counter'][i+1])
            message='Suspicious points detected'
    
    plt.plot(s_points,s_times,'k.')
    
    if message==0:
        message='No suspicious points detected'
        plt.legend(['PPS-Events points', 'Linear Fit'])
    else:
        plt.legend(['PPS-Events points', 'Linear Fit','Suspicious points'])
        
    text(0.7, 0.3, str(message), horizontalalignment='center',
         verticalalignment='center', transform=ax.transAxes,bbox=dict(facecolor='white', alpha=0.5))
        
    plt.savefig('plot5.png')

    
def Poisson_test(event_dic,bars):
    
    dift=event_dic['timestamp'][1:]-event_dic['timestamp'][:-1]      
    plt.figure()
    ax = plt.gca()
    hist,bin_edges=np.histogram(dift*10**(-9),bins=20)
    bin_centres=(bin_edges[1:]+bin_edges[:-1])/2 
    hist=np.array(hist)
    i=0
    while i<len(hist):
        if hist[i]==0:
            text(0.65, 0.75,'Zero values in one or more bins detected'+'\n'+'(Removed for the fitness)', horizontalalignment='center',
                 verticalalignment='center', transform=ax.transAxes,bbox=dict(facecolor='white', alpha=1))
            hist=np.delete(hist,i)
            bin_centres=np.delete(bin_centres,i)
        else:
            i=i+1
        
    plt.plot(bin_centres,np.log(hist),'.',color='blue')
    plt.errorbar(bin_centres,np.log(hist),yerr=np.sqrt(hist)/hist,xerr=None,fmt='none',color='blue')
    regress_par=linregress(bin_centres, np.log(hist))
    fit_fn = np.poly1d(regress_par[0:2])
    plt.plot(bin_centres,fit_fn(bin_centres),'r')
    plt.ylabel('log(n[Number of events])')
    plt.xlabel('$\Delta$t (s)')
    plt.title('Time difference distribution between consecutive events')
    plt.legend(['Numer of events', 'Linear Fit'])
    text(0.3, 0.2, 'Linear Regression Parameters:'+'\n'+'r$^{2}$='+str("{0:.3f}".format(regress_par[2]))+'\n'
             +'Slope:'+str("{0:.2f}".format(regress_par[0]))+'s$^{-1}$'+'\n'+'Intercept:'+str("{0:.2f}".format(regress_par[1])), horizontalalignment='center',
             verticalalignment='center', transform=ax.transAxes,bbox=dict(facecolor='white', alpha=0.5))
    text(0.7, 0.6, 'Estimated rate:'+str("{0:.2f}".format(-regress_par[0])+'Hz'), horizontalalignment='center',
             verticalalignment='center', transform=ax.transAxes,bbox=dict(facecolor='white', alpha=0.5))
    plt.savefig('plot6.png')
    
    #If we want to see the histogram with bars
    if bars==True:
        plt.figure()
        plt.hist(dift,bins=20)
        plt.xlabel('$\Delta$t (s)')
        plt.ylabel('n')
        mean_data=np.median(dift)
        sigma_data=np.std(dift)
        text(0.7, 0.9, 'Median:'+str("{0:.2f}".format(mean_data))+'\n'+'$\sigma$='+str("{0:.9f}".format(sigma_data)), horizontalalignment='center',
             verticalalignment='center', transform=ax.transAxes,bbox=dict(facecolor='white', alpha=0.5))
        plt.title('Time difference distribution between consecutive events')
        plt.savefig('plot7.png')
    
def dif_time_test(event_dic):
    difftime=event_dic['timestamp']-event_dic['system_time']
    ax = plt.gca()
    plt.figure()
    plt.hist(difftime,bins=20)
    plt.xlabel('$\Delta$TimeStamp (UCTS-System) (s)')
    plt.ylabel('n')
    mean_data=np.median(difftime*10**(-9))
    sigma_data=np.std(difftime*10**(-9))
    text(0.7, 0.9, 'Median:'+str("{0:.2f}".format(mean_data))+'s'+'\n'+'$\sigma$='+str("{0:.4f}".format(sigma_data))+'s', horizontalalignment='center',
         verticalalignment='center', transform=ax.transAxes,bbox=dict(facecolor='white', alpha=0.5))
    plt.title('Time difference between UCTS Time Stamp and the System Unix Time Stamp')
    plt.savefig('plot8.png')

def execute_tests(event_dic,bars):
    pdf = FPDF()
    Infovstime(event_dic)
    Poisson_test(event_dic,bars)
    dif_time_test(event_dic)
    if bars==True:
        imagelist=['plot1.png','plot2.png','plot3.png','plot4.png','plot5.png','plot6.png','plot7.png','plot8.png']
    else:
        imagelist=['plot1.png','plot2.png','plot3.png','plot4.png','plot5.png','plot6.png','plot8.png']
    pdf.add_page()
    for image in imagelist:
        pdf.image(image,x=20,y=None,w=170)
    pdf.output('Timing_plots_'+time.strftime("%d")+'_'+time.strftime("%m")+'_'+time.strftime("%Y")+'_'+time.strftime("%H:%M:%S")+'.pdf', "F")


if __name__ == '__main__':
    bars=True
    event_dic=read_file('Eventos04_11_19_11_27_14.txt')
    execute_tests(event_dic,bars)
    

    
   
    
    
    
    
    
    
 

    
    
    
    
