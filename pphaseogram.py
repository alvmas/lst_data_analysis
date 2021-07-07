#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  4 11:27:05 2021

@author: alvarom
"""

import pandas as pd
import math
import astropy as ast
import argparse
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import matplotlib as mpl
import os
from astropy.time import Time
from copy import copy
from decimal import *
from pylab import *
from scipy.stats import linregress
from astropy import units as u
import h5py
from gammapy.stats import WStatCountsStatistic
import warnings
from lstchain.reco.utils import get_effective_time,add_delta_t_key
from scipy.stats import chi2,chisquare, norm



def calculate_CountStats(on_file,off_file=None,factor=None):
    if off_file is None:
        raise ValueError('No off data given for the pulsar analysis')
        
    Non=len(on_file)
    Noff=len(off_file)
    stat = WStatCountsStatistic(n_on=Non, n_off=Noff, alpha=factor)
    yerr=np.sqrt(Non + ((factor** 2) * Noff))
    
    return(stat,yerr,Noff*factor)


class PulsarPeak():
    '''
    A class to manipulate and store information about the statistics of certain regions in a phase list object. 2 types of regions:
        1. Signal region: peaks present in the lightcurve
        2. Background region: background statistics
    
    
    Parameters
    ----------
    lc_data : PulsarPhases class object
        PulsarPhases object from which we extract the pulsar phases and the OFF region statistics
    peak_limits : list
        list of phase edges used to define the region.
        To define a non continuos region (for instance, the sum of two independent peaks), a list of four edges can be provided.
    peak_type : str, optional
        'Background' or 'Signal' type to calculate especific statistics.
    
    Note: To obtain peak statistics an OFF region must be defined in the phase list object first.
    
    
    Attributes
    ----------
    limits: str
        region phase edges
    type: str
        the type of region defined (signal or background)
    deltaP: float
        the total phase range of the region
    phases: numpy array
        list of phases that fall into the region
    
    For 'signal' type only:

    sign: float
        Li&Ma Significance of the excess of events from the peak with respect to the background
    Nex: int
        Number of excess events in the peaks
    yerr: float
        Estimated error of the number of excess events
    sign_ratio: float
        Ratio of the significance and the square root of the time of observation.
    s_n_ratio: float
        Signal to noise ratio
    '''


    def __init__(self, lc_data, peak_limits,peaktype='signal'):
        self.limits=peak_limits
        self.make_stats(lc_data)
        self.type=peaktype
        
    def make_stats(self,lc_data):
        if len(self.limits)>2:
            self.deltaP=(self.limits[1]-self.limits[0])+(self.limits[3]-self.limits[2])
            self.phases=np.concatenate((lc_data.phases[(lc_data.phases>self.limits[0]) & (lc_data.phases<self.limits[1])],lc_data.phases[(lc_data.phases>self.limits[2]) & (lc_data.phases<self.limits[3])]))
        else:
            self.deltaP=self.limits[1]-self.limits[0]
            self.phases=lc_data.phases[(lc_data.phases>self.limits[0]) & (lc_data.phases<self.limits[1])]
        self.number=len(self.phases)
        
        if self.type=='signal':
            stats,yerror,noff=calculate_CountStats(self.phases,off_file=lc_data.OFF.phases,factor=(self.deltaP)/lc_data.OFF.deltaP)
            self.sign=stats.sqrt_ts.item()
            self.Nex=stats.n_sig
            self.yerr=yerror
            self.sign_ratio=self.sign/np.sqrt(lc_data.tobs)
            self.s_n_ratio=self.Nex/noff
            
    def change_limits(new_limits,lc_data):
        self.limits=new_limits
        self.make_stats(lc_data)
        
    
        
class PhaseBinning():
    '''
    A class to define the binning used in the construction of the lightcurve. It allows to change, manipulate and find the phase binning of the histogram
    
    
    Parameters
    ----------
    bins : int, array
        If an integer is provided, it is used as the number fix-width bins in the lightcurve.
        If an array are provided, that list will be used as the bin edges of the lightcurve.
    xmin : float
        Lower edge of the lightcuve binning
    xmax : float
        Higher edge of the lightcuve binning
    
    
    Attributes
    ----------
    nbins : int
        Number of bins used in the lightcurve
    xmin : float
        Lower edge of the lightcuve binning
    xmax : float
        Higher edge of the lightcuve binning
    bins: list
        List of bin edges of the lightcuve
    '''
    
    def __init__(self, bins,xmin=None,xmax=None):
        
        if isinstance(bins, int) :
            self.nbins=bins
            if xmin is not None:
                self.xmin=xmin
            else:
                self.xmin=-0.5  
            if xmax is not None:
                self.xmax=xmax
            else:
                self.xmax=0.5
            self.set_edges()
        
        elif isinstance(bins, (list, np.ndarray)):
            self.edges=bins
            self.nbins=len(bins)-1
            self.xmin=bins[0]
            self.xmax=bins[-1]
    
    
    def getNumEdges(self):
        return(self.nbins+1)
    
    def set_edges(self,nbins=None,xmin=None,xmax=None):
            
        if nbins is not None:
            self.nbins=nbins
            if xmax is not None:
                self.xmax=xmax
            if xmin is not None:
                self.xmin=xmin
        else:
            if xmax is not None:
                self.xmax=xmax
            if xmin is not None:
                self.xmin=xmin
        self.edges=np.linspace(self.xmin,self.xmax,self.nbins+1)
        
                
    def Find_LowHiEdge(self,value):
        for i in range(1,self.getNumEdges()):
            if self.edges[i]>=value:
                return(i-1,i)
        return(None)
        
        
    def Find_CloseEdge(self,value):
        if value<self.xmin:
            return(0)
        if value>self.xmax:
            return(self.getNumEdges())
        for i in range(1,self.getNumEdges()):
            if self.edges[i]>=value:
                if (self.edges[i]-value)<(value-self.edges[i-1]):
                    return(i)
                else:
                    return(i-1)

class PeriodicityTest():
    '''
    A class to apply and store the information of the periodicty tests.
    
    Parameters
    ----------
    lc_data : Lightcurve object
        Object from which we extract the lightcurve and apply the binned tests (chi square)
    p_data : PulsarPhases object
        Object from which we extract the pulsar phases and apply the unbinned tests (Zn, H tests)
    
    
    Attributes
    ----------
    chisqr_test :
        Results of the chi square tests. Format: [Statistic, p_value, nsigmas]
    number : int
        Number of phases used in the analysis
    cos : list of float
        Cosine moments of the list of phases
    sin: list of float
        Sine moments of the list of phases
    Zntest_res:
        Results of the Zn tests. Format: [Statistic, p_value, nsigmas]. Default is n=10
    Htest_res:
        Results of the H test. Format: [Statistic, p_value, nsigmas]
    '''
    
    
    
    def __init__(self,lc_data,p_data):
            self.apply_all_tests(lc_data,p_data)
        
    def apply_all_tests(self,lc_data,p_data):
            self.chisqr_res=lc_data.chi_sqr_pulsar_test()
            self.number,self.cos,self.sin=p_data.moments()
            self.apply_moment_tests(p_data)
            self.resume_stats()
    
    def calculate_moments(self,p_data,n):
            self.number,self.cos,self.sin=p_data.moments(n)
    
    def apply_moment_tests(self,p_data):
            self.Zntest_res=p_data.zn_test(self.cos,self.sin,self.number)
            self.Htest_res=p_data.H_test(self.cos,self.sin,self.number)
            
    def resume_stats(self):
            self.resume=pd.DataFrame(data={'Chi_square_test':self.chisqr_res,'Zn_test':self.Zntest_res,'H_test':self.Htest_res},index=["Statistic", "p-value", "Number of $\sigma$"])



class PulsarPhases():

    '''
    MAIN CLASS FOR THE PULSAR ANALYSIS.
    A class to store the pulsar phases and mjd_times to be used in the LST analysis. This class allows to develop all the timing pular analysis using different supporting classes and subclasses.
    
    Parameters
    ----------
    dataframe : DL2 pandas dataframe
        DL2 LST file after the quality selection. Set daraframe to False if it is not available and want to set the attributes manually.
    pdata : List of float
        List of phases (in case no dataframe is available).
    ptimes: List of float
        List of mjd times (in case no dataframe is available).
    tobservation: float
        Total effective time of observation
    peak_limits_1: tuple
        Edges of the first peak. Set to None if no P1 is present.
    peak_limits_2: tuple
        Edges of the second peak. Set to None if no P2 is present.
    off_limits: tuple
        Edges of the OFF region
    add_stats: boolean
        True to do the quantitative statistics
        
    Attributes
    ----------
    phases : list of float
        List of pulsar phases.
    times : list of float
        List of mjd times
    tobs : float
        Effective time of observation in hours
    OFF: PulsarPeak object
        Information of the OFF region
    P1/P2/P1P2: PulsarPeak object
        Information of the signal regions.
    total_sig: float
        Total significance of the signal region
    '''
    

    def __init__(self, dataframe, pdata=None, ptimes=None, tobservation=None,peak_limits_1=[0.983-1,0.026],peak_limits_2=[0.377,0.422],off_limits=[-0.48,0.87-1],add_stats=True):
    
        if dataframe is None:
                self.phases=np.array(pdata)
                self.times=np.array(ptimes)
                self.tobs=tobservation
        else:
                self.phases=np.array(dataframe['pulsar_phase'].to_list())
                self.times=np.array(dataframe['mjd_time'].to_list())
                self.tobs=get_effective_time(dataframe)[1].value/3600
                
        self.OFF=PulsarPeak(self,off_limits,peaktype='background')
        self.P1=None
        self.P2=None
            
        if peak_limits_1 is not None:
            self.P1=PulsarPeak(self,peak_limits_1)
        else:
            print('No P1 limits. Cant create P1 object')
        if peak_limits_2 is not None:
            self.P2=PulsarPeak(self,peak_limits_2)
        else:
            print('No P2 limits. Cant create P2 object')
        
        if add_stats==True:
            if self.P1 is None:
                self.total_sig=self.P2.sign
            else:
                if self.P2 is None:
                    self.total_sig=self.P1.sign
                else:
                    self.P1P2=PulsarPeak(self,peak_limits_1+peak_limits_2)
                    
        self.histogram=Lightcurve(self)
    
    
    def moments(self,n=25):
        plist=(self.phases+0.5)*2*np.pi
        k=np.arange(1,n+1)
        cos_moment=sum(np.cos(np.outer(plist,k)),axis=0)
        sin_moment=sum(np.sin(np.outer(plist,k)),axis=0)
        
        return(len(plist),cos_moment,sin_moment)
        
        
    def zn_test(self,cos,sin,number,n=10):
        cos_moment=cos[0:n+1]
        sin_moment=sin[0:n+1]
        
        Zn=2/number*sum(np.power(cos_moment,2)+np.power(sin_moment,2))
        pvalue_zn=chi2.sf(float(Zn),2*n)
        sigmas_zn=norm.isf(pvalue_zn, loc=0, scale=1)
        
        return(Zn,pvalue_zn,sigmas_zn)

    def H_test(self,cos,sin,number):
        bn=0.398
        h=[]
        
        for m in np.arange(1,len(cos)):
            h.append(2/number*sum(np.power(cos[0:m],2)+np.power(sin[0:m],2))-4*m+4)
        
        H=max(h)
        pvalue_H=np.exp(-bn*H)
        sigmas_H=norm.isf(float(pvalue_H), loc=0, scale=1)
    
        return(H,pvalue_H,sigmas_H)
        
    def show_peak_results(self):
        peak_results={'P1+P2':self.P1P2.sign,'P1':self.P1.sign,'P2':self.P2.sign}
        return(peak_results)
    
    def sigVStime(self):
        t=[0]
        k=0
        s=0
        
        diff=abs(self.times[1:]-self.times[:-1])
        pos=diff>0.0005
        index=np.where(pos)[0]
        index=list(index)
        index.append(len(self.times)-2)

        sign_t=pd.DataFrame(data={'P1P2':[],'P1':[],'P2':[]})
        excess_t=pd.DataFrame(data={'P1P2':[],'P1':[],'P2':[]})
        error_t=pd.DataFrame(data={'P1P2':[],'P1':[],'P2':[]})
        
        for i in index:
                phasest=self.phases[(phases_t.times>self.times[k]) & (phases_t.times<self.times[i+1])]
                timesp=self.times[(phases_t.times>self.times[k]) & (phases_t.times<self.times[i+1])]
                phases_t=PulsarPhases(phasest,timesp)
                
                s=s+phases_t.tobs
                t.append(s)
                k=i
                
                phases_t.OFF=PulsarPeak(phases_t)
                P1P2=PulsarPeak(phases_t,self.P1P2.peak_limits)
                P1=PulsarPeak(phases_t,self.P1.peak_limits)
                P2=PulsarPeak(phases_t,self.P2.peak_limits)
                
                sign_t=pd.concat([sign_t,pd.DataFrame(data={'P1P2':P1P2.sign,'P1':P1.sign,'P2':P2.sign})])
                excess_t=pd.concat([excess_t,pd.DataFrame(data={'P1P2':P1P2.Nex,'P1':P1.Nex,'P2':P2.Nex})])
                error_t=pd.concat([error_t,pd.DataFrame(data={'P1P2':P1P2.yerr,'P1':P1.yerr,'P2':P2.yerr})])
        
        return(t,sign_t,excess_t,error_t)
        
    def show_periodicity_stats(self):
        try:
            return(self.stats)
        except:
            print('No statistics available. Make sure to include add_statistcis=True if you want to calculate them.')
            
            
            
        
class Lightcurve():

    '''
    A class to construct the lightcurve from the information provided by the rest of classes.
    
    Parameters
    ----------
    pulsar_phases: PulsarPhase object
        Object from which we extract the pulsar phases to construct the lightcurve
    nbins : int, array
        Number of bins used in the lightcurve if int
        Bin edges used in the lightcurve if array
    add_stats: boolean
        True to do the quantitative statistics
    
    Attributes
    ----------
    binning : PhaseBinning object
        Binning information
    stats : PeriodicityTest object
        Information of the statistical tests.
    lc: numpy histogram object
        Information of the pular phase histogram
    '''
    

    def __init__(self, pulsar_phases,nbins=50,add_stats=True):
        self.binning=PhaseBinning(nbins)
        self.create_histogram(pulsar_phases.phases)
        
        if add_stats==True:
            self.stats=PeriodicityTest(self,pulsar_phases)
        else:
            print('Not adding periodicity tests')
    
    def create_histogram(self,phases):
        self.lc=np.histogram(phases, bins=self.binning.edges)
        

    def chi_sqr_pulsar_test(self):
        Nbins=len(self.lc[0])
        mean_signal=np.mean(self.lc[0])
        chi_sqr,p_value=chisquare(self.lc[0],mean_signal)
        sigmas=norm.isf(p_value*10**(-7), loc=0, scale=1)
        
        return(chi_sqr,p_value,sigmas)
    
    
    def show_phaseogram(self,pulsar_phases,phase_limits=[0,1],show_stats=True):
    
        plt.figure(figsize=(15,5))
        
        plt.bar((self.lc[1][1:]+self.lc[1][:-1])/2,self.lc[0],width=1/self.binning.nbins,color='blue',alpha=0.5)
        plt.bar((self.lc[1][1:]+self.lc[1][:-1])/2+np.ones(len(self.lc[1][:-1])),self.lc[0],width=1/self.binning.nbins,color='blue',alpha=0.5)
        plt.errorbar((self.lc[1][1:]+self.lc[1][:-1])/2,self.lc[0],yerr=np.sqrt(self.lc[0]),color='blue',fmt='.')
        plt.errorbar((self.lc[1][1:]+self.lc[1][:-1])/2+np.ones(len(self.lc[1][:-1])),self.lc[0],yerr=np.sqrt(self.lc[0]),color='blue',fmt='.')
        
        plt.fill_between(np.linspace(pulsar_phases.OFF.limits[0],pulsar_phases.OFF.limits[1],150), 0,1600500,facecolor="black",color='black',alpha=0.2,label='OFF')
        plt.fill_between(np.linspace(pulsar_phases.OFF.limits[0]+1,pulsar_phases.OFF.limits[1]+1,150), 0,1600500,facecolor="black",color='black',alpha=0.2)
        
        
        if pulsar_phases.P1 is not None:
            plt.fill_between(np.linspace(pulsar_phases.P1.limits[0],pulsar_phases.P1.limits[1],150), 0, 1600500,facecolor="orange",color='orange',alpha=0.2,label='P1')
            plt.fill_between(np.linspace(pulsar_phases.P1.limits[0]+1,pulsar_phases.P1.limits[1]+1,150), 0, 1600500,facecolor="orange",color='orange',alpha=0.2)
        
        if pulsar_phases.P2 is not None:
            plt.fill_between(np.linspace(pulsar_phases.P2.limits[0],pulsar_phases.P2.limits[1],150), 0,1600500,facecolor="green",color='green',alpha=0.2,label='P2')
            plt.fill_between(np.linspace(pulsar_phases.P2.limits[0]+1,pulsar_phases.P2.limits[1]+1,150), 0,1600500,facecolor="green",color='green',alpha=0.2)
        
        plt.hlines(y=np.mean((self.lc[0][(self.lc[1][:-1]>(pulsar_phases.OFF.limits[0])) & (self.lc[1][1:]<(pulsar_phases.OFF.limits[1]))])),xmin=-0.5,xmax=1.5,linestyle='dashed',color='k')
        plt.text(0.3,max(self.lc[0])+3*np.sqrt(min(self.lc[0])),f'Tobs={pulsar_phases.tobs:.1f} h',fontsize=15,bbox=dict(facecolor='white',edgecolor='black'))
        plt.xlabel('Pulsar phase')
        plt.ylabel('Events')
        plt.ylim(min(self.lc[0])-3*np.sqrt(max(self.lc[0])),max(self.lc[0])+2*np.sqrt(max(self.lc[0])))
        plt.xlim(phase_limits[0],phase_limits[1])
        plt.legend(loc=2,fontsize=15)
        
        
        if show_stats==True:
            plt.text((phase_limits[1]+phase_limits[0])/2,min(self.lc[0])-7*np.sqrt(min(self.lc[0])),f'$\chi^{2}$-test: $\chi^{2}$={pulsar_phases.stats.resume.Chi_square_test[0]:.2f} p_value={"{:.2e}".format(pulsar_phases.stats.resume.Chi_square_test[1])} sign={pulsar_phases.stats.resume.Chi_square_test[2]:.2f}$\sigma$ \n H-test: H={pulsar_phases.stats.resume.H_test[0]:.2f} p_value={"{:.2e}".format(pulsar_phases.stats.resume.H_test[1])} sign={pulsar_phases.stats.resume.H_test[2]:.2f}$\sigma$ \n Z$_{{10}}$-test: Z$_{{10}}$={pulsar_phases.stats.resume.Zn_test[0]:.2f} p_value={"{:.2e}".format(pulsar_phases.stats.resume.Zn_test[1])} sign={pulsar_phases.stats.resume.Zn_test[2]:.2f}$\sigma$ ',color='purple',fontsize=15,bbox=dict(facecolor='white',edgecolor='black'))
            plt.text(phase_limits[0],min(self.lc[0])-7*np.sqrt(min(self.lc[0])),f' P1+P2: Sig(Li&Ma):{pulsar_phases.P1P2.sign:.2f}$\sigma$ \n P1: Sig(Li&Ma):{pulsar_phases.P1.sign:.2f}$\sigma$ \n P2: Sig(Li&Ma):{pulsar_phases.P2.sign:.2f}$\sigma$',color='red',fontsize=15,bbox=dict(facecolor='white',edgecolor='red'))
            
        plt.show()
        
        

        
        
        
        
