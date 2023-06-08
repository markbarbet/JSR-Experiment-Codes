# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 13:12:00 2022

@author: mcb22
"""

import flow_mixture_uncertainties as fmu
import numpy as np
from sklearn.linear_model import LinearRegression
import copy
from scipy import optimize
import matplotlib.pyplot as plt

class calibration_curve():
    
    def __init__(self,mixture):
        
        
        self.mixture=mixture
        self.caldata=mixture.calibration_data
        #self.low_caldata=self.get_low_caldata(self.caldata)
        #self.high_caldata=self.get_high_caldata(self.caldata)
        #self.slope=slope
        #self.intercept=intercept
        for i,species in enumerate(list(self.caldata.keys())):
            for k,col in enumerate(list(self.caldata[species].keys())):
                self.caldata=self.get_fit_parameters(col, species, self.caldata)
                #self.low_caldata=self.get_fit_parameters(col, species, self.low_caldata)
                #self.high_caldata=self.get_fit_parameters(col, species, self.high_caldata)
                
    def add_calibration_mixtures(self,calmixes:list):
        self.cal_mixtures=calmixes
        
    def set_new_caldata(self,caldata):
        self.caldata=caldata
        
    def optimization_linear(self,areas,slope,intercept):
        
        return np.multiply(slope,np.array(areas))+intercept
    
    def get_calmix_std_devs(self,nominal,species):
        tempval=0.0
        for i,mix in enumerate(self.cal_mixtures):
            #print(mix.calibration_nominals)
            #print(nominal)
            if species in mix.calibration_nominals.keys():
                if mix.calibration_nominals[species]==nominal:
                    tempval= mix.mixture_uncertainties[species]*100.0
        if tempval==0.0:
            tempval=2.0
            
        return tempval
            
            
    def monte_carlo_fitting(self,caldata,nTrials=4000):
        #y_stds=copy.deepcopy(caldata)
        for i,spec in enumerate(list(caldata.keys())):
            for k,col in enumerate(list(caldata[spec].keys())):
                x_sets=[]
                caldata[spec][col]['y_std']=[]
                for j,value in enumerate(caldata[spec][col]['known_x']):
                    caldata[spec][col]['y_std'].append(self.get_calmix_std_devs(caldata[spec][col]['known_x'][j], spec))
                    current_val=value
                    if j==0:
                        x_sets.append(1)
                    elif j>0:
                        if current_val==caldata[spec][col]['known_x'][j-1]:
                            x_sets[-1]=x_sets[-1]+1
                        elif current_val!=caldata[spec][col]['known_x'][j-1]:
                            x_sets.append(1)
                aFitPars=np.array([])
                for iTrial in range(nTrials):
                    xTrial=[]
                    yTrial=[]
                    x_set_counter=0
                    #x_set_value=0
                    j_pause_val=0
                    for j,value in enumerate(caldata[spec][col]['areas']):
                        xTrial.append(np.random.normal(np.mean(caldata[spec][col]['areas'][j_pause_val:np.sum(x_sets[0:x_set_counter+1])]),
                                                      np.std(caldata[spec][col]['areas'][j_pause_val:np.sum(x_sets[0:x_set_counter+1])])))
                        if j==np.sum(x_sets[0:x_set_counter+1])-1:
                            x_set_counter=x_set_counter+1
                            j_pause_val=j+1
                        #if spec=='N2O':
                            #print(caldata[spec][col]['slope']*xTrial[-1])
                            #print(caldata[spec][col]['intercept'])
                            #print(np.random.normal(caldata[spec][col]['known_x'][j],caldata[spec][col]['y_std'][j]))
                        #yTrial.append(caldata[spec][col]['slope']*xTrial[-1]+caldata[spec][col]['intercept']+np.random.normal(0,
                                                                                                                              #caldata[spec][col]['y_std'][j]))
                        yTrial.append(np.random.normal(caldata[spec][col]['known_x'][j],caldata[spec][col]['y_std'][j]))
                    xTrial=np.array(xTrial)
                    yTrial=np.array(yTrial)
                    
                    
                    try:
                        vTrial,aCova = optimize.curve_fit(self.optimization_linear,xTrial,yTrial,[caldata[spec][col]['slope'],caldata[spec][col]['intercept']])
                    except:
                        dumdum=1
                        continue
                    
                    if np.size(aFitPars)<1:
                        aFitPars=np.copy(vTrial)
                    else:
                        aFitPars=np.vstack((aFitPars,vTrial))
                        
                caldata[spec][col]['Monte-Carlo-Parameters']=aFitPars
                caldata[spec][col]['slope-standard-dev']=np.std(aFitPars[:,0])
                caldata[spec][col]['intercept-standard-dev']=np.std(aFitPars[:,1])
        #return aFitPars
            
                        
                        
                
    def plot_montecarlo_distributions(self,header):
        
        for i,spec in enumerate(list(self.caldata.keys())):
            for k,col in enumerate(list(self.caldata[spec].keys())):
                plt.figure()
                plt.hist(self.caldata[spec][col]['Monte-Carlo-Parameters'][:,0],bins=50,color='b',edgecolor='k')
                plt.xlabel('Area coefficient a')
                plt.ylabel('N(a)')
                plt.title(spec+':'+col+' Calibration Curve Coefficient distribution')
                plt.savefig(header+'_'+spec+'_'+col+'slope_hist.pdf',dpi=600,bbox_inches='tight')
                
                
                plt.figure()
                plt.hist(self.caldata[spec][col]['Monte-Carlo-Parameters'][:,1],bins=50,color='b',edgecolor='k')
                plt.xlabel('Intercept b')
                plt.ylabel('N(b)')
                plt.title(spec+':'+col+' Calibration Curve Intercept distribution')
                plt.savefig(header+'_'+spec+'_'+col+'intercept_hist.pdf',dpi=600,bbox_inches='tight')
                
                plt.figure()
                plt.scatter(self.caldata[spec][col]['Monte-Carlo-Parameters'][:,0],self.caldata[spec][col]['Monte-Carlo-Parameters'][:,1],alpha=0.5,s=9,edgecolor='none')
                plt.xlabel('Area coefficient a')
                plt.ylabel('Intercept b')
                plt.title(spec+':'+col+' Calibration Curve Coefficient vs. Intercept')
                plt.savefig(header+'_'+spec+'_'+col+'slope_intercept_scatter.pdf',dpi=600,bbox_inches='tight')
        
     
    def linear_fit(self,x,y,caldata,species,module,fit_intercept=True):
        model=LinearRegression(fit_intercept=fit_intercept)
        model.fit(x,y)
        r_sq=model.score(x,y)
        slope=model.coef_[0]
        intercept=model.intercept_
        caldata[species][module]['slope']=slope
        caldata[species][module]['intercept']=intercept
        caldata[species][module]['R-squared']=r_sq
        return caldata
        
    def get_fit_parameters(self,module,species,caldata):
        xdata=np.array(caldata[species][module]['areas']).reshape((-1,1))
        ydata=np.array(caldata[species][module]['known_x'])
        if caldata[species][module]['fitType']=='linear':
            func=getattr(self,caldata[species][module]['fitType']+'_fit')
            caldata=func(xdata,ydata,caldata,species,module)
            
        return caldata
    
    def get_high_caldata(self,caldata):
        percent_error={}
        multipliers={}
        high_caldata=copy.deepcopy(caldata)
        for i, spec in enumerate(list(caldata.keys())):
            percent_error[spec]=np.divide(self.mixture.mixture_uncertainties[spec],self.mixture.component_fractions[spec])
            multipliers[spec]=1+percent_error[spec]
        
        for i,spec in enumerate(list(caldata.keys())):
            for j,col in enumerate(list(caldata[spec].keys())):
                high_caldata[spec][col]['areas']=list(np.array(high_caldata[spec][col]['areas'])*multipliers[spec])
                
        return high_caldata
        
    def get_low_caldata(self,caldata):
        percent_error={}
        multipliers={}
        low_caldata=copy.deepcopy(caldata)
        for i, spec in enumerate(list(caldata.keys())):
            percent_error[spec]=np.divide(self.mixture.mixture_uncertainties[spec],self.mixture.component_fractions[spec])
            multipliers[spec]=1-percent_error[spec]
        
        for i,spec in enumerate(list(caldata.keys())):
            for j,col in enumerate(list(caldata[spec].keys())):
                low_caldata[spec][col]['areas']=list(np.array(low_caldata[spec][col]['areas'])*multipliers[spec])
                
        return low_caldata
            
        
        
        
        