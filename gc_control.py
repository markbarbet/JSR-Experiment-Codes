# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 15:53:20 2020

@author: Mark Barbet
"""


import requests
import json
import os
import pandas as pd
import sys
import numpy as np



class chromatograph():
    
    
    def __init__(self,ip:str='10.10.0.1'):
        self.ip='http://'+ip
        self.gc=requests.session()
        r=self.gc.get(self.ip)
        if str(r)=='<Response [200]>':
            print('Connection with GC successful.')
        else:
            print('Failed to connect to GC.')
        
    def set_calibration_species(self,cal_assignments:dict):
        self.calibration_species=cal_assignments
        self.species_data={'A':{},'B':{},'C':{},'D':{}}
        self.species_area_data={'A':{},'B':{},'C':{},'D':{}}
        for i,col in enumerate(self.calibration_species.keys()):
            #print(self.calibration_species[col])
            for k, spec in enumerate(self.calibration_species[col]):
                self.species_data[col][spec]=[]
                #self.species_area_data={'A':{},'B':{},'C':{},'D':{}}
                self.species_area_data[col][spec]=[]
        
    def bakeout(self,time:int=1200):
        self.get(self.ip+'/v1/scm/sessions/system-manager!cmd.bakeout?duration='+str(time)+'s')
    def get_run_id(self,name='DefaultRunName'):
        id_list=self.gc.get(self.ip+'/v1/runData?text='+name).json()['runs']
        if len(id_list)==1:
            run_id=id_list[0]['$id']
        elif len(id_list)>1:
            index=int(input('Multiple runs with same name for '+name+' - please specify count of file to use, starting with index 0:'))
            run_id=id_list[index]['$id']
        elif len(id_list)<1:
            sys.exit('No run with name '+name+', exiting process.  Check filenames.')
        return run_id
        
    def get_data(self,name='DefaultRunName'):
        
        if name != '':
            data=self.gc.get(self.ip+'/v1/runData/'+self.get_run_id(name)).json()
            return data
        
        elif name=='':
            print('Please provide a file name')
            
    def get_calibrated_species_count(self,calibration_species:dict,data:dict):
        
        for i,column in enumerate(calibration_species.keys()):
            for k,species in enumerate(calibration_species[column]):
                tempdata=data['detectors']['module'+column+':tcd']['analysis']['peaks']
                for a, dic in enumerate(tempdata):
                    #for j,val in enumerate(tempdata[a].keys()):
                        
                        if 'label' in tempdata[a].keys():
                            if tempdata[a]['label']==species:
                                #print(tempdata[a])
                                self.species_data[column][species].append(tempdata[a]['concentration'])
                                #print(a,j,species)
    def get_calibrated_species_areas(self,calibration_species:dict,data:dict):
        
        for i,column in enumerate(calibration_species.keys()):
            for k,species in enumerate(calibration_species[column]):
                tempdata=data['detectors']['module'+column+':tcd']['analysis']['peaks']
                for a, dic in enumerate(tempdata):
                    #for j,val in enumerate(tempdata[a].keys()):
                        
                        if 'label' in tempdata[a].keys():
                            if tempdata[a]['label']==species:
                                #print(tempdata[a])
                                self.species_area_data[column][species].append(tempdata[a]['area'])
                                #print(a,j,species)
    
    
    def get_calibration_data(self,method_name:str,calibration_species:dict):
        #print(self.gc.get((self.ip+'/v1/methods/userMethods/'+method_name)).json())
        calibration=self.gc.get((self.ip+'/v1/methods/userMethods/'+method_name)).json()['peakParameters']['calibration']['detectors']
        #inverted_species_columns={}
        species=[]
        #print(calibration_species)
        for i,col in enumerate(list(calibration_species.keys())):
            for j,spec in enumerate(calibration_species[col]):
                if spec not in species:
                    species.append(spec)
        caldata={}
        #print(calibration_species)
        for i,spec in enumerate(species):
            caldata[spec]={}
            for j,col in enumerate(list(calibration_species.keys())):
                if spec in list(calibration_species[col]):
                    caldata[spec][col]={'areas':[],'known_x':[]}
        
        for i,col in enumerate(list(calibration.keys())):
            module=col.split(':')[0][-1]
            for j,compound in enumerate(calibration[col]['calibrationPeaks']):
                cal_species=compound['compoundName']
                for k,calpoint in enumerate(compound['calibrationPoints']):
                    if 'knownConcentration' in list(calpoint.keys()) and 'area' in list(calpoint.keys()) and cal_species in species:
                        if module in list(calibration_species.keys()):
                            #print('poooop')
                            if cal_species in calibration_species[module]:
                                #print('poooop')
                                caldata[cal_species][module]['areas'].append(calpoint['area'])
                                caldata[cal_species][module]['known_x'].append(calpoint['knownConcentration'])
                if module in list(calibration_species.keys()):
                    if cal_species in calibration_species[module]:
                        caldata[cal_species][module]['fitType']=compound['fitType']
        return caldata
    
                            
    def load_method(self,name:str):
        self.gc.get(self.ip+'/v1/scm/sessions/system-manager!cmd.loadMethod?methodLocation=/v1/methods/userMethods/'+name)                 
    def run_method(self,optional_parameters:dict={}):
        if optional_parameters=={}:
            response=self.gc.get(self.ip+'/v1/scm/system-manager!cmd.run?runWhenReady=true')
        else:
            #jsonData = json.dumps(optional_parameters)
            response = self.gc.post(self.ip+'/v1/scm/system-manager!cmd.run',json=optional_parameters)
    
    def build_run_parameters(self,runWhenReady:bool=True,name:str='DefaultRunName',tags=[]):
        option_dict={}
        option_dict['runWhenReady']=runWhenReady
        option_dict['annotations']={}
        option_dict['annotations']['name']=name
        option_dict['annotations']['tags']=tags
        return option_dict
    
        
    def abort_run(self):
        self.gc.get(self.ip+'/v1/scm/sessions/system-manager!cmd.abort')
        
    def get_species_averages(self):
        self.species_averages={'A':{},'B':{},'C':{},'D':{}}
        for i, col in enumerate(self.species_data.keys()):
            for j, species in enumerate(self.species_data[col].keys()):
                self.species_averages[col][species]=np.mean(np.array(self.species_data[col][species]))
                
    def get_species_areas(self):
        self.species_areas={'A':{},'B':{},'C':{},'D':{}}
        for i, col in enumerate(self.species_area_data.keys()):
            for j, species in enumerate(self.species_area_data[col].keys()):
                self.species_areas[col][species]=np.mean(np.array(self.species_area_data[col][species]))
                
            
        
    
        

    
                            
        