# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 16:43:34 2022

@author: Mark Barbets
"""

import os
import numpy as np
import copy
import gc_control as gc
import gc_calibration as gccu
import pandas as pd
import matplotlib.pyplot as plt




class mixture():
    
    def __init__(self,controller_list):
        self.controllers=controller_list
        self.total_flow=self.calculate_total_flow(self.controllers)
        self.components=self.get_unique_components(self.controllers)
        self.component_fractions=self.get_component_fractions(self.components, self.controllers, self.total_flow)
        self.uncertainty_contributions=self.get_uncertainty_contributions(self.controllers, self.components)
        self.total_flow_contributions=self.get_total_flow_contributions(self.controllers)
        self.controller_dict=self.get_controller_by_id(self.controllers)
        self.partials={}
        for i,specie in enumerate(self.components):
            self.partials[specie]={}
            for h,string in enumerate(self.total_flow_contributions):
                if string not in self.uncertainty_contributions[specie]:
                    self.uncertainty_contributions[specie].append(string)
        self.partials=self.get_partial_derivatives(self.partials,
                                                   self.uncertainty_contributions,
                                                   self.total_flow_contributions,
                                                   self.controller_dict)
        self.sigmas=self.get_sigmas(self.uncertainty_contributions, self.controller_dict)
        self.mixture_uncertainties,self.uncertainty_terms=self.calculate_propogated_uncertainty(self.sigmas, self.partials)
        self.total_flow_uncertainty=0
        for i,spec in enumerate(list(self.mixture_uncertainties.keys())):
            self.total_flow_uncertainty=self.total_flow_uncertainty+(self.mixture_uncertainties[spec])**2.0
            
        self.total_flow_uncertainty=np.sqrt(self.total_flow_uncertainty)
        
    #def add_nominal_X_for_calibration(self,nominal_dict:dict):
    #    self.nominal_X=nominal_dict
    
    def get_residence_time_uncertainty(self,residence_time,volume,Q,sigma_volume,sigma_Q,Tamb=298.15,Pamb=15.0,ReactorTemp=298.15,ReactorPress=15.0,
                                       sigmaTa=0.01,sigmaTr=0.01,sigmaPa=0.01,sigmaPr=0.01,mixing_factor=0.0):
        #term1=(sigma_volume**2.0)/((Q/1000.0)**2.0)
        term1=(sigma_volume**2.0)*((np.divide(ReactorTemp*Pamb,(Q/1000)*ReactorPress*Tamb))**2.0)
        #term2=((sigma_Q/1000.0)**2.0)*(volume**2.0)/((Q/1000)**4.0)
        term2=((sigma_Q/1000.0)**2.0)*((np.divide(-volume*ReactorTemp*Pamb,((Q/1000)**2.0)*ReactorPress*Tamb))**2.0)
        term3=((sigmaTr*ReactorTemp)**2.0)*((np.divide(volume*Pamb,(Q/1000)*ReactorPress*Tamb))**2.0)
        term4=((sigmaPa*Pamb)**2.0)*((np.divide(volume*ReactorTemp,(Q/1000)*ReactorPress*Tamb))**2.0)
        term5=((sigmaPr*ReactorPress)**2.0)*((np.divide(-volume*ReactorTemp*Pamb,(Q/1000)*(ReactorPress**2.0)*Tamb))**2.0)
        term6=((sigmaTa*Tamb)**2.0)*((np.divide(-volume*ReactorTemp*Pamb,(Q/1000)*ReactorPress*(Tamb**2.0)))**2.0)
        term7=(mixing_factor*residence_time)**2.0
        sigma=np.sqrt(term1+term2+term3+term4+term5+term6+term7)
        self.JSR_residence_time_uncertainty=sigma/residence_time
        return sigma/residence_time
    
    def add_IR_measurements(self,file,species,cutoff_time=120,NO_channels=['chemcell','IR'],output_index=1,backpropagate=True):
        
        tempdata=pd.read_csv(file,encoding="ANSI",index_col=False)
        tempdata=tempdata[tempdata['Ellapsed Time(sec)']>cutoff_time]
        IR_averages={}
        IR_stds={}
        
        for i,spec in enumerate(species):
            IR_averages[spec]={}
            IR_stds[spec]={}
            if spec=='NO' and 'chemcell' in NO_channels:
                
                plt.figure()
                plt.plot(tempdata['Ellapsed Time(sec)'],tempdata['NOX (ppm)'])
                plt.title('NOx ppm, condition '+str(output_index))
                a1,b1=np.polyfit(tempdata['Ellapsed Time(sec)'],tempdata['NOX (ppm)'],1)
                plt.xlabel('Time')
                plt.ylabel('ppm')            
                plt.plot(tempdata['Ellapsed Time(sec)'],a1*tempdata['Ellapsed Time(sec)']+b1)
                plt.savefig(spec+'_chemcell_data_id'+str(output_index)+'_unadjusted.pdf',dpi=1200,bbox_inches='tight')
                plt.figure()
                plt.plot(tempdata['Ellapsed Time(sec)'],tempdata['NOX (ppm)']-a1*tempdata['Ellapsed Time(sec)'])
                plt.xlabel('Time')
                plt.ylabel('ppm')
                plt.title('NOx ppm, condition '+str(output_index))
                plt.savefig(spec+'_chemcell_data_id'+str(output_index)+'.pdf',dpi=1200,bbox_inches='tight') 
                if backpropagate:
                    IR_averages[spec]['chemcell']=np.average(tempdata['NOX (ppm)']-a1*tempdata['Ellapsed Time(sec)'])/1000000.0
                    IR_stds[spec]['chemcell']=np.std(tempdata['NOX (ppm)']-a1*tempdata['Ellapsed Time(sec)'])/1000000.0
                elif not backpropagate:
                    IR_averages[spec]['chemcell']=np.average(tempdata['NOX (ppm)'])/1000000.0
                if IR_averages[spec]['chemcell']==0.0:
                    IR_averages[spec]['chemcell']=1e-9
                    
                    
                if not backpropagate:
                    IR_stds[spec]['chemcell']=np.std(tempdata['NOX (ppm)'])/1000000.0
                    dif=np.abs(a1*np.array(tempdata['Ellapsed Time(sec)'])[0]/1000000-a1*np.array(tempdata['Ellapsed Time(sec)'])[-1]/1000000)
                    error=dif
                    IR_stds[spec]['chemcell']=np.sqrt(IR_stds[spec]['chemcell']**2.0+error**2.0)
                
            if spec=='NO2':
                plt.figure()
                plt.plot(tempdata['Ellapsed Time(sec)'],tempdata['HC CONC SPANNED(ppm)'])
                plt.xlabel('Time')
                plt.ylabel('ppm')
                plt.title('NO2 ppm, condition '+str(output_index))
                a2,b2=np.polyfit(tempdata['Ellapsed Time(sec)'],tempdata['HC CONC SPANNED(ppm)'],1)
                plt.plot(tempdata['Ellapsed Time(sec)'],a2*tempdata['Ellapsed Time(sec)']+b2)
                plt.savefig(spec+'_data_id'+str(output_index)+'_unadjusted.pdf',dpi=1200,bbox_inches='tight')
                plt.figure()
                plt.plot(tempdata['Ellapsed Time(sec)'],tempdata['HC CONC SPANNED(ppm)']-a2*tempdata['Ellapsed Time(sec)'])
                plt.xlabel('Time')
                plt.ylabel('ppm')
                plt.title('NO2 ppm, condition '+str(output_index))
                plt.savefig(spec+'_data_id'+str(output_index)+'.pdf',dpi=1200,bbox_inches='tight')
                if backpropagate:
                    IR_averages[spec]['IR']=np.average(tempdata['HC CONC SPANNED(ppm)']-a2*tempdata['Ellapsed Time(sec)'])/1000000.0
                    IR_stds[spec]['IR']=np.std(tempdata['HC CONC SPANNED(ppm)']-a2*tempdata['Ellapsed Time(sec)'])/1000000.0
                else:
                    IR_averages[spec]['IR']=np.average(tempdata['HC CONC SPANNED(ppm)'])/1000000.0
                if IR_averages[spec]['IR']==0.0:
                    IR_averages[spec]['IR']=1e-9
                if not backpropagate:
                    IR_stds[spec]['IR']=np.std(tempdata['HC CONC SPANNED(ppm)'])/1000000.0
                    dif=np.abs(a2*np.array(tempdata['Ellapsed Time(sec)'])[0]/1000000-a2*np.array(tempdata['Ellapsed Time(sec)'])[-1]/1000000)
                    error=dif
                    #print(dif)
                    IR_stds[spec]['IR']=np.sqrt(IR_stds[spec]['IR']**2.0+error**2.0)
                
            if spec=='NO' and 'IR' in NO_channels:
                plt.figure()
                
                plt.plot(tempdata['Ellapsed Time(sec)'],tempdata['CO2 CONC SPANNED(ppm)'])
                plt.xlabel('Time')
                plt.ylabel('ppm')
                plt.title('CO2 ppm, condition '+str(output_index))
                a2,b2=np.polyfit(tempdata['Ellapsed Time(sec)'],tempdata['CO2 CONC SPANNED(ppm)'],1)
                plt.plot(tempdata['Ellapsed Time(sec)'],a2*tempdata['Ellapsed Time(sec)']+b2)
                plt.savefig(spec+'_data_id'+str(output_index)+'_unadjusted.pdf',dpi=1200,bbox_inches='tight')
                plt.figure()
                plt.plot(tempdata['Ellapsed Time(sec)'],tempdata['CO2 CONC SPANNED(ppm)']-a2*tempdata['Ellapsed Time(sec)'])
                plt.xlabel('Time')
                plt.ylabel('ppm')
                plt.title('CO2 ppm, condition '+str(output_index))
                plt.savefig(spec+'_data_id'+str(output_index)+'.pdf',dpi=1200,bbox_inches='tight')
                if backpropagate:
                    IR_averages[spec]['IR']=np.average(tempdata['CO2 CONC SPANNED(ppm)']-a2*tempdata['Ellapsed Time(sec)'])/1000000.0
                    IR_stds[spec]['IR']=np.std(tempdata['CO2 CONC SPANNED(ppm)']-a2*tempdata['Ellapsed Time(sec)'])/1000000.0
                else:
                    IR_averages[spec]['IR']=np.average(tempdata['CO2 CONC SPANNED(ppm)'])/1000000.0
                if IR_averages[spec]['IR']==0.0:
                    IR_averages[spec]['IR']=1e-9
                if not backpropagate:
                    IR_stds[spec]['IR']=np.std(tempdata['CO2 CONC SPANNED(ppm)'])/1000000.0
                    dif=np.abs(a2*np.array(tempdata['Ellapsed Time(sec)'])[0]/1000000-a2*np.array(tempdata['Ellapsed Time(sec)'])[-1]/1000000)
                    error=dif
                    IR_stds[spec]['IR']=np.sqrt(IR_stds[spec]['IR']**2.0+error**2.0)
            
            if spec=='N2O':
                
                plt.figure()
                plt.plot(tempdata['Ellapsed Time(sec)'],tempdata['CO CONC SPANNED(ppm)'])
                plt.xlabel('Time')
                plt.ylabel('ppm')
                plt.title('N2O ppm, condition '+str(output_index))
                a2,b2=np.polyfit(tempdata['Ellapsed Time(sec)'],tempdata['CO CONC SPANNED(ppm)'],1)
                plt.plot(tempdata['Ellapsed Time(sec)'],a2*tempdata['Ellapsed Time(sec)']+b2)
                plt.savefig(spec+'_data_id'+str(output_index)+'_unadjusted.pdf',dpi=1200,bbox_inches='tight')
                plt.figure()
                plt.plot(tempdata['Ellapsed Time(sec)'],tempdata['CO CONC SPANNED(ppm)']-a2*tempdata['Ellapsed Time(sec)'])
                plt.xlabel('Time')
                plt.ylabel('ppm')
                plt.title('N2O ppm, condition '+str(output_index))
                plt.savefig(spec+'_data_id'+str(output_index)+'.pdf',dpi=1200,bbox_inches='tight')
                if backpropagate:
                    IR_averages[spec]['IR']=np.average(tempdata['CO CONC SPANNED(ppm)']-a2*tempdata['Ellapsed Time(sec)'])/1000000.0
                    IR_stds[spec]['IR']=np.std(tempdata['CO CONC SPANNED(ppm)']-a2*tempdata['Ellapsed Time(sec)'])/1000000.0
                else:
                    IR_averages[spec]['IR']=np.average(tempdata['CO CONC SPANNED(ppm)'])/1000000.0
                if IR_averages[spec]['IR']==0.0:
                    IR_averages[spec]['IR']=1e-9
                if not backpropagate:
                    IR_stds[spec]['IR']=np.std(tempdata['CO CONC SPANNED(ppm)'])/1000000.0
                    dif=np.abs(a2*np.array(tempdata['Ellapsed Time(sec)'])[0]/1000000-a2*np.array(tempdata['Ellapsed Time(sec)'])[-1]/1000000)
                    error=dif
                    IR_stds[spec]['IR']=np.sqrt(IR_stds[spec]['IR']**2.0+error**2.0)
                
            self.IR_averages=IR_averages
            self.IR_stds=IR_stds
                
    def get_IR_calibration_uncertainties(self,species,calibration_mixture):
        cal_uncertainties={}
        for i,spec in enumerate(species):
            cal_uncertainties[spec]=calibration_mixture.mixture_uncertainties[spec]/calibration_mixture.component_fractions[spec]
            
        self.IR_cal_uncertainties=cal_uncertainties
        
    def get_combined_exp_averages(self):
        experiment_averages={}
        if hasattr(self, 'IR_averages') and hasattr(self,'species_column_counts'):
            for i, spec in enumerate(list(self.species_column_counts.keys())):
                if spec not in list(experiment_averages.keys()):
                    experiment_averages[spec]={}
                for j,col in enumerate(list(self.species_column_counts[spec].keys())):
                    experiment_averages[spec][col]=self.species_column_counts[spec][col]/100.0
            for i,spec in enumerate(list(self.IR_averages.keys())):
                if spec not in list(experiment_averages.keys()):
                    experiment_averages[spec]={}
                for j,channel in enumerate(list(self.IR_averages[spec].keys())):
                    experiment_averages[spec][channel]=self.IR_averages[spec][channel]
        elif hasattr(self, 'IR_averages') and not hasattr(self,'species_column_counts'):
            experiment_averages=copy.deepcopy(self.IR_averages)
        elif not hasattr(self, 'IR_averages') and hasattr(self,'species_column_counts'):
            experiment_averages=copy.deepcopy(self.species_column_counts)
            for i, spec in enumerate(list(self.species_column_counts.keys())):
                for j,col in enumerate(list(self.species_column_counts[spec].keys())):
                    experiment_averages[spec][col]=experiment_averages[spec][col]/100.0
        
        self.experiment_averages=experiment_averages
        return experiment_averages
    
    def get_total_measurement_uncertainties(self,uncertainty_contributions,caldata={},devices=['GC','IR'],ir_species=[],
                                            precalculated_cal_uncertainties=True,precalculated_noise=True,
                                            precalculated_ir_cal_uncertainties=True,precalculated_ir_noise=True,
                                            NO_IR_channels=['chemcell','IR']):
        sigmas={}
        if 'GC' in devices:
            for i,spec in enumerate(list(caldata.keys())):
                if spec not in list(sigmas.keys()):
                    sigmas[spec]={}
                for j,col in enumerate(list(caldata[spec].keys())):
                    sigmas[spec][col]=0.0
                    
                    for i,name in enumerate(list(uncertainty_contributions[spec][col].keys())):
                        if name=='absolute':
                            sigmas[spec][col]=sigmas[spec][col]+(((uncertainty_contributions[spec][col][name]))**2.0)

                        elif name=='base':
                            if self.species_column_counts[spec][col]/100.0<uncertainty_contributions[spec][col][name]['cutoff']:
                                sigmas[spec][col]=sigmas[spec][col]+(((uncertainty_contributions[spec][col][name]['value']))**2.0)
                            else:
                                pass
                        else:
                            sigmas[spec][col]=sigmas[spec][col]+(((uncertainty_contributions[spec][col][name]/100.0)*(self.species_column_counts[spec][col]/100.0))**2.0)
                        
                    if precalculated_cal_uncertainties:
                        sigmas[spec][col]=sigmas[spec][col]+(self.cal_uncertainties[spec][col]**2.0)
                        
                    if precalculated_noise:
                        sigmas[spec][col]=sigmas[spec][col]+((self.species_column_count_stddevs[spec][col]/100.0)**2.0)
                    sigmas[spec][col]=np.sqrt(sigmas[spec][col])
                    
        if 'IR' in devices and bool(ir_species):
            for i,spec in enumerate(ir_species):
                if spec not in list(sigmas.keys()):
                    sigmas[spec]={}
                for j,channel in enumerate(list(uncertainty_contributions[spec].keys())):
                    if channel in ['chemcell','IR']:
                        sigmas[spec][channel]=0.0
                        for k,name in enumerate(list(uncertainty_contributions[spec][channel].keys())):
                            
                            if name=='absolute':
                                sigmas[spec][channel]=sigmas[spec][channel]+(((uncertainty_contributions[spec][channel][name]))**2.0)
                            elif name=='base':
                                if self.IR_averages[spec][channel]<uncertainty_contributions[spec][channel][name]['cutoff']:
                                    sigmas[spec][channel]=sigmas[spec][channel]+(((uncertainty_contributions[spec][channel][name]['value']))**2.0)
                                else:
                                    pass
                            else:
                                sigmas[spec][channel]=sigmas[spec][channel]+(((uncertainty_contributions[spec][channel][name]/100.0)*(self.IR_averages[spec][channel]))**2.0)
                
                        if precalculated_cal_uncertainties:
                            sigmas[spec][channel]=sigmas[spec][channel]+((self.IR_cal_uncertainties[spec]*self.IR_averages[spec][channel])**2.0)
                            
                            
                        if precalculated_noise:
                        
                            sigmas[spec][channel]=sigmas[spec][channel]+(self.IR_stds[spec][channel]**2.0)
                        sigmas[spec][channel]=np.sqrt(sigmas[spec][channel])
                        
        self.measurement_sigmas=sigmas
    
    
    def add_calibration_nominals(self,nom_dict):
        self.calibration_nominals=nom_dict
        
    def calculate_calibration_sigma(self,area,species,column,caldata):
        term1=((caldata[species][column]['slope-standard-dev']/100.0)**2.0)*(area**2.0)
        term2=((caldata[species][column]['intercept-standard-dev']/100.0)**2.0)
        term3=((caldata[species][column]['slope']/100.0)**2.0)*(self.species_column_area_stddevs[species][column]**2.0)
        print(species,term1,term2,term3)
        sigma=np.sqrt(term1+term2+term3)
        return sigma
    
    def get_calibration_uncertainties(self,caldata):
        cal_uncertainties=copy.deepcopy(self.species_column_areas)
        for i,spec in enumerate(list(caldata.keys())):
            for k,col in enumerate(list(caldata[spec].keys())):
                cal_uncertainties[spec][col]=self.calculate_calibration_sigma(self.species_column_areas[spec][col],spec,col,caldata)
        self.cal_uncertainties=cal_uncertainties
                
    
    def add_gc_measurements(self,gc_file_names, species_columns:dict,ip_address:str):
        #import gc_control as gc
        a=gc.chromatograph(ip=ip_address)
        #a.set_calibration_species(cal_assignments={'A':['N2'],'B':['N2O'],'C':['N2O']})
        a.set_calibration_species(cal_assignments=species_columns)
        species=[]
        #print(species_columns)
        for i,col in enumerate(list(species_columns.keys())):
            for j,spec in enumerate(species_columns[col]):
                if spec not in species:
                    species.append(spec)
        inverted_species_columns={}
        #print(species)
        for i,spec in enumerate(species):
            inverted_species_columns[spec]=[]
            for j,col in enumerate(list(species_columns.keys())):
                if spec in list(species_columns[col]):
                    inverted_species_columns[spec].append(col)
        species_averages={}
        species_areas={}
        #print(inverted_species_columns)
        std_devs={}
        area_std_devs={}
        temp_avs={}
        temp_areas={}
        measurement_averages={}
        measurement_areas={}
        
        species_column_counts={}
        species_column_areas={}
        for k,spec in enumerate(list(inverted_species_columns.keys())):
            temp_avs[spec]=[]
            temp_areas[spec]=[]
            species_column_areas[spec]={}
            species_column_counts[spec]={}
        for i, name in enumerate(gc_file_names):
            a.get_calibrated_species_count(calibration_species=a.species_data,data=a.get_data(name=name))  
            a.get_calibrated_species_areas(calibration_species=a.species_data,data=a.get_data(name=name))
            a.get_species_averages()
            a.get_species_areas()
            species_averages[name]={}
            species_areas[name]={}
            #std_devs[name]={}
            #temp_avs[name]={}
            for k,spec in enumerate(list(inverted_species_columns.keys())):
                species_averages[name][spec]={}
                species_areas[name][spec]={}
            #print(species_averages[name])
            for j,col in enumerate(list(species_columns.keys())):
                for k,spec in enumerate(list(species_columns[col])):
                    #print(col,spec,name)
                    #print(name,spec)
                    #print(species_averages[name][spec])
                    species_averages[name][spec][col]=a.species_averages[col][spec]
                    species_areas[name][spec][col]=a.species_areas[col][spec]
                    temp_areas[spec].append(a.species_areas[col][spec])
                    temp_avs[spec].append(a.species_averages[col][spec])
                    
                    if col in list(species_column_areas[spec].keys()):
                        species_column_areas[spec][col].append(a.species_areas[col][spec])
                        species_column_counts[spec][col].append(a.species_averages[col][spec])
                    else:
                        species_column_areas[spec][col]=[a.species_areas[col][spec]]
                        species_column_counts[spec][col]=[a.species_averages[col][spec]]
            #for k,spec in enumerate(list(inverted_species_columns.keys())):
                #species_averages[name][spec]=np.divide(species_averages[name][spec],float(len(inverted_species_columns[spec])))
        species_column_area_stddevs=copy.deepcopy(species_column_areas)
        species_column_count_stddevs=copy.deepcopy(species_column_counts)
        #print(species_column_counts)
        for i,spec in enumerate(list(species_column_areas.keys())):
            for k,col in enumerate(list(species_column_areas[spec].keys())):
                species_column_area_stddevs[spec][col]=np.std(species_column_areas[spec][col])
                species_column_count_stddevs[spec][col]=np.std(species_column_counts[spec][col])
                species_column_counts[spec][col]=np.mean(species_column_counts[spec][col])
                species_column_areas[spec][col]=np.mean(species_column_areas[spec][col])
        for k,spec in enumerate(list(inverted_species_columns.keys())):
            std_devs[spec]=np.std(temp_avs[spec])
            area_std_devs[spec]=np.std(temp_areas[spec])
            measurement_areas[spec]=np.mean(temp_areas[spec])
            measurement_averages[spec]=np.mean(temp_avs[spec])
        self.measurement_areas=measurement_areas
        self.measurement_averages=measurement_averages
        self.species_averages=species_averages
        self.species_areas=species_areas
        self.std_devs=std_devs
        self.area_std_devs=area_std_devs
        self.species_column_areas=species_column_areas
        self.species_column_counts=species_column_counts
        self.species_column_area_stddevs=species_column_area_stddevs
        self.species_column_count_stddevs=species_column_count_stddevs
        
    def reprocess_gc_data(self,calibration_curve):
        
        std_devs={}
        measurement_averages={}
        temp_avs={}
        for j,spec in enumerate(list(self.measurement_averages.keys())):
            temp_avs[spec]=[]
        
        for i,name in enumerate(list(self.species_areas.keys())):
            for j,spec in enumerate(list(self.species_areas[name].keys())):
                for k,col in enumerate(list(self.species_areas[name][spec].keys())):
                    a=self.species_areas[name][spec][col]
                    newval=0
                    newval=newval+calibration_curve.caldata[spec][col]['slope']*a
                    newval=newval+calibration_curve.caldata[spec][col]['intercept']
                    self.species_averages[name][spec][col]=newval
                    temp_avs[spec].append(newval)
                    
        for j,spec in enumerate(list(self.measurement_averages.keys())):
            std_devs[spec]=np.std(temp_avs[spec])
            measurement_averages[spec]=np.mean(temp_avs[spec])
                    
        self.measurement_averages=measurement_averages
        self.std_devs=std_devs
        
        
    # def add_gc_calibrations(self,gc_file_names, species_columns:dict,ip_address:str):
    #     import gc_control as gc
    #     a=gc.chromatograph(ip=ip_address)
    #     #a.set_calibration_species(cal_assignments={'A':['N2'],'B':['N2O'],'C':['N2O']})
    #     a.set_calibration_species(cal_assignments=species_columns)
    #     species=[]
    #     #print(species_columns)
    #     for i,col in enumerate(list(species_columns.keys())):
    #         for j,spec in enumerate(species_columns[col]):
    #             if spec not in species:
    #                 species.append(spec)
    #     inverted_species_columns={}
    #     #print(species)
    #     for i,spec in enumerate(species):
    #         inverted_species_columns[spec]=[]
    #         for j,col in enumerate(list(species_columns.keys())):
    #             if spec in list(species_columns[col]):
    #                 inverted_species_columns[spec].append(col)
    #     species_averages={}
    #     #print(inverted_species_columns)
    #     for i, name in enumerate(gc_file_names):
    #         a.get_calibration_areas(calibration_species=a.species_data,data=a.get_data(name=name))            
    #         #a.get_species_averages()
    #         a.get_species_areas()
    #         species_averages[name]={}
    #         for k,spec in enumerate(list(inverted_species_columns.keys())):
    #             species_averages[name][spec]=0.0
    #         #print(species_averages[name])
    #         for j,col in enumerate(list(species_columns.keys())):
    #             for k,spec in enumerate(list(species_columns[col])):
    #                 #print(col,spec,name)
    #                 #print(name,spec)
    #                 #print(species_averages[name][spec])
    #                 species_averages[name][spec]=species_averages[name][spec]+a.species_areas[col][spec]
    #         for k,spec in enumerate(list(inverted_species_columns.keys())):
    #             species_averages[name][spec]=np.divide(species_averages[name][spec],float(len(inverted_species_columns[spec])))
    #     self.species_areas=species_averages
    
    def add_gc_calibrations(self,method,species_columns:dict,ip_address):
        a=gc.chromatograph(ip=ip_address)
        a.set_calibration_species(cal_assignments=species_columns)
        calibration_data=a.get_calibration_data(method, species_columns)
        self.calibration_data=calibration_data
        return calibration_data
        
    def calculate_propogated_uncertainty(self,sigmas,partials):
        mixture_uncertainties={}
        uncertainty_terms={}
        for i,specie in enumerate(list(partials.keys())):
            mixture_uncertainties[specie]=0.0
            uncertainty_terms[specie]={}
            for j,string in enumerate(list(partials[specie].keys())):
                mixture_uncertainties[specie]=mixture_uncertainties[specie]+np.multiply(partials[specie][string]**2,
                                                                                        sigmas[specie][string]**2)
                uncertainty_terms[specie][string]=np.multiply(partials[specie][string]**2,sigmas[specie][string]**2)
            mixture_uncertainties[specie]=np.sqrt(mixture_uncertainties[specie])
        return (mixture_uncertainties,uncertainty_terms)
        
        
    def get_sigmas(self,uncertainty_contributions,controller_dict):
        
        sigmas={}
        for i,specie in enumerate(list(uncertainty_contributions.keys())):
            sigmas[specie]={}
            for j,string in enumerate(uncertainty_contributions[specie]):
                if 'X' in string:
                    cid=string.split('d')[-1]
                    #print(string,specie,controller_dict)
                    sigmas[specie][string]=controller_dict[cid].sigma_Xi[specie]
                    
                    
                elif 'Q' in string:
                    cid=string.split('d')[-1]
                    sigmas[specie][string]=controller_dict[cid].sigma_q
        return sigmas
                    
        
    def get_partial_derivatives(self,pds,uncertainty_contributions,total_flow_cont,controllers_dict):
        
        for i,specie in enumerate(list(pds.keys())):
            tempvals={}
            
            for j,string in enumerate(uncertainty_contributions[specie]):
                
                if 'X' in string:
                    #Calculates dX_i/dX_i,j = Q_j/Q_total
                    
                    cid=string.split('d')[-1]
                    tempvals[string]=controllers_dict[cid].flowrate
                    tempvals[string]=np.divide(tempvals[string],self.total_flow)
                    
                elif 'Q' in string:
                    cid=string.split('d')[-1]
                    term1=np.divide(self.get_Xij(cid, specie, controllers_dict),self.total_flow)
                    term2=0.0
                    for k,c in enumerate(list(controllers_dict.keys())):
                        term2=term2+np.multiply(self.get_Xij(c, specie, controllers_dict),controllers_dict[c].flowrate)
                    term2=np.divide(term2,self.total_flow**2.0)
                    tempvals[string]=term1-term2
            pds[specie]=copy.deepcopy(tempvals)
                    
        return pds
    
    def get_Xij(self,cid,specie,controller_dict):
        if specie in list(controller_dict[cid].gas.keys()):
            Xij=controller_dict[cid].gas[specie]
        else:
            Xij=0.0
        return Xij
        
    
    def get_controller_by_id(self,controllers):
        controller_dict={}
        for i,controller in enumerate(controllers):
            controller_dict[str(controller.controller_id)]=controller
        return controller_dict
    
    def calculate_total_flow(self,controllers):
        total=0
        for i,controller in enumerate(controllers):
            total=total+controller.flowrate
        return total
        
    def get_unique_components(self,controllers):
        components=[]
        for i,controller in enumerate(controllers):
            for j,comp in enumerate(list(controller.gas.keys())):
                if comp not in components:
                    components.append(comp)
        return components
    
    def get_component_fractions(self,components,controllers,total_Q):
        
        components=dict.fromkeys(components,0.0)
        for i,controller in enumerate(controllers):
            for j,comp in enumerate(list(controller.gas.keys())):
                components[comp]=components[comp]+np.divide(np.multiply(controller.gas[comp],controller.flowrate),total_Q)
        return components
    
    def get_uncertainty_contributions(self,controllers,components):
        uncertainty_contributions=dict.fromkeys(components,None)
        for i,controller in enumerate(controllers):
            for j,comp in enumerate(list(controller.gas.keys())):
                #if not controller.pure_component:
                    if uncertainty_contributions[comp] is None:
                        uncertainty_contributions[comp]=['X_'+comp+'_id'+str(controller.controller_id)]
                    else:
                        uncertainty_contributions[comp].append('X_'+comp+'_id'+str(controller.controller_id))
                    
        for i,species in enumerate(components):
            for j,controller in enumerate(controllers):
                if species in list(controller.gas.keys()):
                    if uncertainty_contributions[species] is None:
                        uncertainty_contributions[species]=['Q_id'+str(controller.controller_id)]
                    else:
                        uncertainty_contributions[species].append('Q_id'+str(controller.controller_id))
                    print('Adding Q_'+species+' for controller '+str(controller.controller_id)+' to dict for species '+species)
        return uncertainty_contributions
    
    def get_total_flow_contributions(self,controllers):
        contributions=[]
        for i,controller in enumerate(controllers):
            contributions.append('Q_id'+str(controller.controller_id))
        return contributions
        
    #def get_component_uncertainty(self,component):
        
    
    
class controller():
    
    def __init__(self,pure,
                 gas,
                 flowrate,
                 max_range,
                 controller_id,
                 mix_percent_error=0.02,
                 diluent=None,
                 full_scale_error=0.001,
                 control_error=0.005,
                 purity=0.999999):
        #Initialize as follows:
        #set pure to True if mixture is pure component.  Otherwise set false.
        #gas should be a dictionary of species names and their mole fraction
        #flowrate is the total flow rate this controller is set to (L/min)
        #max_range is the largest possible flow rate for this controller with the selected gas mixture (L/min)
        #controller_id is the id given to this controller in LabVIEW
        #Set diluent if this is a mixture: example {'He':1.0} for a mixture with diluent only as helium
        self.pure_component=pure
        self.gas=self.normalize_gas_mixture(gas)
        self.flowrate=flowrate
        self.max_range=max_range
        self.controller_id=controller_id
        self.full_scale_error=full_scale_error
        self.control_error=control_error
        self.mix_percent_error=mix_percent_error
        self.purity=purity
        self.diluent=diluent
        self.sigma_q=self.get_sigma_q(self.max_range, self.flowrate, self.full_scale_error, self.control_error)
        self.sigma_Xi=self.get_sigma_Xi(self.gas,self.mix_percent_error,self.pure_component,self.purity,self.diluent)
        
    def get_sigma_Xi(self,gas,mix_percent_error,pure_component,purity,diluent):
        sigma_Xi={}
        if pure_component:
            sigma_Xi[list(gas.keys())[0]]=np.divide(1.0-purity,1.0)
        else:
            if isinstance(mix_percent_error,dict):
                for i,specie in enumerate(list(gas.keys())):
                    sigma_Xi[specie]=mix_percent_error[specie]*gas[specie]
            else:
                for i,specie in enumerate(list(gas.keys())):
                    sigma_Xi[specie]=mix_percent_error*gas[specie]
        if not pure_component and diluent is not None:
            if len(list(diluent.keys()))==1:
                dil=list(diluent.keys())[0]
                sumhigh=0.0
                sumlow=0.0
                summ=0.0
                for i,spec in enumerate(list(gas.keys())):
                    if dil!=spec:
                        sumhigh=sumhigh+gas[spec]+sigma_Xi[spec]
                        sumlow=sumlow+gas[spec]-sigma_Xi[spec]
                        summ=summ+gas[spec]
                nom=1.0-summ
                high=1.0-sumlow
                low=1.0-sumhigh
                sigma_Xi[list(diluent.keys())[0]]=max(np.abs(nom-high),np.abs(nom-low))
            else:
                pass
        return sigma_Xi
    
    def normalize_gas_mixture(self,gas):
        
        summ=0.0
        for i,item in enumerate(list(gas.keys())):
            summ=summ+gas[item]
        constant=np.divide(1.0,summ)
        for i,item in enumerate(list(gas.keys())):
            gas[item]=np.multiply(gas[item],constant)
        return gas
        
    def get_sigma_q(self,max_range,flowrate,full_scale_error,control_error):
        error=np.multiply(control_error,flowrate)+np.multiply(full_scale_error,max_range)
        return error
        