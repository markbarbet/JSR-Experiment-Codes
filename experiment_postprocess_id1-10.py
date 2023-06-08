# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 17:41:19 2022

@author: Mark Barbet
"""

import flow_mixture_uncertainties as fmu
import gc_calibration_uncertainties as gcu
import pandas as pd
import os
import gc_calibration as gcc
import matplotlib.pyplot as plt
import yaml

diluent='He'
Tamb=293.15
mixing_factor=0.05
output_csv_dir_string='/home/mcb/MSI/data/ID1-10/'
jsr_yaml_template='jsr_template.yaml'
experiment_id_start=1
experiment_file_subdir='N2O_DoE_Oct_2022\\MSI\\error_checking'
devices=['GC','IR']
os.chdir(r'C:\Users\mcb22\OneDrive\Documents\JSR_Experiments')
exp_conds='s_14_experiment_conditions.xlsx'
ip='10.10.0.1'
method='DoE_N2O_2022_M2'
calconds='GC_calibration_uncertainty_conditions.xlsx'
ircalconds='IR_calibration_conditions_uncertainty.xlsx'
montecarlo_plotnames='montecarlo_output'
gc_calibration_dict={'A':['NO','O2','N2'],'C':['N2O']}
ir_species=['NO','NO2']
ir_no_channels=['chemcell']
ir_subdir='N2O_DoE_Oct_2022\\Final_IR_ID1-17'
JSR_volume=8.5e-5
JSR_volume_stddev=2e-6
pressure=1.02069 #atm
pressure_uncertainty=0.01
temperatures=[1050,1100,1050,1050,1050,1100,1050,1050,1100,1050]
temperature_uncertainties=[0.01]*len(temperatures)
residence_times=[0.45,1.5,0.45,0.45,0.45,1.5,0.45,0.45,1.5,0.45]
uncertainty_conts={'NO':{'A':{'drift':1.0},'chemcell':{'linearity':1.0,'absolute':0.000005}},
                   'N2':{'A':{'drift':1.0,'extra':3.0}},
                   'O2':{'A':{'drift':1.0,'extra':3.0}},
                   'N2O':{'C':{'drift':1.0,'extra':5.0}},
                   'NO2':{'IR':{'linearity':1.0,'absolute':0.00005}}}
exp_names=[['DoE_exp_Oct6_id1_run1','DoE_exp_Oct6_id1_run2','DoE_exp_Oct6_id1_run6','DoE_exp_Oct6_id1_run8'],
           ['DoE_exp_Oct11_id2_run1','DoE_exp_Oct11_id2_run2','DoE_exp_Oct11_id2_run3'],
           ['DoE_exp_Oct6_id3_run2','DoE_exp_Oct6_id3_run3'],
           ['DoE_exp_Oct6_id4_run1','DoE_exp_Oct6_id4_run2','DoE_exp_Oct6_id4_run3'],
           ['DoE_exp_Oct6_id5_run2','DoE_exp_Oct6_id5_run3'],
           ['DoE_exp_Oct11_id6_run2','DoE_exp_Oct11_id6_run3'],
           ['DoE_exp_Oct6_id7_run1','DoE_exp_Oct6_id7_run2','DoE_exp_Oct6_id7_run3'],
           ['DoE_exp_Oct6_id8_run2','DoE_exp_Oct6_id8_run4'],
           ['DoE_exp_Oct11_id9_run2','DoE_exp_Oct11_id9_run4'],
           ['DoE_exp_Oct6_id10_run2','DoE_exp_Oct6_id10_run3']]

IR_datasets=["id1 - REDO.csv",
           'id2.csv',
           'id3.csv',
           'id4.csv',
           'id5.csv',
           'id6.csv',
           'id7.csv',
           'id8.csv',
           'id9.csv',
           'id10.csv']

#residence_times=[]
if 'GC' in devices:
    caldatafile=pd.ExcelFile(os.path.join(os.getcwd(),calconds))
    calsheet1=pd.read_excel(caldatafile,'1')
    cal_length=len(calsheet1['ControllerID'])
    calmixtures=[]
if 'IR' in devices:
    ircaldatafile=pd.ExcelFile(os.path.join(os.getcwd(),ircalconds))
    ircalsheet1=pd.read_excel(ircaldatafile,'1')
    ircal_length=len(ircalsheet1['ControllerID'])
    ircalmixtures=[]
datafile=pd.ExcelFile(os.path.join(os.getcwd(),exp_conds))
sheet1=pd.read_excel(datafile,'1')
exp_length=len(sheet1['ControllerID'])
mixtures=[]
exp_length=len(sheet1['ControllerID'])


caldata=[]
data=[]
coldata=[]
std_devs=[]
cal_curves=[]
if 'GC' in devices:
    caldict=gc_calibration_dict
    
    for i in range(cal_length):
        cal_controllers=[]
        nominal_concentrations={}
        for j,sheet in enumerate(caldatafile.sheet_names):
            tempdata=pd.read_excel(caldatafile, str(j+1))
            calconditions={}
            for k,name in enumerate(tempdata.columns):
                if name not in ['ControllerID','Pure','Max_Q','Q'] and 'Nominal' not in name:
                    calconditions[name]=tempdata[name][i]
            if tempdata['Q'][i]!=0.0:
                cal_controllers.append(fmu.controller(pure=tempdata['Pure'][i],gas=calconditions, flowrate=tempdata['Q'][i], max_range=tempdata['Max_Q'][i], controller_id=tempdata['ControllerID'][i]))
                #print(tempdata['Pure'][i],conditions,tempdata['Q'][i],tempdata['Max_Q'][i],tempdata['ControllerID'][i])
            for k,name in enumerate(tempdata.columns):
                if 'Nominal' in name:
                    nominal_concentrations[name.split('_')[-1]]=tempdata[name][i]
        calmix=fmu.mixture(cal_controllers)
        calmix.add_calibration_nominals(nominal_concentrations)
        calmixtures.append(calmix)
        
if 'IR' in devices:
    #caldict={'A':['NO','O2','N2'],'C':['N2O']}
    
    for i in range(ircal_length):
        ircal_controllers=[]
        nominal_concentrations={}
        for j,sheet in enumerate(ircaldatafile.sheet_names):
            tempdata=pd.read_excel(ircaldatafile, str(j+1))
            ircalconditions={}
            for k,name in enumerate(tempdata.columns):
                if name not in ['ControllerID','Pure','Max_Q','Q'] and 'Nominal' not in name:
                    ircalconditions[name]=tempdata[name][i]
            if tempdata['Q'][i]!=0.0:
                ircal_controllers.append(fmu.controller(pure=tempdata['Pure'][i],gas=ircalconditions, flowrate=tempdata['Q'][i], max_range=tempdata['Max_Q'][i], controller_id=tempdata['ControllerID'][i]))
                #print(tempdata['Pure'][i],conditions,tempdata['Q'][i],tempdata['Max_Q'][i],tempdata['ControllerID'][i])
            for k,name in enumerate(tempdata.columns):
                if 'Nominal' in name:
                    nominal_concentrations[name.split('_')[-1]]=tempdata[name][i]
        ircalmix=fmu.mixture(ircal_controllers)
        ircalmix.add_calibration_nominals(nominal_concentrations)
        ircalmixtures.append(ircalmix)

for i in range(exp_length):
    controllers=[]
    for j,sheet in enumerate(datafile.sheet_names):
        tempdata=pd.read_excel(datafile, str(j+1))
        conditions={}
        for k,name in enumerate(tempdata.columns):
            if name not in ['ControllerID','Pure','Max_Q','Q','ResidenceTime','Temperature','Pressure'] and 'Nominal' not in name:
                conditions[name]=tempdata[name][i]
        if tempdata['Q'][i]!=0.0:
            controllers.append(fmu.controller(pure=tempdata['Pure'][i],gas=conditions, flowrate=tempdata['Q'][i], max_range=tempdata['Max_Q'][i], controller_id=tempdata['ControllerID'][i]))
            #print(tempdata['Pure'][i],conditions,tempdata['Q'][i],tempdata['Max_Q'][i],tempdata['ControllerID'][i])
    mix=fmu.mixture(controllers)
    caldata.append(mix.add_gc_calibrations(method, caldict, ip))
    mix.add_gc_measurements(exp_names[i], caldict, ip)
    #coldata.append(mix.species_averages)
    #data.append(mix.measurement_averages)
    #std_devs.append(mix.std_devs)
    mixtures.append(mix)
    mixtures[-1].get_residence_time_uncertainty(residence_times[i],JSR_volume, mixtures[-1].total_flow, JSR_volume_stddev,
                                                mixtures[-1].total_flow_uncertainty,Tamb=Tamb,ReactorTemp=temperatures[i],
                                                sigmaTr=temperature_uncertainties[i],mixing_factor=mixing_factor)
    #residence_times.append(JSR_volume/(mixtures[-1].total_flow/1000.0))
    if 'GC' in devices and 'IR' not in devices:
        cal_curves.append(gcc.calibration_curve(mix))
        mixtures[-1].reprocess_gc_data(cal_curves[-1])
        coldata.append(mix.species_averages)
        data.append(mix.measurement_averages)
        std_devs.append(mix.std_devs)
        cal_curves[-1].add_calibration_mixtures(calmixtures)
        if i==0:
            cal_curves[-1].monte_carlo_fitting(cal_curves[-1].caldata)
            cal_curves[-1].plot_montecarlo_distributions(montecarlo_plotnames)
        else:
            cal_curves[-1].set_new_caldata(cal_curves[0].caldata)
        mixtures[-1].get_calibration_uncertainties(cal_curves[0].caldata)
        mixtures[-1].get_total_measurement_uncertainties(uncertainty_conts,caldata=cal_curves[-1].caldata,devices=devices)
        
    if 'GC' not in devices and 'IR' in devices:
        #cal_curves.append(gcc.calibration_curve(mix))
        #mixtures[-1].reprocess_gc_data(cal_curves[-1])
        #coldata.append(mix.species_averages)
        #data.append(mix.measurement_averages)
        std_devs.append(mix.std_devs)
        #cal_curves[-1].add_calibration_mixtures(calmixtures)
        tempfilename=os.path.join(os.path.join(os.getcwd(),ir_subdir),IR_datasets[i])
        mixtures[-1].add_IR_measurements(tempfilename, ir_species,NO_channels=ir_no_channels,output_index=i+1)
        #mixtures[-1].add_IR_measurements(IR_datasets[i], ir_species,NO_channels=ir_no_channels,output_index=i+1)
        mixtures[-1].get_IR_calibration_uncertainties(ir_species, ircalmixtures[i])
        #mixtures[-1].get_calibration_uncertainties(cal_curves[0].caldata)
        mixtures[-1].get_total_measurement_uncertainties(uncertainty_conts,devices=devices,ir_species=ir_species,NO_IR_channels=ir_no_channels)
    if 'GC' in devices and 'IR' in devices:
        cal_curves.append(gcc.calibration_curve(mix))
        mixtures[-1].reprocess_gc_data(cal_curves[-1])
        coldata.append(mix.species_averages)
        data.append(mix.measurement_averages)
        std_devs.append(mix.std_devs)
        cal_curves[-1].add_calibration_mixtures(calmixtures)
        if i==0:
            cal_curves[-1].monte_carlo_fitting(cal_curves[-1].caldata)
            cal_curves[-1].plot_montecarlo_distributions(montecarlo_plotnames)
        else:
            cal_curves[-1].set_new_caldata(cal_curves[0].caldata)
        mixtures[-1].get_calibration_uncertainties(cal_curves[0].caldata)
        tempfilename=os.path.join(os.path.join(os.getcwd(),ir_subdir),IR_datasets[i])
        mixtures[-1].add_IR_measurements(tempfilename, ir_species,NO_channels=ir_no_channels,output_index=i+1)
        mixtures[-1].get_IR_calibration_uncertainties(ir_species, ircalmixtures[i])
        mixtures[-1].get_total_measurement_uncertainties(uncertainty_conts,caldata=cal_curves[-1].caldata,devices=devices,
                                                         ir_species=ir_species,NO_IR_channels=ir_no_channels)
    mixtures[-1].get_combined_exp_averages()
with open(os.path.join(os.getcwd(),jsr_yaml_template)) as f:
            template = yaml.load(f,Loader=yaml.FullLoader)
            
            for i,mix in enumerate(mixtures):
            
                template['common-properties']['temperature']['value-list']=[float(round(temperatures[i],9))]
                template['common-properties']['temperature']['relative-uncertainty']=float(temperature_uncertainties[i])
                template['common-properties']['pressure']['value']=float(round(pressure,9))
                template['common-properties']['pressure']['relative-uncertainty']=float(pressure_uncertainty)
                #print(template['apparatus'])
                template['apparatus']['residence-time']['value']=float(round(residence_times[i],9))
                template['apparatus']['residence-time']['relative-uncertainty']=float(mixtures[i].JSR_residence_time_uncertainty)
                template['common-properties']['composition']=[]
                mole_sum=0
                for j,species in enumerate(list(mixtures[i].component_fractions.keys())): 
                    template['common-properties']['composition'].append({'species':species,
                                                                         'mole-fraction':float(round(mixtures[i].component_fractions[species],9)),
                                                                         'relative-uncertainty':float(round(mixtures[i].mixture_uncertainties[species]/mixtures[i].component_fractions[species],9))})
                    #mole_sum=mole_sum+conditions[species][i]
                    #Above value of relative uncertainty is not necessary for this procedure, leave unchanged
                #diluent_value=1.0-mole_sum
                #print(diluent_value)
                #template['common-properties']['composition'].append({'species':self.input_options['diluent'],
                #                                                     'mole-fraction':float(round(diluent_value,9)),
                #                                                     'relative-uncertainty':0.05})
                template['common-properties']['assumptions']['thermal-boundary']='isothermal'
                template['common-properties']['assumptions']['mechanical-boundary']='constant pressure'
                template['apparatus']['volume']=float(round(JSR_volume,9))
                template['file-author']={'name':'Burke Lab JSR post-processing: Mark Barbet'}
                template['apparatus']['kind']='JSR'
                template['apparatus']['facility']='Computer Simulated'
                template['apparatus'].pop('institution',None)
                template['common-properties']['assumptions']['equation-of-state']='ideal gas'
                
                #outfilenames.append(self.input_options['yaml_output_name']+str(i+1)+'.yaml')
                
                #templates=[{}]
                #templatescount=0
                #used_species_apps={}
                number_of_yamls=0
                for j,species in enumerate(list(mixtures[i].experiment_averages.keys())):
                    for k,app in enumerate(list(mixtures[i].experiment_averages[species].keys())):
                        if k+1>number_of_yamls:
                            number_of_yamls=k+1
                tempdicts=[[]]
                for count in range(number_of_yamls):
                    template['datapoints']={'mole-fraction':[]}  
                    
                    for j,species in enumerate(list(mixtures[i].experiment_averages.keys())):
                        #if species not in templates[templatescount]:
                            #templates[templatescount][species]=''
                            
                        if len(list(mixtures[i].experiment_averages[species].keys()))<count+1:
                            pass
                        else:
                            template['datapoints']['mole-fraction'].append({'csvfile':os.path.join(output_csv_dir_string,
                                                                            'experiment_id'+str(i+experiment_id_start)+'_'+species+'_'+str(count)+'.csv'),
                                                                            'targets':[{'name':species,
                                                                                        'species':species,
                                                                                        'absolute-uncertainty':0.0,
                                                                                        'relative-uncertainty':0.05}]})
                        if not os.path.exists(os.path.join(os.getcwd(),experiment_file_subdir)):
                            os.makedirs(os.path.join(os.getcwd(),experiment_file_subdir))
                            
                        for k,app in enumerate(list(mixtures[i].experiment_averages[species].keys())):
                            with open(os.path.join(os.path.join(os.getcwd(),experiment_file_subdir),
                                               'experiment_id'+str(i+experiment_id_start)+'_'+species+'_'+str(k)+'.csv'),'w') as f:
                            #lines=['Temperature,'+species+',Relative_Uncertainty\n']
                            
                            
                                lines=['Temperature,'+species+',Relative_Uncertainty\n']
                                lines.append(str(round(temperatures[i],9))+','+str(round(mixtures[i].experiment_averages[species][app],9))+','+str(round(mixtures[i].measurement_sigmas[species][app]/mixtures[i].experiment_averages[species][app],9))+'\n')
                                lines[-1]=lines[-1].rstrip('\n')
                                f.writelines(lines)
                                
                    template['datapoints']['concentration']=[]
                    template['datapoints']['concentration'].append({'csvfile':None,
                                                             'targets':[{'name':None,
                                                                         'species':None,
                                                                         'absolute-uncertainty':None,
                                                                         'relative-uncertainty':None}]})
                    
                    with open(os.path.join(os.path.join(os.getcwd(),experiment_file_subdir),'experiment_id'+str(i+experiment_id_start)+'_'+str(count)+'.yaml'),'w') as f:
                            yaml.safe_dump(template, f,default_flow_style=False)