#%%

import thermo as therm
import matplotlib.pyplot as plt
from scipy import interpolate
import numpy as np
import pandas as pd
# from CoolProp.CoolProp import PropsSI
import cantera as ct
import csv
import os

save_figs = True
save_csv = True

# filedir = 'C:\\Users\\HP USER\\Google Drive\\Burke Group\\JSR Experiments\\New_Pathway_NO\\' # input directory
filedir = r'F:\Google Drive\Burke Group\Mark\DesignOfExperiments\grid\updated_data'
#filedir = r'F:\Google Drive\Burke Group\JSR Experiments\tempfol'
filedir = r'F:\Google Drive\Burke Group\JSR Experiments\N2O_DoE_2022'
filedir = r'C:\Users\mcb22\OneDrive\Documents\JSR_Experiments\JSR_training'
filedir=os.path.abspath(filedir)
conds_file='conditions.csv'
file_headings=['experiment','ratio','temperature','pressure','residence-time']

experiment_name = 'ch4_exp'

model = 'F:\\Google Drive\\Burke Group\\codes\\Mechanisms\\glarborg\\reduced_Glarborg_v1.cti'
model = 'G:\My Drive\Mechanisms\\glarborg\\reduced_Glarborg_v1.cti'
methods =  ['COOLPROP', 'DIPPR_PERRY_8E', 'VDI_PPDS',
'VDI_TABULAR', 'GHARAGHEIZI', 'YOON_THODOS', 'STIEL_THODOS', 'LUCAS_GAS']


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# JSR OPERATIONAL SPACE
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#------------------------------------------------------------------------------
# Dimensions of the Reactor
#------------------------------------------------------------------------------
Reactor = 0.029     # Burke Lab
Nozzle = 300        # Burke Lab 
#Reactor = 0.070    # Test (Matras R1)
#Nozzle = 1000      # Test (Matras R1)

# Ambient Temperature [K]
Tamb = 298.15 

#ResidenceTime  
# tau = [1]        
tau = [1] # input residence times as list

# Reactor Temperature [K]
# reactorTemp = [1200]
# reactorTemp = np.arange(700, 1300, 100) # input temperature as list
# reactorTemp = [700, 750, 800, 825, 850, 875, 900, 925, 950, 1000, 1050, 1100, 1150, 1180]
reactorTemp = [1050]

# Reactor Pascals [Pa]
#reactorPres=1.02069*ct.one_atm

species = ['N2', 'He', 'Ar', 'O2', 'CH4', 'N2O', 'NO', 'H2', 'NH3', 'NO2', 'NH3', 'H2O'] # input species as list

diluent = 'Ar'

conds=pd.read_csv(os.path.join(filedir,conds_file))
listnames=list(conds.columns)
if 'ID' in listnames:
    listnames.remove('ID')
conds=conds[listnames]


# moles_list = [ # calibration Feb 1

    # {'N2O' : 0, 'H2' : 0.00025, 'O2' : 0, 'Ar' : 0.99975},    
    # {'N2O' : 0, 'H2' : 0.0005, 'O2' : 0, 'Ar' : 0.9995},
    # {'N2O' : 0, 'H2' : 0.001, 'O2' : 0, 'Ar' : 0.999},
    # {'N2O' : 0, 'H2' : 0.0025, 'O2' : 0, 'Ar' : 0.9975},
    # {'N2O' : 0, 'H2' : 0.005, 'O2' : 0, 'Ar' : 0.995},
    # {'N2O' : 0, 'H2' : 0.01, 'O2' : 0, 'Ar' : 0.99},

    # {'N2O' : 0, 'H2' : 0, 'O2' : 0.00025, 'Ar' : 0.99975},
    # {'N2O' : 0, 'H2' : 0, 'O2' : 0.0005, 'Ar' : 0.9995},
    # {'N2O' : 0, 'H2' : 0, 'O2' : 0.00075, 'Ar' : 0.99925},
    # {'N2O' : 0, 'H2' : 0, 'O2' : 0.001, 'Ar' : 0.999},
    # {'N2O' : 0, 'H2' : 0, 'O2' : 0.0015, 'Ar' : 0.9985},
    # {'N2O' : 0, 'H2' : 0, 'O2' : 0.002, 'Ar' : 0.998},
    # {'N2O' : 0, 'H2' : 0, 'O2' : 0.004, 'Ar' : 0.996},

    # {'N2O' : 0, 'H2' : 0.01, 'O2' : 0.01, 'Ar' : 0.98},
    # {'N2O' : 0, 'H2' : 0.01, 'O2' : 0.005, 'Ar' : 0.985},
    # {'N2O' : 0, 'H2' : 0.01, 'O2' : 0.0025, 'Ar' : 0.9875},

    # {'N2O' : 0.005, 'H2' : 0, 'O2' : 0, 'Ar' : 0.995}
    
    # {'N2O' : 0, 'H2' : 0.008, 'O2' : 0.008, 'Ar' : 0.984}
    
    # {'N2O' : 0, 'H2' : 0.005, 'O2' : 0.00125, 'Ar' : 0.99375}    
    
    # {'N2O' : 0, 'H2' : 0.005, 'O2' : 0.0005, 'Ar' : 0.9945}    
    
    # ]

# moles_list = [ # experiment Feb 4

#     {'N2O' : 0.005, 'H2' : 0.01, 'O2' : 0.004, 'Ar' : 0.981},     
    
#     ]


# moles_list = [ # experiment Feb 9
           
    
#     # {'O2' : 0.00055},
#     # {'O2' : 0.00080},
#     # {'O2' : 0.00105},
#     # {'O2' : 0.00130},
#     # {'O2' : 0.00155},
#     # {'O2' : 0.00180},
#     # {'O2' : 0.00205}

#     {'N2' : ((1-0.0003947)*0.33) + ((1-0.0004296)*0.33), 'NO2' : (0.0003947) * 0.33, 'H2O' : 0.01, 'NH3' : (0.0004296) * 0.33},
#     {'N2' : ((1-0.0003947)*0.33333333333333333333333333333333) + ((1-0.0004296)*0.333333333333333333333333333333), 'NO2' : (0.0003947) * 0.3333333333333333333333333333, 'H2O' : 0, 'NH3' : (0.0004296) * 0.333333333333333333333333333333333333333333}
#     ]
    
# moles_list = [ # experiment Feb 11

#     {'Ar' : 0.9902, 'N2O' : 0.005, 'H2' : 0.004, 'O2' : 0.0008},     

#     ]
moles_list=[]
reactorTemp=[]
reactorPres=[]
tau=[]
for i,te in enumerate(conds['temperature']):
    tempdic={}
    for n, name in enumerate(conds.columns):
        
        if name not in file_headings:
            tempdic[name]=conds[name][i]
    moles_list.append(tempdic)
    reactorTemp.append(te)
    reactorPres.append(conds['pressure'][i]*ct.one_atm)
    tau.append(conds['residence-time'][i])
    
    

# Adds the appropriate diluent balance to the reactants
# for moles in moles_list:
#     if sum(moles.values()) < 1:
#         moles[diluent] = 1 - sum(moles.values())

#------------------------------------------------------------------------------
# Species Information/Thermo Data
#------------------------------------------------------------------------------
# Stoichiometric Air/Methane (example)
#  (1)CH4 +  (2)(02 + 3.76N2) = (2)H2O + (1)CO2 + (7.52)N2
#------------------------------------------------------------------------------


dat = {}
ID = 1

for i, m in enumerate([*moles_list[0]][:-1]):
    dat[str(m)] = []

dat['T (K)'] = []
dat['T (C)'] = []
dat['P (psi)'] = []
dat['tau (s)'] = []
dat['Q_total (L/min)'] = []

for i, m in enumerate([*moles_list[0]]):
    dat[str('Q_' + m + ' (L/min)')] = []
    
dat['trngl'] = []

visc_list=['N2O','NO2','NO']
def ideal_density(MW,P,T):
    rho=MW*P/(8.314*T)
    return rho

for a, te in enumerate(reactorTemp):
    
            print("Running condition "+str(a))
        
            
            
            #gas = ct.Solution(model)
            #gas.TPX = te, reactorPres[a], moles_list[a]
            
        
            #print('Running for: '+str(moles_list[a])+' at '+str(te)+'K and tau ='+str(tau[a])+'s')    
        
            # Adimensional Constant
            # *(we need to find more data on this constant. thers is supposedly a paper
            # that provides A as a function of temperature) 
            A = 2./np.pi  #0.3                   #np.pi/2 
            A= 0.0022894*te-0.38259
            
            n = np.zeros(len(moles_list[a]))
            M = np.zeros(len(moles_list[a]))
            Cp = np.zeros(len(moles_list[a]))
            rho_amb = np.zeros(len(moles_list[a]))
            rho = np.zeros(len(moles_list[a]))
            nu = np.zeros(len(moles_list[a]))
            
            for i, s in enumerate(species):
                
                vars()[s+'_amb'] = therm.Chemical(s,T=Tamb)
                vars()[s] = therm.Chemical(s,T=te,P=reactorPres[a])
            
            for j, m in enumerate([*moles_list[a]]):
                
                if m in species:
                    
                    # Mole Fractions []
                    n[j] = moles_list[a][m]
                    
                    # Molecular Mass [g/mol]
                    M[j] = vars()[m].MW
                        
                    # Cp [J/(K*g)]
                    Cp[j] = vars()[m].Cpg/1000    
                    
                    # Ambient Mass Density [g/m^3]
                    try:
                        rho_amb[j] = vars()[m+'_amb'].rhog*1000
                    except:
                        rho_amb[j] = ideal_density(vars()[m+'_amb'].MW,reactorPres[a],te)
                    
                    # Mass Density [g/m^3]
                    if m!='NO2' and m!='NO':
                        rho[j] = vars()[m].rhog*1000
                    
                    # Dynamic Viscosity [(g/(m*s^2))*s]
                    if m in visc_list:
                        if m=='N2O':
                            temp_mu=[9.356,10.28,11.3,12.31,13.3,14.28,15.15,15.24,
                                16.19,17.12,23.22,27.18,30.84,34.24,37.41,40.38,40.96]

                            temptemps=[182,200,220,240,260,280,298.15,300,320,340,480,580,680,780,880,980,1000]
                            f=interpolate.interp1d(temptemps,temp_mu,bounds_error=False,fill_value="extrapolate")
                            nu[j]=f(te)/1000000
                        elif m=='NO2':
                            nu[j] = therm.Chemical('Ar',T=te,P=reactorPres[a]).mug
                            rho[j] = therm.Chemical('Ar',T=te,P=reactorPres[a]).rhog
                        elif m=='NO':
                            nu[j] = therm.Chemical('Ar',T=te,P=reactorPres[a]).mug
                            
                    else:
                        
                        nu[j] = vars()[m].mug*1000
                    
                    #print(nu[j])
                    #temp_visc = therm.viscosity.ViscosityGas(m)                    
                    #nu[j] = temp_visc.calculate(te, methods[7]) #* 1000
            
            # Mass Fractions []
            # *(mole fraction below if given mass fraction) 
            y = (n*M)/(sum(n*M))
            #n = (y/M)/sum(y/M)
            
            # Universal Gas Constant [J/(K*mol)]
            R_uni = 8.3144598 
    
            #------------------------------------------------------------------------------
            # Mixture Information/Thermo Data
            #------------------------------------------------------------------------------
            # Molecular Mass [g/mol]
            M_mix = sum(n*M)
            # Gas Constant [J/(K*g)]
            R_mix = (R_uni)/M_mix
            # Mass Density [g/m^3]
            rho_mix_amb = sum(n*rho_amb)
            rho_mix = sum(n*rho)
            # Dynamic Viscosity [g/(m*s)]
            nu_mix = sum(n*nu)
            # nu_mix = gas.viscosity
            # Speed of Sound (T,P) [m/s]
            k_mix = (sum(y*Cp))/((sum(y*Cp)-R_mix))
            c_mix = np.sqrt(k_mix*R_mix*te*1000)
            
            #------------------------------------------------------------------------------
            # Mass/Volumetric Flow Rates
            #------------------------------------------------------------------------------
            ## Volumetric Flow Rate at MFC's [L/min]
            #Qamb = 0.5
            #
            ## Volumetric Flow Rate in Reactor [m^3/s]
            #Q = (Qamb/(1000*60))/(rho_mix/rho_mix_amb)
            #
            ## Mass Flow Rate
            #MFR = Q*rho_mix
            #
            ## Residence Time [s]
            #tau = (np.pi*(Reactor**3))/(3*Q)
            
            # Residence Time [s]
            
            # Volumetric Flow Rate in Reactor [m^3/s]
            Q = 0.000082/tau[a]
            
            # Mass Flow Rate
            MFR = Q*rho_mix 
            
            # Volumetric Flow Rate at MFC's
            # PV=mRT
            Qamb = Q*(rho_mix/rho_mix_amb)
            
            #------------------------------------------------------------------------------
            # The Three Limits
            #------------------------------------------------------------------------------
            # TURBULENCE
            # The jets flowing from the nozzles must be turbulent. This imposes an upper
            # limit of the residence time. Putting this in terms of Reactor Radius (R)
            # and Nozzle Diameter (d) and making it an equality...
            dt=[]
            Rt=[]
            for x in range(1,Nozzle*10000,1):
                x = x*10**(-6)
                dt.append(x)
                Rt.append(((230*tau[a]*nu_mix*x)/(rho_mix*A))**(1/3))
            
            # RECYCLING RATE
            # This is the ratio of residence time to the time it takes to carry all gas in 
            # the reactor from one nozzle to the adjacent one. Experimentally observed that
            # recycling rate must be >30 for good gas-phase mixing. This creates a lower
            # limit on the ratio of Reactor Radius to Nozzle Diameter. Writing in terms of 
            # Reactor Radius (R) and Nozzle Diameter (d) and making it an equality...   
            dr=[]
            Rr=[]
            for x in range(1,Nozzle*10000,1):
                x = x*10**(-6)
                dr.append(x)
                Rr.append((60*x)/(np.pi*A))
                
            # SONIC LIMIT
            # Velocity of the gas exiting the nozzle must be less than the speed of sound
            # at the T and P of the reaction. This will create a lower limit on the 
            # residence time. Putting this in terms of Reactor Radius (R) and Nozzle
            ds=[]
            Rs=[]
            for x in range(1,Nozzle*10000,1):
                x = x*10**(-6)
                ds.append(x)    
                Rs.append(((3*tau[a]*(x**2)*c_mix)/4)**(1/3))
                
            #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            # PLOTTING THE RESULTS
            #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            fig,ax = plt.subplots()
            #------------------------------------------------------------------------------
            # Finding the Corners of the Triangle to Bound an Operational Space
            #------------------------------------------------------------------------------
            point1d = (920*nu_mix)/(3*rho_mix*A*c_mix)
            point1R = ((3*tau[a]*(point1d**2)*c_mix)/4)**(1/3)
            #print(point1d)
            #print(point1R)
            
            point2d = (3*tau[a]*c_mix*((np.pi)**3)*(A**3))/864000
            point2R = (60*point2d)/(np.pi*A)
            #print(point2d)
            #print(point2R)
            
            point3d = ((230*tau[a]*nu_mix*(np.pi**3)*(A**2))/(rho_mix*(60**3)))**(1/2)
            point3R = ((230*tau[a]*nu_mix*point3d)/(rho_mix*A))**(1/3)
            #print(point3d)
            #print(point3R)
            #------------------------------------------------------------------------------
            # Setting Plot Boundaries 
            #------------------------------------------------------------------------------
            ax.set_xlim([point1R - 0.3*point1R,point2R + 0.3*point2R])                                    
            ax.set_ylim([point1d - 0.3*point1d,point2d + 0.3*point2d])                                
            
            #------------------------------------------------------------------------------
            # Setting the aspect ratio of the plot
            #------------------------------------------------------------------------------
            plt.gca().set_aspect(aspect=0.5,adjustable='box')
            
            #------------------------------------------------------------------------------
            # Combining Turbulence and Recycling lines to make shading possible 
            #------------------------------------------------------------------------------
            line4 = np.maximum(Rt,Rr)
            # Bounding this new line using the triangle corners
            line4 = line4[np.logical_and(line4>point1R,line4<point2R)]
            
            #------------------------------------------------------------------------------
            # Shading the triangle and locating the plot text
            #------------------------------------------------------------------------------
            step=[]
            ds_cut=[]
            Rs_cut=[]
            dr_cut=[]
            Rr_cut=[]
            dt_cut=[]
            Rt_cut=[]
            d=0
            
            ax.loglog(Reactor,Nozzle*10**(-6),"o",color='k')
            ax.loglog(Rt,dt,color='b', linestyle='--')
            ax.loglog(Rr,dr,color='r', linestyle='-.')
            ax.loglog(Rs,ds,color='g', linestyle='-')
                
            #------------------------------------------------------------------------------
            plt.title('JSR Operational Limits for '+str(moles_list[a])[:], color='k')
            ax.set_ylabel('Nozzle Diameter [m]', color='k')
            ax.set_xlabel('Reactor Radius [m]', color='k')
            
            # ax = fig.add_subplot(1, 1, 1)
            
            ax.text(0.05, 0.93, r'$\dot{Q}$ (Ln/min) = %.3f'%(Qamb*1000*60*(reactorPres[a]/ct.one_atm)*(273.15/298.15))+' L/min',color='k', size = 8, transform=ax.transAxes, horizontalalignment='left', verticalalignment='center')
            
            ax.text(0.65, 0.55, 'T = %.2f'%te+' K',color='k', size = 8, transform=ax.transAxes, horizontalalignment='left', verticalalignment='center')
            ax.text(0.65, 0.49, r'$\tau$ = %.3f'%tau[a]+' s',color='k', size = 8, transform=ax.transAxes, horizontalalignment='left', verticalalignment='center')
            ax.text(0.65, 0.43, r'$\dot{m}$ = %.3f'%MFR+' g/s',color='k', size = 8, transform=ax.transAxes, horizontalalignment='left', verticalalignment='center')
            ax.text(0.65, 0.37, r'$\dot{Q}$ at MFC = %.3f'%(Qamb*1000*60)+' L/min',color='k', size = 8, transform=ax.transAxes, horizontalalignment='left', verticalalignment='center')
            ax.text(0.65, 0.31, r'$\dot{Q}$ in R = %.3f'%(Q*1000*60)+' L/min',color='k', size = 8, transform=ax.transAxes, horizontalalignment='left', verticalalignment='center')
            
            plt.text(0.22, 0.07, r'$R_{reactor}$= %.3f'%(Reactor)+', $R_{nozzle}$=%.4f'%(Nozzle*10**(-6)),color='k', size = 8, transform=ax.transAxes, horizontalalignment='left', verticalalignment='center')
          
            plt.legend(['Reactor','Turbulence Limit','Recycling Limit','Sonic Limit'],
                       fancybox=True, loc='lower right', prop={'size': 6})
            # Adding tick marks to all sides of the plot 
            ax.tick_params(axis='both',which='both',direction='in',width=1.25,top=True,
                            right=True)
            #plt.show()
            # PLOT DATA
        
            if save_figs:
                fileout = filedir +'\\'+ experiment_name +'_triangle_ID' + str(ID) +'.pdf'
                plt.savefig(fileout,dpi=1200,bbox_inches='tight')
            
            # WRITE DATA
    
            for i, m in enumerate([*moles_list[a]][:-1]):
                dat[str(m)].append(moles_list[a][m]*1000000)   
            
            dat['T (K)'].append(te)
            dat['T (C)'].append(te-273.15)
            dat['P (psi)'].append(np.round(reactorPres[a]*0.000145038,2))
            dat['tau (s)'].append(tau[a])
            dat['Q_total (L/min)'].append(Qamb*1000*60*(reactorPres[a]/ct.one_atm)*(273.15/298.15))
     
            for i, m in enumerate([*moles_list[a]]):
                dat[str('Q_' + m + ' (L/min)')].append(Qamb*1000*60*(reactorPres[a]/ct.one_atm)*(273.15/298.15)*list(moles_list[a].values())[i])  
            
            if Nozzle*10**(-6) < interpolate.interp1d(Rt, dt)(Reactor) and Nozzle*10**(-6) < interpolate.interp1d(Rr, dr)(Reactor) and Nozzle*10**(-6) > interpolate.interp1d(Rs, ds)(Reactor):
                dat['trngl'].append('Y')
            else:
                dat['trngl'].append('N')
                    
            
            ID += 1

df = pd.DataFrame(data=dat)
df.index = np.arange(1, len(df)+1)
if save_csv == True:
    df.to_csv(os.path.join(filedir, experiment_name + '_datasheet.csv'))

#plt.show()
 