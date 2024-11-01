# Script for post-processing results from STOMP Wallula chemistry model
# to feed to PEST for parameter estimation, modified from Katie Muller's
# analysis scripts

# Katie added (1) Fake Secondary Mineral Observations to nudge optimization towards desired secondary mineral ppts;
# (2) adjusted obs weights to be normalized obs weights (i.e., sum to 1) for improved model goodness of fit comparisions;
# (3) Error check to see if simulation made it to the final time,
# if not writes all simulation outputs to be 100 to send PEST away from that parameter space.

import pandas as pd
from pathlib import Path
import numpy as np
import subprocess
import os

# Modify these species to match observations written in PEST input file.

# List Aqueous Species for Comparision
# aq_species = ['fe', 'mn', 'ca', 'mg', 'k', 'na', 'si', 'hco3', 'pH']
aq_species = ['fe', 'mn', 'ca', 'mg', 'k', 'na', 'si', 'pH']
discount_species = ['fe', 'k', 'na', 'si']

# List of Minerals in Simulation (needs to be all inculsive list for aqueous chemistry calculation to be correct) 
minerals = ['anatase', 'ankerite', 'aragonite', 'beidellite-ca', 'beidellite-k','beidellite-mg', 'calcite', 'chalcedony', \
            'clinopyroxene', 'dawsonite', 'dolomite', 'glass', 'magnetite','plagioclase', 'rhodochrosite', 'siderite']

# List STOMP nodes where observation data is available
# OFT = 1,1,9 = node 625
# IRFT1 = 1,1,13 = node 937
# IRFT2 = 1,1,20 = node 1483
nodes = [625, 937, 1483]

p = subprocess.Popen('out2csv output.csv', shell=True, stdin=subprocess.PIPE, 
                      stdout=subprocess.PIPE, universal_newlines=True)
newline = os.linesep # [1]
commands = ['a', 'a']
p.communicate( newline.join( commands))

# File Names
f_simulation_output = 'output.csv'

# Simulation Pressure Variable Name
sim_pressure_name = ['Pressure Differential for Well #3, MPa']
pressure_obs_name = 'Walulla-Post-Injection-Hydraulic-Test-Field-Data.csv'
chem_obs_name = 'Water_Chem.xlsx'

# Directories
model_dir = Path('./')
data_dir = Path('./')

# DataFrames
sim_output = pd.read_csv(model_dir / f_simulation_output)
sim_col_names = sim_output.keys()
chem_obs = pd.read_excel(data_dir / chem_obs_name)
pressure_obs = pd.read_csv(data_dir / pressure_obs_name)

# Hydraulic Testing Times:
x_min = 705.6806/365.25 #yr (start of test)
x_max = 709.3576/365.25 #day (end of test)

# Check if simulation completed
if sim_output.TM.iloc[-1] < x_max:
    sim_error = True
    print('Simulation did not reach required time for PEST comparision')
else:
    sim_error = False
####

sim_col_names = sim_output.keys()

# Loop over simulation columns and extract column headers for each node and aq_species for comparision to observations
sim_output_no_min = sim_output
for j, name in enumerate(minerals):
        for element in sim_col_names:
            try:
                index = element.index(name)
                sim_output_no_min = sim_output_no_min.drop(columns = element)
            except ValueError:
                pass
sim_col_names_no_min = sim_output_no_min.keys()

model_output = []
observations = []
observation_name = []
weights = []
node_weight = []
species_weight = []

for i, node in enumerate(nodes):
    res = []
    for j, name in enumerate(aq_species):
        res = []
        for element in sim_col_names_no_min:
            try:
                index = element.index(name)
                if element.index(str(node)):
                    res.append(element)
            except ValueError:
                pass
        
        # Sum species to get the total component concentration
        total_conc = sim_output_no_min[res].sum(axis=1)

        filtered_obs = np.where((chem_obs['STOMP_Node'] == node) & (chem_obs['Aq_Species'] == '{}'.format(name)))
        filtered_obs = chem_obs.loc[filtered_obs]
        
        observations.append(np.log(filtered_obs['Conc [M]'][filtered_obs['Time [yr]']>0]).tolist())
        
        for each_time in enumerate(np.log(filtered_obs['Conc [M]'][filtered_obs['Time [yr]']>0]).tolist()):
            observation_name.append('{} {} {}'.format(node, name, each_time))
        
        model_output.append(np.log(np.interp(filtered_obs['Time [yr]'][filtered_obs['Time [yr]']>0],sim_output['TM'],total_conc)).tolist())
        obs_weight = 1.0 / np.abs(np.sum(np.log(filtered_obs['Conc [M]'][filtered_obs['Time [yr]']>0]).tolist()))
        if node == 625:
            node_weight = 0.3
        else:
            node_weight = 1.0
        if element in discount_species:
            species_weight = 0.25
        else:
            species_weight = 1.0
        weight = np.ones(len(np.log(filtered_obs['Conc [M]'][filtered_obs['Time [yr]']>0]).tolist()))*obs_weight*node_weight*species_weight
        weights.append(weight.tolist())

# Secondary Minerals 
res = []
secondary_mineral_name = []
mols_secondary_minerals = []    

for j, name in enumerate(minerals):
    res = []
    for element in sim_col_names:
        try:
            index = element.index(name)
            if element.index(str(node)):
                res.append(element)
        except ValueError:
            pass        
    intial_conc = sim_output.iloc[:1]     

    if (intial_conc[res].iloc[0] == 0).any():
        secondary_mineral_name.append(name)
        mols = sim_output.loc[(sim_output.TM == 1.93205)][res].values[0]
        mols_secondary_minerals.append(mols)
        
normalized_secondary = mols_secondary_minerals/sum(mols_secondary_minerals)
#####    
# Secondary minerals -- proxy observation points
proxy_secondary_obs = [0.05, 0.22, 0.25, 0.05, 0.05, 0.05, 0, 0.05, 0, 0, 0.05, 0.23]

for i, item in enumerate(normalized_secondary):
    observation_name.append('normalized_secondary_'+ str(secondary_mineral_name[i]))
    observations.append(item)

wt = np.ones(len(normalized_secondary))
weights.append(wt.tolist())    
        
# Now do the pressure data
# Convert from simulation output of psi to MPa
psi_to_MPa = 0.00689476

x_min = 705.6806/365.25 #yr
initial_pressure_value = np.interp(x_min, sim_output['TM'],sim_output['PW3CW-1'])
sim_pressure_diff_MPa = sim_output['PW3CW-1'].subtract(initial_pressure_value, axis = 0)*psi_to_MPa
sim_time_day = sim_output['TM']*365.25
obs_time_day = pressure_obs['Full Time [yr]']*365.25
x_min = 705.6806 #day
x_max = 709.3576 #day

pressure_obs = pressure_obs.drop(index=0)
observations.append(pressure_obs['Post Injection Pressure Differential [Mpa]'].tolist())

for item in pressure_obs['Post Injection Pressure Differential [Mpa]']:
    observation_name.append('pressure')

model_output.append(np.interp(pressure_obs['Full Time [yr]'],sim_output['TM'],sim_pressure_diff_MPa).tolist())
obs_weight = 1.0 / (np.abs(np.sum(pressure_obs['Post Injection Pressure Differential [Mpa]'])))
node_weight = 1.0
weight = np.ones(len(pressure_obs['Post Injection Pressure Differential [Mpa]'].tolist())) * obs_weight * node_weight
weights.append(weight.tolist())

observations = np.concatenate(observations)
model_output = np.concatenate(model_output)
weights = np.concatenate(weights)

with open("observation_data_processed.txt", "w") as txt_file:
    for concentration in observations:
        txt_file.write('%.20f' % concentration + '\n')
i = 0
with open("trial.ins","w") as txt_file:
    txt_file.write('pif @' + '\n')
    for concentration in observations:
        i += 1
        txt_file.write('l1 [obs'+str(i) + ']1:22' + '\n')

with open("trial.pst","w") as txt_file:
    block = """pcf
* control data
norestart  estimation
# Number of parameters, number of observations, number of parameter groups, number of prior information, number of observation groups
    2     """ + str(i) + """     1     0     1
# Number of pairs of input template files, number of model output reading instruction files, 
    1     1 single point   1   0   0
  5.0   2.0   0.3  0.03    10
  3.0   3.0 0.001  0
  0.1
   30  0.01     3     3  0.01     3
    1     1     1
* parameter groups
dgr         relative 0.01  0.0  switch  2.0 parabolic
* parameter data
glassrate          none  relative    7.17e-07      7.0e-8   7.0e-6 dgr             1.0000        0.0000      1
clinorate          none  relative    4.13e-05      4.0e-6   4.0e-4 dgr             1.0000        0.0000      1
* observation groups
obsgroup
* observation data
# Obs Name  Observation        Weight    Group
"""
    i = 0
    txt_file.write(block)
    
    # Normalize all weights 
    total_weight = weights.sum()
    normalized_wt = weights/total_weight
    
    for concentration in observations:
        i += 1
        txt_file.write('obs'+str(i) + '   %.20f' % concentration + '    '+ str(normalized_wt[i-1]) +'    obsgroup' + '\n')
    block = """* model command line
./trial.sh
* model input/output
input.tpl  input
trial.ins wallula_output_processed.txt
* prior information"""
    txt_file.write(block)

with open("wallula_output_processed.txt", "w") as txt_file:
    if sim_error:
        model_output_error = np.ones(len(model_output.tolist()))*100
        for concentration in model_output_error:
            txt_file.write('   %.20f' % concentration + '\n')
    else:
        for concentration in model_output:
            txt_file.write('   %.20f' % concentration + '\n')
        
i=0
with open("obs_with_names.txt", "w") as txt_file:
    for concentration in observations:
        i += 1
        txt_file.write('obs_name :' + str(observation_name[i-1]) + 'obs'+ str(i) + '   %.20f' % concentration + '    '+'wt:' + str(weights[i-1])+'normalized_wt:' + str(normalized_wt[i-1]) + '\n')