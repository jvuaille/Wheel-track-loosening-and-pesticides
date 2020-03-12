# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 08:17:40 2019

@author: cbn978
"""

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from collections import defaultdict
import sys
sys.path.append(r'C:\Users\cbn978\OneDrive - Københavns Universitet\PESTCAST\Git\PyDaisy\pydaisy')
sys.path.insert(0,r'path_to_scripts_folder')
from Functions import weigh_field_areas_sources_stats
from Functions import weigh_field_areas_total_loads


outputs_folder_path = r'C:\Users\cbn978\OneDrive - Københavns Universitet\PhD\Future Cropping\Daisy\sim\parallel_runs\1000y_6_1004_1007'

sowt=[0.7, 0.85] # surface outside WT
siwt=[0.3, 0.15] # surface inside WT
sad=[0.05,0.1] # surface above the drains
sawd=[0.95, 0.90] # surface away from the drains

stat_endpoints=['average','std_error','median','10th_perc','90th_perc']
#chemicals=['glyphosate','phenmedipham', 'metamitron', 'diflufenican',
#           'clomazone', 'fluopyram', 'ethofumesate', 'desamino-metamitron','pyraclostrobin',       
#            'iodosulfuron-methyl-sodium', 'metsulfuron-methyl',  'AE-B107137',
#            'AE-05422291', 'bromoxynil-octanoate', 'bromoxynil', 'dibromo-benzamide',
#             'BAS-500-6', 'BAS-500-7',  'AMPA',  'methyl_carbamate']

chemicals=['glyphosate','phenmedipham','metamitron',
           'clomazone', 'AMPA',  'methyl_carbamate',
           'ethofumesate','desamino-metamitron']

# store the extracted and weighed data 
dict_Pest_sources = defaultdict()
dict_Pest_ally_sources = defaultdict()
dict_Water_sources = defaultdict()
# loop on drain spacing
for icasespace, case_space in enumerate([20,10]):
    # loop on detention capacity
    for threshold in [5,20]:
        # loop on WT area
        for icasetrack, case_track in enumerate([30,15]):
                # water
                field_water = pd.read_csv(outputs_folder_path+'\\'+'global_wflow_sources_field_%d_%d.csv'%(
                        threshold, case_space), index_col = 0)
                drainconnected_water = pd.read_csv(outputs_folder_path+'\\'+'global_wflow_sources_drain_ditch_%d_%d.csv'%(
                        threshold, case_space), index_col = 0)
                # pesticides
                for chemical in chemicals:
                    # 90th percentile load from columns drain connected or not 
                    field_pest = pd.read_csv(outputs_folder_path+'\\'+'global_scores_sources_field_%d_%d.csv'%(
                            threshold, case_space), index_col = 0)
                    drainconnected_pest = pd.read_csv(outputs_folder_path+'\\'+'global_scores_sources_drain_ditch_%d_%d.csv'%(
                            threshold, case_space), index_col = 0)
                    
                    # 300 yearly loads for 3 pesticides
                    if chemical=='glyphosate' or chemical=='phenmedipham' or chemical=='metamitron' :
                        field_pest_ally = pd.read_csv(outputs_folder_path+'\\'+'global_loads_sources_field_%d_%d.csv'%(
                                threshold, case_space), index_col = 0)
                        drainconnected_pest_ally = pd.read_csv(outputs_folder_path+'\\'+'global_loads_sources_drain_ditch_%d_%d.csv'%(
                                threshold, case_space), index_col = 0)
                        
                    # weighed 90th percentile
                    for i,statendpoint in enumerate(stat_endpoints):
                        for source in ['Soil drain', 'Biopore drain', 'Runoff']:
                            dict_Pest_sources[chemical+'_'+source+'_'+statendpoint+'_%d_%d_%d'%(
                                    threshold,case_space,case_track)]=(
                            weigh_field_areas_sources_stats(sowt[icasetrack]*sawd[icasespace],
                                                   siwt[icasetrack]*sawd[icasespace], 
                                                   sowt[icasetrack]*sad[icasespace], 
                                                   siwt[icasetrack]*sad[icasespace], 
                                                   sawd[icasespace], sad[icasespace], 
                                                   field_pest[source+'_'+chemical], 
                                                   drainconnected_pest[source+'_'+chemical])[i])
                        for source in list(drainconnected_water.columns):
                            dict_Water_sources[source+'_'+statendpoint+'_%d_%d_%d'%(
                                    threshold,case_space,case_track)]=(
                            weigh_field_areas_sources_stats(sowt[icasetrack]*sawd[icasespace], 
                                                            siwt[icasetrack]*sawd[icasespace], 
                                                            sowt[icasetrack]*sad[icasespace], 
                                                            siwt[icasetrack]*sad[icasespace], 
                                                            sawd[icasespace], sad[icasespace], 
                                                            field_water[source], 
                                                            drainconnected_water[source])[i])
                    # weighed yearly loads       
                    for source in ['Soil drain', 'Biopore drain', 'Runoff']:
                        if chemical=='glyphosate' or chemical=='metamitron' or chemical=='phenmedipham':
                            dict_Pest_ally_sources[chemical+'_'+source+'_%d_%d_%d'%(
                                    threshold,case_space,case_track)]=(
                            weigh_field_areas_total_loads(sowt[icasetrack]*sawd[icasespace],
                                                   siwt[icasetrack]*sawd[icasespace], 
                                                   sowt[icasetrack]*sad[icasespace], 
                                                   siwt[icasetrack]*sad[icasespace], 
                                                   sawd[icasespace], sad[icasespace], 
                                                   field_pest_ally[source+'_'+chemical], 
                                                   drainconnected_pest_ally[source+'_'+chemical]))
                    

# =============================================================================
# Water balance
# =============================================================================
statistic=  'average'# '90th_perc'
hatch_dict={5:'', 10:"", 20: '///'}
collect,balance=defaultdict(list), defaultdict(list)
for source in list(drainconnected_water.columns):
    fig, ax = plt.subplots(nrows=1, ncols=len([5,20])*4, sharex=True, sharey=True, figsize=(10,5)) 
    fig.text(0.47, 0.04, 'Treatment', fontsize=12)
    fig.text(0.02, 0.5, '%s amount of %s (mm)'%(statistic.title() , source.lower()), va='center', rotation='vertical', fontsize=12) #statistic.split('_')[0],statistic.split('_')[1]  
    my_WT_names=['6-m working width','3-m working width','6-m working width','3-m working width']
    my_x_WT=[0.125,0.34, 0.555, 0.77]
    my_d_names=['20-m drain spacing', '10-m drain spacing']
    my_x_D=[0.215, 0.65]
    for n in range(4):
        fig.text(x=my_x_WT[n],s=my_WT_names[n], y=0.92, fontsize=12)
    for n in range(2):
        fig.text(x=my_x_D[n],s=my_d_names[n], y=0.97, fontsize=12)
    spots=np.arange(0,6)
    n=0
    for icasespace, case_space in enumerate([20,10]): #
        for icasetrack, case_track in enumerate([15,30]):           
            for ithreshold, threshold in enumerate([5,20]): #10,
                source_water=list(dict_Water_sources[source+'_'+statistic+'_%d_%d_%d'%(threshold,case_space,case_track)].values())
                treatments_water=[source_water[0],source_water[1]+source_water[2],
                                  source_water[1]+source_water[3],source_water[1]+source_water[4],
                                  source_water[1]+source_water[5],source_water[1]+source_water[6]]
                
                balance[source+'_%d_%d_%d'%(threshold,case_space,case_track)]=treatments_water

                ax[n].bar(spots, treatments_water, width=0.35, color='lightgray', hatch=hatch_dict[threshold])
                ax[n].set_xticks(np.arange(0,6))
                ax[n].set_xticklabels(['OLOH','OLOH + ILD00','OLOH + ILD08','OLOH + ILD16','OLOH + IHD00','OLOH + IHD08'], fontsize=12, rotation=90)
                ax[n].tick_params(labelsize=12, axis='both', which='major')
                n+=1
    
    plt.ylim(bottom=0)
    plt.subplots_adjust(left=0.07,bottom=0.38, right=0.95,
            top=0.90, wspace=0.25)
    plt.savefig(outputs_folder_path+'\\'+'Plot_%s_%s_thresholds_spacing_WT.png'%(source, statistic))
    my_index=['OLOH','OLOH + ILD00','OLOH + ILD08','OLOH + ILD16','OLOH + IHD00','OLOH + IHD08']
    pd.DataFrame(data=balance, index=my_index).to_csv(outputs_folder_path+'\\WaterBalance_%s.csv'%statistic)  
                

# =============================================================================
# pesticides
# =============================================================================
spraying_rates={'glyphosate': 540, 'phenmedipham': 64*3, 'metamitron': 525*3,
                'fluopyram':125 ,'pyraclostrobin':37.5, 'diflufenican':10,
                'ethofumesate':50*3, 'iodosulfuron-methyl-sodium':2.5, 
           'bromoxynil-octanoate':30, 'AMPA':540, 'methyl_carbamate': 64*3, 'desamino-metamitron': 525*3,
           'AE-B107137':10, 'AE-05422291':10,'metsulfuron-methyl':2.5,'BAS-500-6':37.5, 
           'BAS-500-7':37.5, 'bromoxynil':30, 'dibromo-benzamide':30, 'clomazone':36}



# cumulated percentage of years 

chemical='phenmedipham'  # 'glyphosate'
fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False, figsize=(6,6)) #10,
treatments=['OLOH','OLOH + ILD00','OLOH + ILD08','OLOH + ILD16','OLOH + IHD00','OLOH + IHD08']
mymap=['green','purple','orange','mediumseagreen','crimson','tomato']
mymarker=['o','*','+','^','x','.']

lines=[]
case_space=20
case_track=30 
threshold=5
loads=defaultdict(list)
first=True
for n, source in enumerate(['Soil drain','Biopore drain','Runoff']):
    pest=dict_Pest_ally_sources[chemical+'_%s_%d_%d_%d'%(source,
            threshold,case_space,case_track)]
    for it,t in enumerate(['NW_only','WT_0','WT_8','WT_16','WT_0_HM','WT_8_HM']):
        if first:
            loads[t]=[0 for i in range(len(pest[t]))]
        loads[t]=loads[t]+(pest[t])
    first = False
for it,t in enumerate(['NW_only','WT_0','WT_8','WT_16','WT_0_HM','WT_8_HM']):
    cum_loads=np.unique(loads[t], return_counts=True)
    line1, = ax[0].plot(np.cumsum(cum_loads[1]/np.sum(cum_loads[1])*100),(cum_loads[0]
    /spraying_rates[chemical]*100), marker=mymarker[it],markersize=5, color=mymap[it],label=treatments[it])
    lines.append(line1)
    ax[1].plot(np.cumsum(cum_loads[1]/np.sum(cum_loads[1])*100),(cum_loads[0]
    /spraying_rates[chemical]*100), marker=mymarker[it],markersize=5, color=mymap[it],label=treatments[it])
    lines.append(line1)
    ax[0].tick_params(labelsize=15, axis='both', which='major')    
    ax[1].tick_params(labelsize=15, axis='both', which='major')  
ax[0].grid(True)
ax[1].grid(True)
#ax[0].set_ylim(bottom=0.1, top=2)
#ax[1].set_ylim(bottom=0, top=0.05)
ax[0].set_ylim(bottom=1, top=5)
ax[1].set_ylim(bottom=0, top=0.1)
ax[0].spines['bottom'].set_visible(False)
ax[1].spines['top'].set_visible(False)
ax[0].xaxis.tick_top()
ax[0].tick_params(labeltop=False)
ax[1].xaxis.tick_bottom()
d = .01
kwargs = dict(transform=ax[0].transAxes, color='k', clip_on=False)
ax[0].plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax[0].plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
kwargs.update(transform=ax[1].transAxes)  
ax[1].plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax[1].plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
fig.text(0.01, 0.5,'Load of %s'%chemical+' (% applied)',fontsize=15,va='center', rotation='vertical')
plt.xlabel('Cumulated % of years',fontsize=15)
#ax[0].set_yticks([0.5,1,1.5,2]) 
#ax[1].set_yticks([0,0.01,0.02,0.03,0.04,0.05]) 
plt.xlim(left=50, right=101)
plt.xticks([50,75,90,100]) 
ax[0].legend(fontsize=15, loc='upper left')
plt.subplots_adjust(hspace=0.2, left=0.17, right=0.96, top=0.96)
plt.savefig(outputs_folder_path+'\\'+'Plot_5_20_30_cumulated_%s_3origins.png'%chemical)


chemical='metamitron'
fig, ax = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=False, figsize=(6,6)) #10,
treatments=['OLOH','OLOH + ILD00','OLOH + ILD08','OLOH + ILD16','OLOH + IHD00','OLOH + IHD08']
mymap=['green','purple','orange','mediumseagreen','crimson','tomato']
mymarker=['o','*','+','^','x','.']
lines=[]
case_space=20
case_track=30 
threshold=5
loads=defaultdict(list)
first=True
for n, source in enumerate(['Soil drain','Biopore drain','Runoff']):
    
    pest=dict_Pest_ally_sources[chemical+'_%s_%d_%d_%d'%(source,
            threshold,case_space,case_track)]
    for it,t in enumerate(['NW_only','WT_0','WT_8','WT_16','WT_0_HM','WT_8_HM']):
        if first:
            loads[t]=[0 for i in range(len(pest[t]))]
        loads[t]=loads[t]+(pest[t])
    first = False
for it,t in enumerate(['NW_only','WT_0','WT_8','WT_16','WT_0_HM','WT_8_HM']):
    cum_loads=np.unique(loads[t], return_counts=True)
    line1, = ax.plot(np.cumsum(cum_loads[1]/np.sum(cum_loads[1])*100),(cum_loads[0]
    /spraying_rates[chemical]*100), marker=mymarker[it],markersize=5, color=mymap[it],label=treatments[it])
    lines.append(line1)
    ax.tick_params(labelsize=15, axis='both', which='major')  
ax.grid(True)
ax.set_ylim(bottom=0, top=1.6)
fig.text(0.01, 0.5,'Load of %s'%chemical+' (% applied)',fontsize=15,va='center', rotation='vertical')
plt.xlabel('Cumulated % of years',fontsize=15)
plt.xlim(left=50, right=101)
plt.xticks([50,75,90,100]) 
ax.legend(fontsize=15, loc='upper left')
plt.subplots_adjust(hspace=0.2, left=0.17, right=0.96, top=0.96)
plt.savefig(outputs_folder_path+'\\'+'Plot_5_20_30_cumulated_%s_3origins.png'%chemical)


# =============================================================================
# 90th percentile loads of glyphosate, phenmedipham and metamitron
# reference setting
# =============================================================================
chemicals=['glyphosate', 'phenmedipham', 'metamitron']
statistic='90th_perc'
hatch_dict={5:'', 10:"", 20: '/'}

fig, ax = plt.subplots(nrows=1, ncols=3, sharex=True, sharey=False, figsize=(9,7)) #10,
fig.text(0.48, 0.09, 'Treatment', fontsize=15)
fig.text(0.01, 0.65, '$\mathregular{%s^{%s}}$ %s. load of compound ('%(statistic.split('_')[0].split('th')[0],
         statistic.split('_')[0].split('90')[1],statistic.split('_')[1])+"% applied)", 
va='center', rotation='vertical', fontsize=15) 

my_x_thresh=np.arange(0.12, 0.9, 0.105)
spots=np.arange(0,6)
n=0
chemical_name=0.95
chemical_name_x=[0.18, 0.455, 0.775]
lines=[]
for ichemical,chemical in enumerate(chemicals):
    fig.text(x=chemical_name_x[ichemical],s=chemical.title(), y=chemical_name, fontsize=15) 
    for icase_space, case_space in enumerate([20]):   
        for icasetrack, case_track in enumerate([30]): 
            for ithreshold, threshold in enumerate([5]): 
                matrix_pest=list(dict_Pest_sources[chemical+'_Soil drain_'+statistic+'_%d_%d_%d'%(threshold,case_space,case_track)].values())
                treatments_matrix=[matrix_pest[0],matrix_pest[1]+matrix_pest[2],
                                 matrix_pest[1]+matrix_pest[3],matrix_pest[1]+
                                 matrix_pest[4],matrix_pest[1]+matrix_pest[5],
                                 matrix_pest[1]+matrix_pest[6]]
                treatments_matrix_perc_applied=np.array(treatments_matrix)/spraying_rates[chemical]*100
                biop_pest=list(dict_Pest_sources[chemical+'_Biopore drain_'+statistic+'_%d_%d_%d'%(threshold,case_space,case_track)].values())
                treatments_biop_pest=[biop_pest[0],biop_pest[1]+biop_pest[2],
                                   biop_pest[1]+biop_pest[3],biop_pest[1]+
                                   biop_pest[4],biop_pest[1]+biop_pest[5],
                                   biop_pest[1]+biop_pest[6]]
                treatments_biop_pest_perc_applied=np.array(treatments_biop_pest)/spraying_rates[chemical]*100
                runo_pest=list(dict_Pest_sources[chemical+'_Runoff_'+statistic+'_%d_%d_%d'%(threshold,case_space,case_track)].values())
                treatments_runo_pest=[runo_pest[0],runo_pest[1]+runo_pest[2],
                                runo_pest[1]+runo_pest[3],runo_pest[1]+
                                runo_pest[4],runo_pest[1]+runo_pest[5],
                                runo_pest[1]+runo_pest[6]]
                treatments_runo_pest_perc_applied=np.array(treatments_runo_pest)/spraying_rates[chemical]*100
                line1 = ax[n].bar(spots, treatments_matrix_perc_applied, width=0.4, 
                  color='darkcyan', label='Matrix', hatch=hatch_dict[threshold])
                line2 = ax[n].bar(spots, treatments_biop_pest_perc_applied, width=0.4, 
                  color='silver', label='Biopores', hatch=hatch_dict[threshold], bottom=treatments_matrix_perc_applied)
                line3 = ax[n].bar(spots, treatments_runo_pest_perc_applied, width=0.4, 
                  color='crimson', label='Runoff', hatch=hatch_dict[threshold], bottom=treatments_biop_pest_perc_applied+treatments_matrix_perc_applied)
                  
                ax[n].tick_params(labelsize=15, axis='both', which='major')
                lines.append(line1)
                lines.append(line2)
                lines.append(line3)
                
                ax[n].set_ylim(bottom=0, top=[0.05, 0.3, 2][n])
                ax[n].set_yticks([[0.0, 0.025, 0.05],[0.00, 0.15,0.30],[0,1,2]][n]) 
                ax[n].set_xticks(np.arange(0,6))
                ax[n].set_xticklabels(['OLOH','OLOH + ILD00','OLOH + ILD08','OLOH + ILD16','OLOH + IHD00','OLOH + IHD08'], rotation=90, fontsize=15)
              
    n+=1
plt.subplots_adjust(left=0.14,bottom=0.38, right=0.95,
        top=0.92, wspace=0.35, hspace=0.5)
fig.legend((lines[0], lines[1], lines[2]),('Matrix drainage','Biopore drainage','Surface runoff'),fontsize=15,loc=(0.12,0.010), ncol=3) 
plt.savefig(outputs_folder_path+'\\'+'Plot_5_20_30_%s_%s_%s_%s_3origins.png'%(statistic, chemicals[0],chemicals[1],chemicals[2])) 



# =============================================================================
# 90th percentile loads of glyphosate, phenmedipham and metamitron
# sensitivity analysis
# =============================================================================
chemicals=['glyphosate', 'phenmedipham',  'metamitron']
#            ,'fluopyram','diflufenican',, 'ethofumesate','clomazone'
#             ,'pyraclostrobin','AMPA',  'methyl_carbamate','desamino-metamitron'
#            'iodosulfuron-methyl-sodium', 'metsulfuron-methyl',  'AE-B107137',
#            'AE-05422291', 'bromoxynil-octanoate', 'bromoxynil', 'dibromo-benzamide',
#             'BAS-500-6', 'BAS-500-7',  ]

statistic='90th_perc'
hatch_dict={5:'', 10:"", 20: ''}

risk_inc_WT, risk_dec_loose, total_loads, runoff = defaultdict(list), defaultdict(list), defaultdict(list), defaultdict(list)
fig, ax = plt.subplots(nrows=3, ncols=3, sharex=True, sharey='row', figsize=(7,10)) #10,
fig.text(0.48, 0.06, 'Treatment', fontsize=12)
fig.text(0.02, 0.5, '$\mathregular{%s^{%s}}$ %s. load of compound ('%(statistic.split('_')[0].split('th')[0],
         statistic.split('_')[0].split('90')[1],statistic.split('_')[1])+"% applied)", va='center', rotation='vertical', fontsize=12) #(g $\mathregular{ha^{-1}}$)
my_WT_names=['6-m working width','3-m working width','3-m working width']
my_x_WT=[0.15,0.44, 0.72]
my_d_names=['20-m drain spacing', '10-m drain spacing']
my_x_D=[0.3, 0.72]
for n in range(3):
    fig.text(x=my_x_WT[n],s=my_WT_names[n], y=0.95, fontsize=12)
for n in range(2):
    fig.text(x=my_x_D[n],s=my_d_names[n], y=0.98, fontsize=12)
spots=np.arange(0,6)
n=0
chemical_name=[0.91,0.67,0.43]
chemical_name_x=[0.48, 0.46,0.48] 
lines=[]
for ichemical,chemical in enumerate(chemicals):
    fig.text(x=chemical_name_x[ichemical], y=chemical_name[ichemical],s=chemical.title(), fontsize=12) 
    l=0
    for icase_space, case_space in enumerate([20, 10]): 
        if case_space == 20:
            ww=[15,30]
        else:
            ww=[30]
        for icasetrack, case_track in enumerate(ww): 
            for ithreshold, threshold in enumerate([20]): 
                matrix_pest=list(dict_Pest_sources[chemical+'_Soil drain_'+statistic+'_%d_%d_%d'%(threshold,case_space,case_track)].values())
                treatments_matrix=[matrix_pest[0],matrix_pest[1]+matrix_pest[2],
                                 matrix_pest[1]+matrix_pest[3],matrix_pest[1]+
                                 matrix_pest[4],matrix_pest[1]+matrix_pest[5],
                                 matrix_pest[1]+matrix_pest[6]]
                treatments_matrix_perc_applied=np.array(treatments_matrix)/spraying_rates[chemical]*100
                
                biop_pest=list(dict_Pest_sources[chemical+'_Biopore drain_'+statistic+'_%d_%d_%d'%(threshold,case_space,case_track)].values())
                treatments_biop_pest=[biop_pest[0],biop_pest[1]+biop_pest[2],
                                   biop_pest[1]+biop_pest[3],biop_pest[1]+
                                   biop_pest[4],biop_pest[1]+biop_pest[5],
                                   biop_pest[1]+biop_pest[6]]
                treatments_biop_pest_perc_applied=np.array(treatments_biop_pest)/spraying_rates[chemical]*100
 
                runo_pest=list(dict_Pest_sources[chemical+'_Runoff_'+statistic+'_%d_%d_%d'%(threshold,case_space,case_track)].values())
                treatments_runo_pest=[runo_pest[0],runo_pest[1]+runo_pest[2],
                                runo_pest[1]+runo_pest[3],runo_pest[1]+
                                runo_pest[4],runo_pest[1]+runo_pest[5],
                                runo_pest[1]+runo_pest[6]]
                treatments_runo_pest_perc_applied=np.array(treatments_runo_pest)/spraying_rates[chemical]*100
#                ax[n,l].grid(True, linestyle='--')
                line1 = ax[n,l].bar(spots, treatments_matrix_perc_applied, width=0.4, 
                  color='darkcyan', label='Matrix', hatch=hatch_dict[threshold])
                line2 = ax[n,l].bar(spots, treatments_biop_pest_perc_applied, width=0.4, 
                  color='silver', label='Biopores', hatch=hatch_dict[threshold], bottom=treatments_matrix_perc_applied)
                line3 = ax[n,l].bar(spots, treatments_runo_pest_perc_applied, width=0.4, 
                  color='crimson', label='Runoff', hatch=hatch_dict[threshold], 
                  bottom=treatments_biop_pest_perc_applied+treatments_matrix_perc_applied)
                  
                ax[n,l].tick_params(labelsize=12, axis='both', which='major')
                lines.append(line1)
                lines.append(line2)
                lines.append(line3)
                
                # risk increased by wheeling for matrix, biopore drain flow and runoff
                ## risk from total 90th perc. load
                risk_inc_WT[chemical+'_%d_%d_%d_total'%(threshold,case_space,case_track)]=((
                        treatments_biop_pest_perc_applied[1::3]+treatments_matrix_perc_applied[1::3]+
                        treatments_runo_pest_perc_applied[1::3])/(treatments_biop_pest_perc_applied[0]+
                        treatments_matrix_perc_applied[0]+treatments_runo_pest_perc_applied[0])-1)*100
                
                ## risk from runoff
                risk_inc_WT[chemical+'_%d_%d_%d_runoff'%(threshold,case_space,case_track)]=(
                           treatments_runo_pest_perc_applied[1::3]/
                           treatments_runo_pest_perc_applied[0]-1)*100  
                        
                # risk decreased by loosening for matrix, biopore drain flow and runoff        
                ## risk from total 90th perc. load
                risk_dec_loose[chemical+'_%d_%d_%d_total'%(threshold,case_space,case_track)]=(1-(
                        treatments_biop_pest_perc_applied[np.array([2,3,5])]+
                        treatments_matrix_perc_applied[np.array([2,3,5])]+
                        treatments_runo_pest_perc_applied[np.array([2,3,5])])/(
                        treatments_biop_pest_perc_applied[np.array([1,1,4])]+
                        treatments_matrix_perc_applied[np.array([1,1,4])]+
                        treatments_runo_pest_perc_applied[np.array([1,1,4])]))*100
                        
                ## risk from runoff
                risk_dec_loose[chemical+'_%d_%d_%d_runoff'%(threshold,case_space,case_track)]=(1-
                           treatments_runo_pest_perc_applied[np.array([2,3,5])]/
                           treatments_runo_pest_perc_applied[np.array([1,1,4])])*100            
                
                # total loads for all scenarios
                total_loads[chemical+'_%d_%d_%d'%(threshold,case_space,case_track)]=(treatments_biop_pest_perc_applied+
                        treatments_matrix_perc_applied+
                        treatments_runo_pest_perc_applied)
                
                # runoff loads for all scenarios
                runoff[chemical+'_%d_%d_%d'%(threshold,case_space,case_track)]=treatments_runo_pest_perc_applied
                
                ax[n,l].set_ylim(bottom=0, top=[0.05, 0.3, 2][n]) # [2,1]  [0.002, 0.008, 0.2]
                ax[n,l].set_yticks([[0.0, 0.025, 0.05],[0.00, 0.15,0.3],[0,1,2]][n]) #  [0,1,2],[0,0.5,1]  [0.0, 0.001, 0.002],[0.00, 0.004,0.008],[0,0.1,0.2]
                if n==2:
                    ax[n,l].set_xticks(np.arange(0,6))
                    ax[n,l].set_xticklabels(['OLOH','OLOH + ILD00','OLOH + ILD08','OLOH + ILD16','OLOH + IHD00','OLOH + IHD08'], rotation=90, fontsize=12)
                  
            l+=1
    n+=1

plt.subplots_adjust(left=0.15,bottom=0.22, right=0.95,
        top=0.90, wspace=0.2, hspace=0.2)
fig.legend((lines[0], lines[1], lines[2]),('Matrix drainage','Biopore drainage','Surface runoff'), loc=(0.09,0.015),fontsize=12, ncol=3)
plt.savefig(outputs_folder_path+'\\'+'Plot_20_sensitivity_%s_%s_%s_%s_3origins.png'%(statistic, chemicals[0],chemicals[1],chemicals[2]))  #chemical

my_index=['OLOH + ILD00','OLOH + IHD00']
pd.DataFrame(data=risk_inc_WT, index=my_index).to_csv(outputs_folder_path+'\\risk_increase_20_WT_%s_%s_%s_%s.csv'%(statistic, chemicals[0],chemicals[1],chemicals[2]))  
my_index=['OLOH + ILD08','OLOH + ILD16','OLOH + IHD08']
pd.DataFrame(data=risk_dec_loose, index=my_index).to_csv(outputs_folder_path+'\\risk_decrease_loos_20_%s_%s_%s_%s.csv'%(statistic, chemicals[0],chemicals[1],chemicals[2]))   
my_index=['OLOH','OLOH + ILD00','OLOH + ILD08','OLOH + ILD16','OLOH + IHD00','OLOH + IHD08']
pd.DataFrame(data=total_loads, index=my_index).to_csv(outputs_folder_path+'\\total_loads_treatments_%s_%s_%s_%s.csv'%(statistic, chemicals[0],chemicals[1],chemicals[2]))   
pd.DataFrame(data=runoff, index=my_index).to_csv(outputs_folder_path+'\\runoff_loads_treatments_%s_%s_%s_%s.csv'%(statistic, chemicals[0],chemicals[1],chemicals[2]))   




#statistic='90th_perc'
#hatch_dict={5:'', 10:"", 20: ''}
#risk=defaultdict()
##chemicals=['fluopyram','pyraclostrobin', 'diflufenican','AE-B107137', 'AE-05422291',
##           'iodosulfuron-methyl-sodium', 'metsulfuron-methyl','BAS-500-6', 'BAS-500-7',
##           'bromoxynil-octanoate', 'bromoxynil', 'dibromo-benzamide']
#for chemical in chemicals:
#    fig, ax = plt.subplots(nrows=1, ncols=len([5,20])*4, sharex=True, sharey=True, figsize=(10,5)) #10,
#    fig.text(0.47, 0.04, 'Treatment', fontsize=12)
#    fig.text(0.02, 0.5, '%s %s loss of compound (g $\mathregular{ha^{-1}}$)'%(statistic.split('_')[0],statistic.split('_')[1]), va='center', rotation='vertical', fontsize=12) #
#    spots=np.arange(0,6)
#    for icasetrack, case_track in enumerate([15,30]):
#        for icasedrain, case_drain in enumerate([5,10]):            
#            for ithreshold, threshold in enumerate([5,10,20]): #
#                field=[list(dict_Pest[chemical+'_'+statistic+'_%d_%d_%d'%(threshold,case_track,case_drain)].values())[k] for k in range(7)]
#                field_loosening=[field[0],field[1]+field[2],field[1]+field[3],field[1]+field[4],field[1]+field[5],field[1]+field[6]]
#                field_loosening_perc_applied=np.array(field_loosening)/spraying_rates[chemical]*100
#                drainco=[list(dict_Pest[chemical+'_'+statistic+'_%d_%d_%d'%(threshold,case_track,case_drain)].values())[k] for k in np.arange(7,14)]
#                drainco_loosening=[drainco[0],drainco[1]+drainco[2],drainco[1]+drainco[3],drainco[1]+drainco[4],drainco[1]+drainco[5],drainco[1]+drainco[6]]
#                drainco_loosening_perc_applied=np.array(drainco_loosening)/spraying_rates[chemical]*100
#                
#                ax[ithreshold+icasedrain*2+icasetrack*2].bar(spots, field_loosening_perc_applied, width=0.35, color='lightgray', label='Away from drain lines', hatch=hatch_dict[threshold])
#                ax[ithreshold+icasedrain*2+icasetrack*2].bar(spots, drainco_loosening_perc_applied, width=0.35, color='darkgray', label='Above drain lines', hatch=hatch_dict[threshold], bottom=field_loosening)
#                ax[ithreshold+icasedrain*2+icasetrack*2].set_xticks(np.arange(0,6))
#                ax[ithreshold+icasedrain*2+icasetrack*2].set_xticklabels(['OLOH','OLOH + ILD00','OLOH + ILD08','OLOH + ILD16','OLOH + IHD00','OLOH + IHD08'], rotation=90)
#                ax[ithreshold+icasedrain*2+icasetrack*2].tick_params(labelsize=12, axis='both', which='major')
#                
#                loosening=np.array(field_loosening)+np.array(drainco_loosening)
#                risk[chemical+'_%d_%d_%d'%(case_drain,case_track,threshold)]=(1-loosening[1:]/loosening[0])*100 # percentage improvement
#
#                fig.text(x=[0.21,0.625][icasetrack], y=0.85, s='%d '%(threshold)+'mm', fontsize=12)
#            fig.text(x=[0.21,0.625][icasetrack], y=0.90, s='Portion of wheel tracks: %d '%(case_track)+'%', fontsize=12)
#        fig.text(x=[0.21,0.625][icasedrain], y=0.95, s='Portion above drain lines: %d '%(case_drain)+'%', fontsize=12)
#        plt.ylim(bottom=0, top=100)
#        fig.suptitle(t='Active ingredient: %s'%chemical, y=0.99, x=0.5)
#        plt.subplots_adjust(left=0.15,bottom=0.38, right=0.95,
#                top=0.80, wspace=0.25, hspace=0.05)
#        plt.savefig(outputs_folder_path+'\\'+'Plot_%s_%s_above_drain_area_5_and_10_drains.png'%(statistic, chemical))
#        my_index=['OLOH + ILD00','OLOH + ILD08','OLOH + ILD16','OLOH + IHD00','OLOH + IHD08']
#        pd.DataFrame(data=risk, index=my_index).to_csv(outputs_folder_path+'\\risk_improvement_%s.csv'%statistic)  
#         


                
#statistic='median'
#hatch_dict={5:'', 10:"", 20: ''}
#fig, ax = plt.subplots(nrows=1, ncols=len([5,10,20])*2, sharex=True, sharey=True, figsize=(10,5)) #
#fig.text(0.47, 0.04, 'Treatment', fontsize=12)
#fig.text(0.02, 0.5, '%s amount of water transported by runoff\, drain-connected biopores or matrix drainage (mm $\mathregular{y^{-1}}$)'%(statistic.title()), va='center', rotation='vertical', fontsize=12)
#spots=np.arange(0,6)
#for icasedrain, case_drain in enumerate([5,10]):     
#    for ithreshold, threshold in enumerate([5,10,20]): #
#        field=[list(dict_Water[statistic+'_%d_%d'%(threshold,case_drain)].values())[k] for k in range(7)]
#        field_loosening=[field[0],field[1]+field[2],field[1]+field[3],field[1]+field[4],field[1]+field[5],field[1]+field[6]]
##        field_stde=[list(dict_Water['std_error'+'_%d_%d'%(threshold,case_drain)].values())[k] for k in range(7)]
#        drainco=[list(dict_Water[statistic+'_%d_%d'%(threshold,case_drain)].values())[k] for k in np.arange(7,14)]
#        drainco_loosening=[drainco[0],drainco[1]+drainco[2],drainco[1]+drainco[3],drainco[1]+drainco[4],drainco[1]+drainco[5],drainco[1]+drainco[6]]
##        drainco_stde=[list(dict_Water['std_error'+'_%d_%d'%(threshold,case_drain)].values())[k] for k in np.arange(7,14)]
#        ax[ithreshold+icasedrain*2+icasetrack*2].bar(spots, field_loosening, width=0.35, color='lightgray', label='Away from drain lines', hatch=hatch_dict[threshold])
##        ax[ithreshold].errorbar(spots, field, field_stde, linestyle='None', marker='^',color='black')
#
#        ax[ithreshold+icasedrain*2+icasetrack*2].bar(spots, drainco_loosening, width=0.35, color='darkgray', label='Above drain lines', hatch=hatch_dict[threshold], bottom=field_loosening)
##        ax[ithreshold].errorbar(spots, np.array(field)+np.array(drainco), drainco_stde, linestyle='None', marker='^',color='black')
#        ax[ithreshold+icasedrain*2+icasetrack*2].set_xticks(np.arange(0,6))
#        ax[ithreshold+icasedrain*2+icasetrack*2].set_xticklabels(['OL,OH','OL,OH + ILD00','OL,OH + ILD08','OL,OH + ILD16','OL,OH + IHD00','OL,OH + IHD08'], fontsize=12, rotation=90)
#        ax[ithreshold+icasedrain*2+icasetrack*2].tick_params(labelsize=12, axis='both', which='major')
#    fig.text(x=[0.21,0.625][icasedrain], y=0.95, s='Portion above drain lines: %d '%([5,10][icasedrain])+'%', fontsize=12)
#
#plt.ylim(bottom=0)
#plt.subplots_adjust(left=0.12,bottom=0.38, right=0.95,
#        top=0.9, wspace=0.25)
#plt.savefig(outputs_folder_path+'\\'+'Plot_%s_water_balance_above_drain_area_5_and_10_drains.png'%(statistic))
#


