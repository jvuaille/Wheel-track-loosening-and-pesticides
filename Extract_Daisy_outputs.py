# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 11:19:44 2019

@author: cbn978
"""

from collections import defaultdict
import numpy as np
import pandas as pd
from datetime import datetime

path_folder = r''
import sys
sys.path.append(r'path_to_pydaisy')
from pydaisy.Daisy import DaisyDlf
sys.path.insert(0,r'path_to_scripts_folder')
from Functions import extract
from Functions import extract_inDict_sum_variable_sources
from Functions import stats_on_source_list
from Functions import stats_on_sum_list, load_source_list

## pesticides
chemicals=['glyphosate','phenmedipham','metamitron',
           'clomazone', 'AMPA',  'methyl_carbamate',
           'ethofumesate','desamino-metamitron']

## Daisy outputs for pesticides and water
sources=['Leak matrix','Leak biopores', 'Soil drain', 'Biopore drain','Runoff']
wsources=['Precipitation','Actual evapotranspiration','Matrix percolation',
          'Biopore percolation','Matrix drain flow','Biopore drain flow','Runoff']

## Daisy run in parallel: 3 runs per setting
# Each run starts at a different date
my_start_dates={0:1966, 1: 2299, 2:2632} 
# different settings for different detention capacity and drain spacing
sensitivity_analysis={'5_10':['1_5_10m','2_5_10m','3_5_10m'], 
                      '5_20':['1_5_20m','2_5_20m','3_5_20m'], 
                      '20_10':['1_20_10m','2_20_10m','3_20_10m'],
                      '20_20':['1_20_20m','2_20_20m','3_20_20m']}

folder_name = 'myfoldername'
all_spraying_dates_check, wet_spraying_dates = [], [] 

# years to remove: glyphosate sprayed on wet soil
wet_spraying_dates=[2047, 2113, 2119, 2143, 2179, 2209, 2215, 2221, 2302, 2332,
                    2338, 2359, 2362, 2386, 2398, 2410, 2416, 2446, 2476, 2506,
                    2536, 2542, 2548, 2590, 2701, 2713, 2749, 2785, 2803, 2815,
                    2821, 2893, 2929, 2947, 2953]
prefix= 'myprefixfilename'

# loop on the column type: drain connected or not
for col in ['field','drain']: 
    # loop on drain spacing
    for space in [10,20]: 
        # loop on detention capacity
        for threshold in [5,20]: 
            # store the 90 th percentile loads of pesticides and water
            global_scores_sources = defaultdict(list)
            global_wflow_sources = defaultdict(list)
            # store the yearly loads of pesticides
            global_loads_sources=defaultdict(list)
            # loop on treatments
            for depth in ['NW','D0','D8','D16','D0_HM','D8_HM']: #['OLOH','ILD00','ILD08','ILD16','IHD00','IHD08']
                print(threshold, space, col, depth)
                # store the yearly load of pesticides and water
                yearly_leaching_sources, yearly_wflow_sources = defaultdict(list), defaultdict(list)
                # find simulated data in corresponding folder
                for my_start_year,folder in enumerate(sensitivity_analysis['%d'%threshold+'_'+'%d'%space]): 
                    
                    print('extraction from folder nb:',my_start_year)
                    # store the pesticide outputs from Daisy
                    field_chemicals=defaultdict()
                    for chem in chemicals:
                        field_chemicals[chem]=DaisyDlf(path_folder+'\\'+folder+
                                       '\\'+'%s-%s_%s_field_%s.dlf'%(prefix, 
                                       depth, col, chem)).Data 
                        field_chemicals[chem]['Date']=field_chemicals[chem].index
                        field_chemicals[chem].index=np.arange(0,len(field_chemicals[chem]))  
                        
                    # use the harvest file to select find starting and ending 
                    # index of 3-month monitoring period each year
                    harvest=DaisyDlf(path_folder+'\\'+folder+'\\'+'%s-harvest.dlf'%prefix)
                    harvest=harvest.Data
                    harvest['Date']=harvest.index
                    harvest.index=np.arange(0,len(harvest))
                    start_harvest=(np.where(harvest['Date'] == datetime.strptime(
                    '%d-%02d-%02d 00:00:00'%(my_start_dates[my_start_year], 9, 
                    25), "%Y-%m-%d %H:%M:%S"))[0][0])
                    # store the starting and ending index every year,
                    # and corresponding water outputs
                    starting_periods, ending_periods = defaultdict(), defaultdict()
                    field_water, starting_periods[depth], ending_periods[depth]=extract(
                    path_folder,folder+'\\'+'%s-%s_%s_field_water'%(prefix, depth,col), 4, 10, 7, 10) 
                    # list of successive crops (on treatment column) from harvest file
                    index_crops = (harvest['column'][start_harvest:][harvest['column']
                    [start_harvest:] == depth+'_'+col].index)
                    crops=harvest['crop'][start_harvest:].loc[index_crops]
                    print('Nb years:',len(crops))
                    # loop over the successive crops
                    for i in range(len(crops)):
                        field_water_period = field_water[starting_periods[depth][i]:ending_periods[depth][i]+1]
                        if not field_water_period['Date'].values[0].year in wet_spraying_dates:
                            if crops.values[i] == 'Sugar Beet':
                                yearly_wflow_sources=extract_inDict_sum_variable_sources(
                                yearly_wflow_sources, wsources, '', field_water_period, leaching=False)
                                for chem in chemicals:
                                    field_chemicals_period = field_chemicals[chem][starting_periods[depth][i]:ending_periods[depth][i]+1]
                                    # activate code below to find wet spraying years
#                                    spraying_days_index=field_chemicals_period['Spray'][field_chemicals_period['Spray'] > 0].index
#                                    spraying_days=field_chemicals_period['Date'].loc[spraying_days_index]
#                                for d in spraying_days_index:
#                                    if field_water_period['Surface water'][d] > 0.5:
#                                        wet_spraying_dates.append(field_water_period['Date'][d].year)
#                                        all_spraying_dates_check.append(spraying_days)
                                    yearly_leaching_sources=extract_inDict_sum_variable_sources(
                                            yearly_leaching_sources, sources, chem, field_chemicals_period)
                                
                global_scores_sources=stats_on_source_list(yearly_leaching_sources, 
                                                           global_scores_sources, sources, chemicals)
                global_wflow_sources=stats_on_sum_list(yearly_wflow_sources, 
                                                       global_wflow_sources, wsources)   
                global_loads_sources=load_source_list(yearly_leaching_sources, 
                                                      global_loads_sources, sources, chemicals)   
            # save as csv files                                       
            pd.DataFrame(data=global_scores_sources, index=['OLOH','ILD00',
            'ILD08','ILD16','IHD00','IHD08']).to_csv(path_folder+'\\'+
            folder_name+'\\global_scores_sources_%s_%s_%s.csv'%(col, threshold, space))  
            pd.DataFrame(data=global_wflow_sources, index=['OLOH','ILD00',
            'ILD08','ILD16','IHD00','IHD08']).to_csv(path_folder+'\\'+
            folder_name+'\\global_wflow_sources_%s_%s_%s.csv'%(col, threshold, space))  
            pd.DataFrame(data=global_loads_sources, index=['OLOH','ILD00',
            'ILD08','ILD16','IHD00','IHD08']).to_csv(path_folder+'\\'+
            folder_name+'\\global_loads_sources_%s_%s_%s.csv'%(col, threshold, space))  