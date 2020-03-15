# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 17:32:49 2020

@author: cbn978
"""
from collections import defaultdict
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from datetime import datetime
import sys
sys.path.append(r'C:\Users\cbn978\Documents\PESTCAST\Git\PyDaisy\pydaisy')

from pydaisy.Daisy import DaisyDlf

# =============================================================================
# file extraction
# =============================================================================
def extract(my_path_folder,file, my_start_month, my_start_day, my_end_month, my_end_day):
    my_extraction=DaisyDlf(my_path_folder+'\\'+'%s.dlf' %file)
    my_extraction=my_extraction.Data
    years=np.unique([my_extraction.index[i].year for i in range(len(my_extraction.index))])
    my_extraction['Date']=my_extraction.index
    my_extraction.index=np.arange(0,len(my_extraction))
    my_start, my_end = [], []
    for i in range(len(years)-1):
        my_start.append(my_extraction['Date'][my_extraction['Date'] == datetime.strptime('%0d-%02d-%02d 00:00:00'%(years[i+1], my_start_month, my_start_day), "%Y-%m-%d %H:%M:%S")].index[0])
        my_end.append(my_extraction['Date'][my_extraction['Date'] == datetime.strptime('%0d-%02d-%02d 00:00:00'%(years[i+1], my_end_month, my_end_day), "%Y-%m-%d %H:%M:%S")].index[0])            
    return my_extraction, my_start, my_end
            
      
def stats_on_sum_list(my_dict_inputs, my_dict_outputs, my_lists):
    for my_list in my_lists:     
        if len(my_dict_inputs[my_list]) > 0:
            my_dict_outputs[my_list].append([
                                   np.median(my_dict_inputs[my_list]),
                                   np.percentile(my_dict_inputs[my_list],10),
                                   np.percentile(my_dict_inputs[my_list],90)])
    return my_dict_outputs

def stats_on_source_list(my_dict_inputs, my_dict_outputs, my_sources, my_chemicals):
    for my_chemical in my_chemicals:
        for my_source in my_sources:   
            if len(my_dict_inputs[my_source+'_'+my_chemical]) > 0:      
                my_dict_outputs[my_source+'_'+my_chemical].append([
                                       np.median(my_dict_inputs[my_source+'_'+my_chemical]),
                                       np.percentile(my_dict_inputs[my_source+'_'+my_chemical],10),
                                       np.percentile(my_dict_inputs[my_source+'_'+my_chemical],90)])
    return my_dict_outputs

def load_source_list(my_dict_inputs, my_dict_outputs, my_sources, my_chemicals):
    for my_chemical in my_chemicals:
        for my_source in my_sources:   
            if len(my_dict_inputs[my_source+'_'+my_chemical]) > 0:                                          
                my_dict_outputs[my_source+'_'+my_chemical].append(my_dict_inputs[my_source+'_'+my_chemical])

    return my_dict_outputs
    
def extract_inDict_sum_variable_sources(my_dict, my_variables, my_dict_drawer, 
                                        my_extracted_file, leaching=True):
    if leaching:
            for my_variable in my_variables:
                my_dict[my_variable+'_'+my_dict_drawer].append(np.sum(my_extracted_file[my_variable]))
    else:
        for my_variable in my_variables:
            my_dict[my_variable].append(np.sum(my_extracted_file[my_variable]))
    return my_dict


# =============================================================================
# analysis
# =============================================================================
def weigh_field_areas_sources_stats(perc_field, perc_WT, perc_drainco, perc_overlap, perc_awd, perc_ad, my_field_file, my_drainco_file): 
    dict_areas_med, dict_areas_10p, dict_areas_90p = defaultdict(), defaultdict(), defaultdict()
    for stati,stat in enumerate([dict_areas_med, dict_areas_10p, dict_areas_90p]):        
        stat['OLOH_only']=eval(my_field_file.values[0])[stati]*perc_awd + eval(my_drainco_file.values[0])[stati]*perc_ad 
        stat['OLOH']=eval(my_field_file.values[0])[stati]*perc_field + eval(my_drainco_file.values[0])[stati]*perc_drainco
        stat['ILD00']=eval(my_field_file.values[1])[stati]*perc_WT + eval(my_drainco_file.values[1])[stati]*perc_overlap
        stat['ILD08']=eval(my_field_file.values[2])[stati]*perc_WT + eval(my_drainco_file.values[2])[stati]*perc_overlap
        stat['ILD16']=eval(my_field_file.values[3])[stati]*perc_WT + eval(my_drainco_file.values[3])[stati]*perc_overlap
        stat['IHD00']=eval(my_field_file.values[4])[stati]*perc_WT + eval(my_drainco_file.values[4])[stati]*perc_overlap
        stat['IHD08']=eval(my_field_file.values[5])[stati]*perc_WT + eval(my_drainco_file.values[5])[stati]*perc_overlap    
    return  dict_areas_med, dict_areas_10p, dict_areas_90p


def weigh_field_areas_total_loads(perc_field, perc_WT, perc_drainco, perc_overlap, 
                              perc_awd, perc_ad, my_field_file, my_drainco_file): 
    dict_ally_loads = defaultdict()
    
    nw_only_loads=(np.array(eval(my_field_file.values[0]))*perc_awd + 
                    np.array(eval(my_drainco_file.values[0]))*perc_ad)
    dict_ally_loads['OLOH_only']=nw_only_loads
    
    nw_loads=(np.array(eval(my_field_file.values[0]))*perc_field + 
             np.array(eval(my_drainco_file.values[0]))*perc_drainco)
        
    wt0_loads=(np.array(nw_loads)+np.array(eval(my_field_file.values[1]))*perc_WT + 
              np.array(eval(my_drainco_file.values[1]))*perc_overlap)
    dict_ally_loads['ILD00']=wt0_loads
    
    wt8_loads=(np.array(nw_loads)+np.array(eval(my_field_file.values[2]))*perc_WT + 
                    np.array(eval(my_drainco_file.values[2]))*perc_overlap)
    dict_ally_loads['ILD08']=wt8_loads
    
    wt16_loads=(np.array(nw_loads)+np.array(eval(my_field_file.values[3]))*perc_WT + 
                    np.array(eval(my_drainco_file.values[3]))*perc_overlap)
    dict_ally_loads['ILD16']=wt16_loads
    
    wt0hm_loads=(np.array(nw_loads)+np.array(eval(my_field_file.values[4]))*perc_WT + 
                    np.array(eval(my_drainco_file.values[4]))*perc_overlap)
    dict_ally_loads['IHD00']=wt0hm_loads
    
    wt8hm_loads=(np.array(nw_loads)+np.array(eval(my_field_file.values[5]))*perc_WT + 
                    np.array(eval(my_drainco_file.values[5]))*perc_overlap)
    dict_ally_loads['IHD08']=wt8hm_loads
     
    return dict_ally_loads
