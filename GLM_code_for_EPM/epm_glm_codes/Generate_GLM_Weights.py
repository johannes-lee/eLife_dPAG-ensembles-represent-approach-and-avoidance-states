from helpers import *
from load_data import *
from GLM_Functions import *
import scipy.io as sio
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import hdf5storage
from functools import reduce
import statsmodels.api as sm

## Generate GLM weights for Rat and Shock assays, for both the beh-GLM with four behaviors only, and full GLM including all of the variables. 

## The keys in the data file are the mouse IDs. The first axis of the values is the cell ID, and the second axis is the variable. 

## The order of the variables is "Approach, Escape, Freeze, Stretch" in "Weights_Rat1_4vars" and "Weights_Shock_4vars", and it is "Distance to Rat, Mouse Speed, Rat Speed, Angle, Approach, Escape, Freeze, Stretch" in  "Weights_Rat1_allvars" and "Distance to Grid, Mouse Speed, Angle, Approach, Escape, Freeze, Stretch" in  "Weights_Shock_allvars"

## The dimensions of the data should be (N is the number of cells): 
## N*4 in "Weights_Rat1_4vars"
## N*4 in "Weights_Shock_4vars"
## N*8 in "Weights_Rat1_allvars"
## N*7 in "Weights_Shock_allvars"
                                            

########################################################

duration = 4
num_base_1 = 1
num_base_2 = 1 # nBases2 need to be odd 
bases_1 = make_raised_cosine_bases(duration, num_base_1)
bases_2 = make_raised_cosine_bases(duration, num_base_2)

#########################################################

data = {}
assay = 'Rat1'
assay_1 = 'rat1'
mouse = ['230', '355', '358', '362', '673', '674', '816', '825']
mouse_id = np.array([0, 1, 3, 4, 5, 6, 7])
foldnames, assay_idx, num_assay = load_assay(assay)

for i in mouse_id:
    c_raw, num_neuron, num_time, num_con, num_beh, num_var, X = load_data_distance_compare(foldnames[assay_idx[i]], assay, False)
    calc1, calc2, Cells1, Cells2 = getcalcs(foldnames, str(mouse[i]), 'epm', assay_1)
    imp = get_glm_weights_4vars(X, calc2, num_con, num_beh, num_base_2, bases_2)
    data[str(mouse[i])] = imp
    
sio.savemat('Weights_Rat1_4vars', data)

##########################################################

data = {}
assay = 'Shock'
assay_1 = 'shock'
mouse = ['230', '355', '358', '362', '673', '674', '816', '825']
mouse_id = np.array([0, 1, 2, 3, 4, 5, 6, 7])
foldnames, assay_idx, num_assay = load_assay(assay)

for i in mouse_id:
    c_raw, num_neuron, num_time, num_con, num_beh, num_var, X = load_data_distance_compare(foldnames[assay_idx[i]], assay, False)
    calc1, calc2, Cells1, Cells2 = getcalcs(foldnames, str(mouse[i]), 'epm', assay_1)
    imp = get_glm_weights_4vars(X, calc2, num_con, num_beh, num_base_2, bases_2)
    data[str(mouse[i])] = imp
    
sio.savemat('Weights_Shock_4vars', data)

##########################################################

data = {}
assay = 'Rat1'
assay_1 = 'rat1'
mouse = ['230', '355', '358', '362', '673', '674', '816', '825']
mouse_id = np.array([0, 1, 3, 4, 5, 6, 7])
foldnames, assay_idx, num_assay = load_assay(assay)

for i in mouse_id:
    c_raw, num_neuron, num_time, num_con, num_beh, num_var, X = load_data_distance_compare(foldnames[assay_idx[i]], assay, False)
    calc1, calc2, Cells1, Cells2 = getcalcs(foldnames, str(mouse[i]), 'epm', assay_1)
    imp = get_glm_weights_allvars(X, calc2, num_con, num_beh, num_base_1, num_base_2, bases_1, bases_2, num_var)
    data[str(mouse[i])] = imp
    
sio.savemat('Weights_Rat1_allvars', data)

##########################################################

data = {}
assay = 'Shock'
assay_1 = 'shock'
mouse = ['230', '355', '358', '362', '673', '674', '816', '825']
mouse_id = np.array([0, 1, 2, 3, 4, 5, 6, 7])
foldnames, assay_idx, num_assay = load_assay(assay)

for i in mouse_id:
    c_raw, num_neuron, num_time, num_con, num_beh, num_var, X = load_data_distance_compare(foldnames[assay_idx[i]], assay, False)
    calc1, calc2, Cells1, Cells2 = getcalcs(foldnames, str(mouse[i]), 'epm', assay_1)
    imp = get_glm_weights_allvars(X, calc2, num_con, num_beh, num_base_1, num_base_2, bases_1, bases_2, num_var)
    data[str(mouse[i])] = imp
    
sio.savemat('Weights_Shock_allvars', data)