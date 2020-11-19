from helpers import *
from GLM_helpers import *
from load_data import *
import scipy.io as sio
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import hdf5storage
from scipy.stats import pearsonr, spearmanr, sem, t, ttest_ind, entropy, wilcoxon, ks_2samp, ttest_1samp, ranksums, zscore
from functools import reduce
from sklearn.linear_model import LinearRegression, Lasso
import statsmodels.api as sm
from sklearn.metrics import mean_squared_error, mutual_info_score, r2_score, explained_variance_score, auc, confusion_matrix
from sklearn.utils import resample
from sklearn.feature_selection import f_regression, mutual_info_regression
from sklearn.preprocessing import StandardScaler
import seaborn as sns
from scipy import integrate
from numpy.polynomial.polynomial import polyfit
from sklearn import mixture
from sklearn.decomposition import PCA
import random




# Get indices of a type of 'assay' across all folds in 'foldnames'
def get_assay_fold_index(foldnames, assay):
    assay_idx = []
    for i in range(len(foldnames)):
        if os.path.basename(foldnames[i]) == assay:
            assay_idx.append(i)       
    return assay_idx

# Get all the folds, indices of a type of 'assay' in those folds, and the number of such assays
def load_assay(assay):
    foldnames = getfoldnames('../')
    #foldnames = getfoldnames('/Users/liujh/Project/data/dPAG_2019_03/')
    assay_idx = get_assay_fold_index(foldnames, assay)
    #assay_fold = foldnames[assay_idx[0]]
    #toyrat_idx = get_assay_fold_index(foldnames, 'ToyRat')
    #toyrat_fold = foldnames[toyrat_idx[1]]
    num_assay = len(assay_idx)
    return foldnames, assay_idx, num_assay


# Load behavioral and neural data
def load_data_distance_compare(fold, assay, compare):
    if assay == 'Rat1' or assay == 'Rat2' or assay == 'Rat3':
        trk = sio.loadmat(fold + '/track_interp2.mat')
        beh = sio.loadmat(fold + '/behav2.mat')
        neu = sio.loadmat(fold + '/Calcium_pca2.mat')
        app = beh['approachIndicesMS']
        esc = beh['escapeIndicesMS']
        fre = beh['freezeIndicesMS']
        stc = beh['stretchIndicesMS']
        c_raw = neu['C_raw']
        dist = trk['distanceMouseRatMS']
        mvel = trk['mouseVelMS']
        rvel = trk['ratVelMS']
        adiff = trk['angleDiffMouseHeadDirRatMS']
        num_neuron = c_raw.shape[1]
        num_time = c_raw.shape[0]
        if compare == True:
            num_con = 3
            num_beh = 4
            num_var = num_con + num_beh
            X = np.zeros((num_time, num_var))
            dist = dist.reshape(-1)
            X[:, 0] = dist[0:num_time]
            X[:, 1] = mvel[0:num_time].reshape(-1)
            X[:, 2] = adiff[0:num_time].reshape(-1)
            X[:, 3] = app[0:num_time].reshape(-1)
            X[:, 4] = esc[0:num_time].reshape(-1)
            X[:, 5] = fre[0:num_time].reshape(-1)
            X[:, 6] = stc[0:num_time].reshape(-1)
        else:
            num_con = 4
            num_beh = 4
            num_var = num_con + num_beh
            X = np.zeros((num_time, num_var))
            dist = dist.reshape(-1)
            X[:, 0] = dist[0:num_time]
            X[:, 1] = mvel[0:num_time].reshape(-1)
            X[:, 2] = rvel[0:num_time].reshape(-1)
            X[:, 3] = adiff[0:num_time].reshape(-1)
            X[:, 4] = app[0:num_time].reshape(-1)
            X[:, 5] = esc[0:num_time].reshape(-1)
            X[:, 6] = fre[0:num_time].reshape(-1)
            X[:, 7] = stc[0:num_time].reshape(-1)
        return c_raw, num_neuron, num_time, num_con, num_beh, num_var, X
    elif assay == 'ToyRat':
        trk = sio.loadmat(fold + '/track_interp2.mat')
        beh = sio.loadmat(fold + '/behav2.mat')
        neu = sio.loadmat(fold + '/Calcium_pca2.mat')
        app = beh['approachIndicesMS']
        esc = beh['escapeIndicesMS']
        fre = beh['freezeIndicesMS']
        stc = beh['stretchIndicesMS']
        c_raw = neu['C_raw']
        dist = trk['distanceMouseToyRatMS']
        mvel = trk['mouseVelMS']
        adiff = trk['angleDiffMouseHeadDirRatMS']
        num_neuron = c_raw.shape[1]
        num_time = c_raw.shape[0]
        num_con = 3
        num_beh = 4
        num_var = num_con + num_beh
        X = np.zeros((num_time, num_var))
        
        dist = dist.reshape(-1)
        for i in range(num_time):
            if np.isnan(dist[i]) == True:
                dist[i] = dist[i-1]
        X[:, 0] = dist[0:num_time]
        X[:, 1] = mvel[0:num_time].reshape(-1)
        X[0:len(adiff), 2] = adiff[0:num_time].reshape(-1)
        X[:, 3] = app[0:num_time].reshape(-1)
        X[:, 4] = esc[0:num_time].reshape(-1)
        X[:, 5] = fre[0:num_time].reshape(-1)
        X[:, 6] = stc[0:num_time].reshape(-1)

        return c_raw, num_neuron, num_time, num_con, num_beh, num_var, X
    
    elif assay == 'Shock' or assay == 'ShockExt' or assay == 'ShockExt2' or assay == 'ShockExt3' or assay == 'ShockHab':
        trk = sio.loadmat(fold + '/track_interp2.mat')
        beh = sio.loadmat(fold + '/behav2.mat')
        neu = sio.loadmat(fold + '/Calcium_pca2.mat')
        app = beh['approachIndicesMS']
        esc = beh['escapeIndicesMS']
        fre = beh['freezeIndicesMS']
        stc = beh['stretchIndicesMS']
        c_raw = neu['C_raw']
        dist = trk['distanceShockgridMS']
        mvel = trk['mouseVelMS']
        adiff = np.pi - np.absolute(trk['mouseAngleMS'])
        num_neuron = c_raw.shape[1]
        num_time = c_raw.shape[0]
        num_con = 3
        num_beh = 4
        num_var = num_con + num_beh
        X = np.zeros((num_time, num_var))
        
        dist = dist.reshape(-1)
        for i in range(num_time):
            if np.isnan(dist[i]) == True:
                dist[i] = dist[i-1]
        for i in range(num_time):
            if np.isnan(adiff[i]) == True:
                if i == 0:
                    adiff[i] = 0
                else:
                    adiff[i] = adiff[i-1]
        X[:, 0] = dist[0:num_time]
        X[:, 1] = mvel[0:num_time].reshape(-1)
        X[:, 2] = adiff[0:num_time].reshape(-1)
        X[:, 3] = app[0:num_time].reshape(-1)
        X[:, 4] = esc[0:num_time].reshape(-1)
        X[:, 5] = fre[0:num_time].reshape(-1)
        X[:, 6] = stc[0:num_time].reshape(-1)

        return c_raw, num_neuron, num_time, num_con, num_beh, num_var, X
    
    
    
# Get raised cosine bases with number 'nBases'   
def make_raised_cosine_bases(duration, nBases):
    total_duration = int(duration * np.power(1.5, nBases-1))
    Bases = np.zeros((total_duration, nBases))
    Bases[0, 0] = 1 
    for i in range(nBases-1):
        period = duration * 2
        t = int(duration)
        x = np.arange(t)
        cosx = np.sin(x * 2 * np.pi / period)
        #Bases[int(t/2):int(t/2)+len(cosx), i+1] = cosx
        Bases[0:len(cosx), i+1] = cosx  
        duration = duration * 1.5
    return Bases


def glm_weights_cell_beh(X, c_raw, num_beh, num_base_2, bases_2, neuron_idx, beh_idx):
    
    X_ori = X
    num_time = X.shape[0]
    num_var = num_beh
    num_base_all = num_beh * num_base_2
    X = StandardScaler().fit_transform(X)
    X_1 = np.zeros((num_time, num_base_all))
    for i in range(num_beh):
        for j in range(num_base_2):
            X_1[:, i*num_base_2+j] = np.convolve(X[:, i], bases_2[:, j])[0:num_time] 
    X_1 = sm.add_constant(X_1, prepend=False)
    
    
    beh_bool = X_ori[:, beh_idx] == 1
    X_new = X_1[beh_bool, :]
    y = c_raw[beh_bool, neuron_idx]
    glm_model = sm.GLM(y, X_new, family=sm.families.Gaussian())
    glm_results = glm_model.fit()
    w = glm_results.params
    
    return w, glm_results, X_new


def glm_weights_cell_4vars(X, c_raw, num_beh, num_base_2, bases_2, neuron_idx):
    
    num_time = X.shape[0]
    num_var = num_beh
    num_base_all = num_beh * num_base_2
    w_all = np.zeros((num_base_all+1, 4))
    for i in range(num_beh):
        if 1 not in X[:, i]:
            continue
        w, glm_results, X_new = glm_weights_cell_beh(X, c_raw, num_beh, num_base_2, bases_2, neuron_idx, i)
        w_all[:, i] = w
    
    return w_all


def get_glm_weights_4vars(X, c_raw, num_con, num_beh, num_base_2, bases_2):
    imp_2 = np.zeros((c_raw.shape[1], num_beh))
    X_beh = X[:, num_con:]
    c_new = np.zeros(c_raw.shape)
    for i in range(c_raw.shape[1]):
        c_new[:, i] = zscore(c_raw[:, i])
    for i in range(c_raw.shape[1]):
        #w_all, y_cov_app, y_cov_esc, y_cov_fre, y_cov_stc = glm_training_onset(X_beh, c_raw, num_beh, num_base_2, bases_2, i)
        w_all = glm_weights_cell_4vars(X_beh, c_new, num_beh, num_base_2, bases_2, i)
        for j in range(num_beh):
            beh_idx = X_beh[:, j]
            if np.nonzero(beh_idx)[0].size:
                imp_2[i, j] = w_all[j*num_base_2, j]
            else:
                imp_2[i, j] = 0

            
    return imp_2

def glm_weights_cell_allvars(X, c_raw, num_con, num_beh, num_base_1, num_base_2, bases_1, bases_2, neuron_idx):
    
    X_ori = X
    num_time = X.shape[0]
    num_var = num_beh + num_con
    num_base_all = num_con * num_base_1 + num_beh * num_base_2
    w_all = np.zeros((num_base_all+1, 8))
    X = StandardScaler().fit_transform(X)
    X_1 = np.zeros((num_time, num_base_all))
    for i in range(num_var):
        for j in range(num_base_2):
            X_1[:, i*num_base_2+j] = np.convolve(X[:, i], bases_2[:, j])[0:num_time] 
    X_1 = sm.add_constant(X_1, prepend=False)
    X_new = X_1
    y = c_raw[:, neuron_idx]
    glm_model = sm.GLM(y, X_new, family=sm.families.Gaussian())
    glm_results = glm_model.fit()
    w = glm_results.params
    #w_all[:, i] = w
    return w

def get_glm_weights_allvars(X, c_raw, num_con, num_beh, num_base_1, num_base_2, bases_1, bases_2, num_var):
    imp_2 = np.zeros((c_raw.shape[1], num_con+num_beh))
    c_new = np.zeros(c_raw.shape)
    for i in range(c_raw.shape[1]):
        c_new[:, i] = zscore(c_raw[:, i])
    for i in range(c_raw.shape[1]):
        w = glm_weights_cell_allvars(X, c_new, num_con, num_beh, num_base_1, num_base_2, bases_1, bases_2, i)
        imp_2[i, :] = w[0:num_con+num_beh]


            
    return imp_2