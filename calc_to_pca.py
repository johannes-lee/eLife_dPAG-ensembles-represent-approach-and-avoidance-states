## this script extracts neural data from each folder and use PCA to convert it to its PC components. interpolate_tracking.py must be run first.

import os
import scipy.io as sio
import numpy as np
import sklearn.decomposition as skd

from helpers import getfoldnames

foldnames = getfoldnames('../')

for fold in foldnames:
    if 'Calcium_pca2.mat' in os.listdir(fold):
        print('already pca\'ed: ' + fold)
        pass
    else:
        if 'track_interp2.mat' in os.listdir(fold):
            if not 'neural_data.mat' in os.listdir(fold):
                print('missing neural_data.mat: ' + fold)
                continue
            
            track = sio.loadmat(fold + '/track_interp2.mat')
            neur = sio.loadmat(fold + '/neural_data.mat')
            length = track['mouse_positionMS'].shape[0]
            c = {}
            c['C_raw'] = neur['neural_data']['C_raw'][0][0].T
            c['C'] = neur['neural_data']['C'][0][0].T
            lengthC = c['C_raw'].shape[0]
            
            c['C_raw'] = c['C_raw'][3:3+length]
            c['C'] = c['C'][3:3+length]
            
            pca = skd.PCA()
            c['C_pca'] = pca.fit_transform(c['C_raw'])
            c['pca_components'] = pca.components_
            c['pca_ev'] = pca.explained_variance_
            c['pca_evr'] = pca.explained_variance_ratio_
            
            #uncomment this next line when running
            sio.savemat(fold + '/Calcium_pca2.mat', c)
            print('just pca\'ed: ' + fold)
        else:
            print('~' + fold)
            pass