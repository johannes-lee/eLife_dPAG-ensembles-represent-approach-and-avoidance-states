## adjusts indices to remove the first 3 frames (from redoing interpolation) in 'BehaviorMS.mat' and saves it in 'behav2.mat'

import os
import scipy.io as sio
import numpy as np
import sklearn.decomposition as skd

from helpers import getfoldnames

foldnames = getfoldnames('../')

for fold in foldnames:
    if 'behav2.mat' in os.listdir(fold):
        print('already adjusted behav: ' + fold)
    else:
        try:
            behav = sio.loadmat(fold + '/BehaviorMS.mat')
            for item in list(behav):
                try:
                    if 'Frame' in item:
                        behav[item] -= 3
                        while np.min(behav[item][0]) < 0:
                            behav[item] = behav[item][1:]
                    elif 'Indices' in item:
                        behav[item] = behav[item][3:-3]
                except:
                    pass
            sio.savemat(fold + '/behav2.mat', behav)
            print('adjusted behav for: ' + fold)
        except:
            print('missing BehaviorMS: ' + fold)