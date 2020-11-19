# Uses linear interpolation to interp the variables in interpvars, saves the data in 'track_interp.mat'

import os
from scipy.interpolate import interp1d
import numpy as np
import scipy.io as sio

from helpers import getfoldnames

foldnames = getfoldnames('../')

interpvars = ['mouse_positionMS', 'mouseVelMS', 'mouseAngleMS', 'rat_positionMS', 'ratAngleMS', 'angleDiffMouseHeadDirRatMS',
             'ratVelMS', 'distanceMouseRatMS', 'mouseVelScaledMS', 'ratVelScaledMS', 'distanceMouseRatScaledMS', 'distanceShockgridMS',
             'distanceShockgridScaledMS', 'angleDiffShockgridMS', 'distanceMouseToyRatMS']
assays = ['EPM', 'OpenField', 'Rat_1', 'Rat_2', 'Rat_3', 'Shock', 'Shock_ext1']

for fold in foldnames:
    if 'track_interp2.mat' in os.listdir(fold):
        print('already interpolated: ' + fold)
        pass
    else:
        track = sio.loadmat(fold + '/track.mat')
        for var in interpvars:
            try:
                track[var] = track[var][3:-3]
                length = track[var].shape[0]
                valind = np.arange(length)[np.all(~np.isnan(track[var]), axis = 1)]
                interpolator = interp1d(valind, track[var][valind], axis = 0)
                track[var] = interpolator(np.arange(np.max(valind)))
            except:
                pass
        
        if 'mouse_positionMS' in track:
            
            #if '673' in fold:
            #    if 'epm' in fold.lower():
            #        mpos = np.array(track['mouse_positionMS'])
            #        track['mouse_positionMS'][:, 0] = mpos[:, 1]
            #        track['mouse_positionMS'][:, 1] = mpos[:, 0]
            
            #uncomment next line when running
            sio.savemat(fold + '/track_interp2.mat', track)
            print('interpolated: ' + fold)
        else:
            print('not interpolated: ' + fold)