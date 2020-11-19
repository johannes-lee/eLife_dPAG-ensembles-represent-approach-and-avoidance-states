# This script reformats 'Tracking.mat' from a struct contatining all the variables to a .mat file 'track.mat' containing all the variables. This should be run before all other scripts.

import os
import scipy.io as sio
import numpy as np

from helpers import getfoldnames

foldnames = getfoldnames('../')

for fold in foldnames:
    if 'track_interp2.mat' in os.listdir(fold):
        print('already reformatted: ' + fold)
    else:
        try:
            t = sio.loadmat(fold + '/Tracking.mat')['Tracking']
            t2 = {}
            for k in range(len(t.dtype.names)):
                t2[t.dtype.names[k]] = t[t.dtype.names[k]][0][0]
            sio.savemat(fold + '/track.mat', t2)
            print('reformatted: ' + fold)
        except:
            print('could not open: ' + fold)

#def reformat():
#    try:
#        t = sio.loadmat('Tracking.mat')['Tracking']
#        try:
#            sio.loadmat('track.mat')
#            print('already reformatted: ' + os.getcwd())
#        except:
#            t2 = {}
#            for k in range(len(t.dtype.names)):
#                t2[t.dtype.names[k]] = t[t.dtype.names[k]][0][0]
#            #uncomment next line when running
#            sio.savemat('track.mat', t2)
#            print('reformatted: ' + os.getcwd())
#    except:
#        for name in os.listdir():
#            try:
#                os.chdir(name)
#                reformat()
#                os.chdir('../')
#            except:
#                pass
#    return
#
#def reformat():
#    
#
#reformat()