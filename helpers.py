import os
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from matplotlib.collections import LineCollection
from sklearn.linear_model import LinearRegression as LR
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import ranksums, ttest_1samp
from scipy.stats import pearsonr
from scipy.interpolate import interp2d
from sklearn.decomposition import PCA
from scipy.stats import sem, zscore
import pdb
from inds_outliers import *

#set font size for pyplot
plt.rcParams.update({'font.size': 16})

# Use getfoldnames to return all sub-directories containing a file 'track.mat'
# Input should be empty, or a string and/or list
# If string input: returns all sub-dirs inside the given directory string input
# If list input: returns sub-dirs inside current dir appended to the input list
# If string and list input: returns sub-dirs inside input string dir appended to input list
# If empty input: returns sub-dirs inside current dir
def getfoldnames(*args):
    if len(args) > 2:
        raise Exception('Too many inputs')
    elif len(args) == 1:
        if type(args[0]) == str:
            if os.path.isdir(args[0]):
                startdir = os.getcwd()
                foldnames_out = getfoldnames2(args[0])
                foldnames_out.sort()
                os.chdir(startdir)
                return foldnames_out
            else:
                raise Exception('The input string is not a valid directory')
        elif type(args[0]) == list:
            foldnames_out = getfoldnames2(args[0])
            foldnames_out.sort()
            return foldnames_out
        else:
            raise Exception('The format of the input is incorrect. Input should be empty, a string, a list, or a string and a list')
    elif len(args) == 2:
        if type(args[0]) == str and type(args[1]) == list:
            if os.path.isdir(args[0]):
                startdir = os.getcwd()
                os.chdir(args[0])
                foldnames_out = getfoldnames2(args[1])
                foldnames_out.sort()
                os.chdir(startdir)
                return foldnames_out
            else:
                raise Exception('The input string is not a valid directory')
        else:
            raise Exception('The format of the input is incorrect. Input should be empty, a string, a list, or a string and a list')
    else:
        foldnames_out = getfoldnames2()
        foldnames_out.sort()
        return foldnames_out

# getfoldnames2 is an internal recursive function used by getfoldnames
# basically the same as getfoldnames but doesn't support multiple inputs AND doesn't stay in the same directory when run
def getfoldnames2(*args):
    if len(args) > 1:
        raise Exception('too many inputs')
    elif len(args) > 0:
        if type(args[0]) == str:
            os.chdir(args[0])
            folds = []
        elif type(args[0]) == list:
            folds = args[0]
        else:
            raise Exception('The format of the input is incorrect. Input should be empty, a string, a list, or a string and a list')
    else:
        folds = []
    if 'Tracking.mat' in os.listdir():
        folds.append(os.getcwd())
    else:
        for name in os.listdir():
            if os.path.isdir(name):
                os.chdir(name)
                getfoldnames(folds)
                os.chdir('../')
    return folds

# getmousenums returns the numbers of the mice whose data is found in:
# sub-directories of the current directory (no input)
# sub-directories of the target directory (input directory string)
def getmousenums(*args):
    if len(args) == 0:
        mousenums_out = getmousenums2([])
        mousenums_out.sort()
        return mousenums_out
    elif len(args) == 1:
        if type(args[0]) == str:
            if os.path.isdir(args[0]):
                startdir = os.getcwd()
                os.chdir(args[0])
                mousenums_out = getmousenums2([])
                mousenums_out.sort()
                os.chdir(startdir)
                return mousenums_out
            else:
                raise Exception('The input string is not a valid directory')
        else:
            raise Exception('The format of the input is incorrect. Input should be empty or a string')
    raise Exception('The format of the input is incorrect. Input should be empty or a string')

# internal recursive function for getmousenums
def getmousenums2(folds):
    for name in os.listdir():
        try:
            int(name)
            folds.append(name)
        except:
            try:
                os.chdir(name)
                getmousenums(folds)
                os.chdir('../')
            except:
                pass
    return folds

def labelfoldbymouse(folds, mousenums):
    labels = ['']*len(folds)
    for i, fold in enumerate(folds):
        for mouse in mousenums:
            if mouse in fold:
                labels[i] = mouse
                break
    return labels

def labelfoldbyassay(folds):
    assays = ['EPM', 'OpenField', 'Rat_1', 'Rat_2', 'Rat_3', 'Shock', 'Shock_ext']
    assaysout = ['epm', 'openfield', 'rat1', 'rat2', 'rat3', 'shock', 'shockext']
    labels = ['']*len(folds)
    for i, fold in enumerate(folds):
        for k, assay in enumerate(assays):
            if assay in fold:
                labels[i] = assaysout[k]
                break
    return labels

# plotlinecolor plots a 2D line of xy_in from indices start to end (xy_in[start:end]) with color according the data[start:end] on the axis ax_in
# xy_in is assumed to be 2-dimensional (not checked)
# to plot the whole array, use start = 0, end = None
# if interpolate == True: interpolates xy_in and data between points based on interptype
#     when using interpolate, end-start must be 3 or greater
# Set the colormap using cmap
# Set alpha value with alpha
# if useraw == True: don't rescale data based on 10th and 90th percentiles.
def plotlinecolor(start, end, xy_in, data, ax_in, interpolate, **kwargs):
    xy = xy_in.astype('float64')
    cmap = 'jet'
    interptype = 'quadratic'
    alpha = 1.0
    percentiles = (0, 100)
    vlim = None
    lw = 2
    if ax_in == None:
        plt.figure(figsize=(15, 15))
        ax_in = plt.gca()
    for key, value in kwargs.items():
        if key == 'cmap':
            cmap = value
        if key == 'interptype':
            interptype = value
        if key == 'alpha':
            alpha = value
        if key == 'percentiles':
            if type(value) == tuple:
                if len(value) == 2:
                    percentiles = value
        if key == 'vlim':
            if type(value) == tuple:
                if len(value) == 2:
                    vlim = value
        if key == 'lw':
            lw = value
    if interpolate:
        #length = end - start must equal 3 or more
        xlim = (xy[start:end, 0].min(), xy[start:end, 0].max())
        ylim = (xy[start:end, 1].min(), xy[start:end, 1].max())
        if percentiles[0] != 0 and percentiles[1] != 100:
            data10_90 = (np.percentile(data[start:end], percentiles[0]), np.percentile(data[start:end], percentiles[1]))
            data_thresh = data[start:end]*(data[start:end] >= data10_90[0])*(data[start:end] <= data10_90[1]) + (
                data10_90[0])*(data[start:end] < data10_90[0]) + (data10_90[1])*(data[start:end] > data10_90[1])
        else:
            data_thresh = data[start:end]
        
        xy_dif = np.diff(xy[start:end], axis = 0)
        xrange = (xlim[1] - xlim[0])
        yrange = (ylim[1] - ylim[0])
        
        if xrange > 0:
            xy_dif[:, 0] = xy_dif[:, 0] / xrange
        else:
            xy_dif[:, 0] = 0
        if yrange > 0:
            xy_dif[:, 1] = xy_dif[:, 0] / yrange
        else:
            xy_dif[:, 1] = 0
        dist_dif = np.linalg.norm(xy_dif, axis = 1)
        data_thresh_dif = np.diff(data_thresh)
        
        length = xy[start:end].shape[0]
        ind = np.arange(length)
        
        dist_dif_max = dist_dif.max()
        dind = 1
        if dist_dif_max > 0:
            dind = 0.01 / dist_dif_max
        dind = np.minimum(dind, 1)
        ind_out = np.arange(0, length - 1, dind)
        interpolator = interp1d(ind, xy[start:end], axis = 0, kind = interptype)
        xy_interp = interpolator(ind_out)
        interpolator = interp1d(ind, data[start:end], axis = 0, kind = interptype)
        data_plot = interpolator(ind_out)
        
        points = xy_interp.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        if percentiles[0] != 0 and percentiles[1] != 100:
            norm = plt.Normalize(data10_90[0], data10_90[1])
        elif vlim != None:
            norm = plt.Normalize(*vlim)
        else:
            norm = plt.Normalize(data[start:end].min(), data[start:end].max())
        lc = LineCollection(segments, cmap = cmap, norm = norm, alpha = alpha)
        lc.set_array(data_plot)
        lc.set_linewidth(lw)
        ax_in.add_collection(lc)
        plt.plot()
        return norm
    else:
        points = xy[start:end].reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        data_plot = data[start:end]
        if percentiles[0] != 0 and percentiles[1] != 100:
            data10_90 = (np.percentile(data[start:end], percentiles[0]), np.percentile(data[start:end], percentiles[1]))
            norm = plt.Normalize(data10_90[0], data10_90[1])
        elif vlim != None:
            norm = plt.Normalize(*vlim)
        else:
            norm = plt.Normalize(data[start:end].min(), data[start:end].max())
        lc = LineCollection(segments, cmap = cmap, norm = norm, alpha = alpha)
        lc.set_array(data_plot)
        lc.set_linewidth(lw)
        ax_in.add_collection(lc)
        dind = 1
        plt.plot()
        return norm
    return

#like plotlinecolor but with no 10th and 90th percentile rescaling, and no interpolation, and always plotting jet
def plotlinecolor1(start, end, pos, data, ax, alpha):
    points = pos[start:end].reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    data_plot = data[start:end]
    norm = plt.Normalize(np.min(data_plot), np.max(data_plot))
    lc = LineCollection(segments, cmap = 'jet', norm = norm, alpha = alpha)
    lc.set_array(data_plot)
    ax.add_collection(lc)
    plt.plot()
    return

# returns the data parameters from input dictionaries track and behav
# data variables are:
# 0: approach
# 1: escape
# 2: freeze
# 3: stretch
# 4: mouse velocity
# 5: threat velocity (if applicable)
# 6: mouse-threat distance
# 7: mouse-threat angle
# 8: mouse X position
# 9: mouse Y position
def getdata(track, behav):
    behavFrame = ['approachFrameMS', 'escapeFrameMS', 'freezeFrameMS', 'stretchFrameMS']
    length = track['mouse_positionMS'].shape[0]
    data = np.zeros((length, 10))
    data[:, 8:10] = track['mouse_positionMS']
    if 'mouseVelMS' in track:
        data[:, 4] = track['mouseVelMS'][:, 0]
        varNames = ['--', '--', '--', '--', 'mouse velocity', '--', '--', 
            '--', 'mouse X position', 'mouse Y position']
    if 'ratVelMS' in track:
        data[:, 5] = track['ratVelMS'][:, 0]
        data[:, 6] = track['distanceMouseRatMS'][:, 0]
        data[:, 7] = track['angleDiffMouseHeadDirRatMS'][:, 0]
        varNames = ['approach', 'escape', 'freeze', 'stretch', 'mouse velocity', 'rat velocity', 
                    'mouse-rat distance', 'mouse-rat angle', 'mouse X position', 'mouse Y position']
    elif 'distanceShockgridMS' in track:
        data[:, 6] = track['distanceShockgridMS'][:, 0]
        data[:, 7] = track['mouseAngleMS'][0:length, 0]
        varNames = ['approach', 'escape', 'freeze', 'stretch', 'mouse velocity', '--', 
                         'mouse-shockgrid distance', 'mouse-shockgrid angle', 'mouse X position', 
                         'mouse Y position']
    else:
        varNames = ['--', '--', '--', '--', '--', '--', '--', 
            '--', 'mouse X position', 'mouse Y position']
    for i, frame in enumerate(behavFrame):
        if frame in behav:
            for k in range(behav[frame].shape[0]):
                data[behav[frame][k, 0]:behav[frame][k, 1], i] = 1
    return data, varNames

def plot_pcs(pos, datapca, pcs_to_plot, *args):
    numPCs = np.size(pcs_to_plot)
    width = int(np.sqrt(numPCs))
    height = int(np.ceil(numPCs / width))
    plt.figure(figsize=(7*width, 7*height))
    for k in range(numPCs):
        plt.subplot(height, width, k+1)
        ax = plt.gca()
        plotlinecolor(0, None, pos, datapca[:, k], ax, 1, True)
        plt.xlabel('X position')
        plt.ylabel('Y position')
        plt.title('PC ' + str(k + 1))
    plt.tight_layout()
    if len(args) > 0:
        savename = args[0]
    if savename.find('.png') >= 0:
        plt.savefig(savename)
    plt.close()
    return

def getshareddata(CR, inds, neurs):
    if type(neurs) == list:
        cells = getsharedcells(CR, inds)
        if type(neurs[0]) == dict:
            pass
        elif type(neurs[0]) == np.ndarray:
            data = []
            for i, ind in enumerate(inds):
                data.append(neurs[i][:, cells[i]])
            return data
    else:
        return False

def getsharedcells(CR, inds):
    valinds = np.ones(CR.shape[0]).astype('bool')
    inds = np.array(inds).tolist()
    if type(inds) == list:
        for ind in inds:
            valinds *= CR[:, ind] != 0
        cells = []
        for ind in inds:
            cells.append(CR[valinds, ind] - 1)
        return cells
    else:
        return False

def lrpredict(in1, out1, in2, **kwargs):
    returnclass = False
    for key, value in kwargs.items():
        if key == 'returnclass':
            returnclass = value
    lr = LR()
    lr.fit(in1, out1)
    out2 = lr.predict(in2)
    if returnclass == True:
        return out2, lr
    return out2

def scatter3d(x, *args, **kwargs):
    ax = None
    fig = None
    y = None
    z = None
    cmap = 'jet'
    c = None
    s = 1
    if len(args) > 0:
        y = args[0]
        if len(args) > 1:
            z = args[1]
    for key, value in kwargs.items():
        if key == 'ax':
            ax = value
        if key == 'fig':
            fig = value
        if key == 'cmap':
            cmap = value
        if key == 'c':
            c = value
        if key == 's':
            s = value
    if fig == None:
        fig = plt.figure(figsize=(15, 15))
        ax = fig.add_subplot(111, projection = '3d')
        if y != None and z != None:
            ax.scatter(x, y, z, s = s, cmap = cmap, c = c)
        else:
            ax.scatter(x[:, 0], x[:, 1], x[:, 1], s = s, cmap = cmap, c = c)
    return

def plot3d(x, y, z, **kwargs):
    return

#return folder containing data for mousenum and assay.
def getfold(foldnames, mousenum, assay):
    if type(mousenum) == int:
        mousenum = str(mousenum)
    assaylength = len(assay)
    for i, fold in enumerate(foldnames):
        fold = fold.lower()
        if mousenum in fold:
            if assay == fold[-assaylength:] and fold[-(assaylength + 1)] == '\\':
                return foldnames[i]
                break
    return False

def getdicts(foldnames, mousenum, assay):
    mfold = getfold(foldnames, mousenum, assay)
    try:
        d1 = sio.loadmat(mfold + '/track_interp2.mat')
    except:
        d1 = {}
    try:
        d2 = sio.loadmat(mfold + '/Calcium_pca2.mat')
    except:
        d2 = {}
    try:
        d3 = sio.loadmat(mfold + '/behav2.mat')
    except:
        d3 = {}
    return [d1, d2, d3]


#plot the heatmap of data using the xy points; ax is the pyplot axis.
def plotheatmap(xy, data, **kwargs):
    figsize=(15, 15)
    cmap = 'jet'
    threshold = 0.00
    avg = True
    HH = None
    WW = None
    step = 20
    interpstep = 1
    colorbar = False
    offset = 2
    returnxy = False
    for key, item in kwargs.items():
        if key == 'ax':
            ax = item
        if key == 'avg':
            avg = item
        if key == 'HH':
            HH = item
        if key == 'WW':
            WW = item
        if key =='step':
            step = item
        if key == 'interpstep':
            interpstep = item
        if key == 'colorbar':
            colorbar = item
        if key == 'offset':
            offset = item
        if key == 'returnxy':
            returnxy = item
    if not 'ax' in kwargs:
        plt.figure(figsize=figsize)
        ax = plt.gca()
        
    
    data = data.flatten()
    data = data - np.min(data)
    xbin = np.arange(step*int(np.min(xy[:, 0]) / step) - offset*step, step*np.ceil(np.max(xy[:, 0]) / step) + offset*step, step)
    ybin = np.arange(step*int(np.min(xy[:, 1]) / step) - offset*step, step*np.ceil(np.max(xy[:, 1]) / step) + offset*step, step)
    avgs = -0.1*np.max(data)*np.ones((xbin.size, ybin.size))
    for i, x in enumerate(xbin):
        for j, y in enumerate(ybin):
            bindata = data[(x <= xy[:, 0])*(xy[:, 0] < x + step)*(y <= xy[:, 1])*(xy[:, 1] < y + step)]
            if bindata.size > 0:
                if avg == True:
                    avgs[i, j] = np.mean(bindata)
                else:
                    avgs[i, j] = np.sum(bindata)
    xinterp = np.arange(np.min(xy[:, 0]) - offset*step, np.max(xy[:, 0]) + offset*step, interpstep)
    yinterp = np.arange(np.min(xy[:, 1]) - offset*step, np.max(xy[:, 1]) + offset*step, interpstep)
    
    interpolater = interp2d(ybin+step/2, xbin+step/2, avgs, kind = 'cubic')
    avgsinterp = interpolater(yinterp, xinterp).T
    avgsinterp2 = avgsinterp
    avgsinterp2[avgsinterp < np.max(avgsinterp)*threshold] = np.nan
    
    #print(xinterp)
    
    if HH != None and WW == None:
        WW = HH
    if HH == None and WW != None:
        HH = WW
    if HH != None and WW != None:
        mask = np.zeros_like(avgsinterp2)
        H = avgsinterp2.shape[0]
        W = avgsinterp2.shape[1]
        HH = 15
        WW = 15
        #print(W, WW)
        #print(H, HH)
        for i in range(HH, H-HH):
            for j in range(WW, W-WW):
                if np.mean(np.isnan(avgsinterp2[i-HH:i+HH, j-WW])) > 0.99:
                    if np.mean(np.isnan(avgsinterp2[i-HH:i+HH, j+WW])) > 0.99:
                        if np.mean(np.isnan(avgsinterp2[i-HH, j-WW:j+WW])) > 0.99:
                            if np.mean(np.isnan(avgsinterp2[i+HH, j-WW:j+WW])) > 0.99:
                                mask[i-HH:i+HH, j-WW:j+WW] = np.nan
                                #mask[i-int(HH/2):i+int(HH/2), j-int(WW/2):j+int(WW/2)] = np.nan
                                #mask[i, j] = np.nan
        avgsinterp2 += mask
    
    im = ax.imshow(avgsinterp2, cmap = cmap)
    plt.ylim(0, avgsinterp2.shape[0])
    if colorbar:
        plt.colorbar(im)
    if returnxy:
        return avgsinterp2, xinterp, yinterp
    return avgsinterp2

# data is df/f data, and o1, o2, c1, and c2 are the time indices in each arm (open and closed)
def celltype(data, o1, o2, c1, c2):
    o =np.hstack((o1, o2))
    c = np.hstack((c1, c2))
    r_c1o1 = ranksums(data[c1], data[o1])
    r_c1o2 = ranksums(data[c1], data[o2])
    r_c2o1 = ranksums(data[c2], data[o1])
    r_c2o2 = ranksums(data[c2], data[o2])
    if r_c1o1[0] < 0 and r_c1o1[1] < 0.05 and r_c1o2[0] < 0 and r_c1o2[1] < 0.05:
        if r_c2o1[0] < 0 and r_c2o1[1] < 0.05 and r_c2o2[0] < 0 and r_c2o2[1] < 0.05:
            return 'o'
    if r_c1o1[0] > 0 and r_c1o1[1] < 0.05 and r_c1o2[0] > 0 and r_c1o2[1] < 0.05:
        if r_c2o1[0] > 0 and r_c2o1[1] < 0.05 and r_c2o2[0] > 0 and r_c2o2[1] < 0.05:
            return 'c'
    return 'n'

def celltypebehav(data, behav):
    o1 = behav['o1Ind'].flatten()
    o2 = behav['o2Ind'].flatten()
    c1 = behav['c1Ind'].flatten()
    c2 = behav['c2Ind'].flatten()
    o =np.hstack((o1, o2))
    c = np.hstack((c1, c2))
    r_c1o1 = ranksums(data[c1], data[o1])
    r_c1o2 = ranksums(data[c1], data[o2])
    r_c2o1 = ranksums(data[c2], data[o1])
    r_c2o2 = ranksums(data[c2], data[o2])
    if r_c1o1[0] < 0 and r_c1o1[1] < 0.05 and r_c1o2[0] < 0 and r_c1o2[1] < 0.05:
        if r_c2o1[0] < 0 and r_c2o1[1] < 0.05 and r_c2o2[0] < 0 and r_c2o2[1] < 0.05:
            return 'o'
    if r_c1o1[0] > 0 and r_c1o1[1] < 0.05 and r_c1o2[0] > 0 and r_c1o2[1] < 0.05:
        if r_c2o1[0] > 0 and r_c2o1[1] < 0.05 and r_c2o2[0] > 0 and r_c2o2[1] < 0.05:
            return 'c'
    return 'n'

def getcalc(foldnames, mouse, assay, threshold=True, denoise=True, returnIndices=False):
    track, neur, behav = getdicts(foldnames, mouse, assay)
    try:
        calc = neur['C_raw']
        if denoise:
            if len(rmpcas[assay][mouse]) > 0:
                pca = PCA()
                cpca = pca.fit_transform(calc)
                cpca[:, rmpcas[assay][mouse]] = 0
                calc = pca.inverse_transform(cpca)
        if threshold:
            var = np.var(calc, axis = 0)
            n = calc.shape[1]
            order = np.argsort(var)
            vmax = np.max(np.var(calc[:, order[0:n - outliers[assay][mouse]]], axis = 0))
            keepind = (var > 0.1*vmax)
            calc = calc[:, keepind]
        if returnIndices:
            return calc, keepind
        return calc
    except:
        pass
    return

def getmpos(foldnames, mouse, assay):
    track, neur, behav = getdicts(foldnames, mouse, assay)
    return track['mouse_positionMS']

def getcalcs(foldnames, mouse, assay1, assay2, returnInds = False, denoise=True):
    CR = sio.loadmat('../' + mouse + '/coreg.mat')['CR']
    ind1 = Inds[assay1][mouse]
    ind2 = Inds[assay2][mouse]
    if ind1 == None or ind2 == None:
        return False

    np.sum((CR[:, ind1] > 0)*(CR[:, ind2] > 0))
    Cells1 = CR[(CR[:, ind1] != 0) * (CR[:, ind2] != 0), ind1] - 1 #Python indices of shared cells between Rat 2 and Rat1 in Rat 2 assay
    Cells2 = CR[(CR[:, ind1] != 0) * (CR[:, ind2] != 0), ind2] - 1 #vice versa

    track1, neur1, behav1 = getdicts(foldnames, mouse, assay1)
    calc1 = getcalc(foldnames, mouse, assay1, False, denoise=denoise)
    calc1 = calc1[:, Cells1]
    mpos1 = track1['mouse_positionMS']
    track2, neur2, behav2 = getdicts(foldnames, mouse, assay2)
    calc2= getcalc(foldnames, mouse, assay2, False, denoise=denoise)
    calc2 = calc2[:, Cells2]
    mpos2 = track2['mouse_positionMS']
    #print(calc1.shape, calc2.shape)
    n1 = neur1['C_raw'].shape[1]
    n2 = neur2['C_raw'].shape[1]

    keep1 = np.zeros(neur1['C_raw'].shape[1])
    keep2 = np.zeros(neur2['C_raw'].shape[1])
    keep1[Cells1] = 1
    keep2[Cells2] = 1
    var1 = np.var(neur1['C_raw'], axis = 0)
    var2 = np.var(neur2['C_raw'], axis = 0)
    order1 = np.argsort(var1)
    order2 = np.argsort(var2)

    vmax1 = np.max(np.var(neur1['C_raw'][:, order1[0:n1-outliers[assay1][mouse]]], axis = 0))
    vmax2 = np.max(np.var(neur2['C_raw'][:, order2[0:n2-outliers[assay2][mouse]]], axis = 0))
    var1b = np.var(calc1, axis = 0)
    var2b = np.var(calc2, axis = 0)
    
    keep1[var1 < 0.1*vmax1] = 0
    keep2[var2 < 0.1*vmax2] = 0
    
    keep1 = keep1[Cells1]
    keep2 = keep2[Cells2]
    keep = keep1*keep2
    keep = keep.astype('bool')

    calc1 = calc1[:, keep]
    calc2 = calc2[:, keep]
    
    if returnInds:
        return calc1, calc2, (keep, Cells1, Cells2)
    return calc1, calc2

def timenorm(calc, kernel_length = 750, window='rect', keep_var = True):
    if window == 'hamming':
        #technically a half-hamming
        kern = np.hamming(int(2*kernel_length + 1))[-kernel_length:]
        #kern = np.hamming(kernel_length)
    elif window == 'rect' or window == 'rectangle':
        kern = np.ones(kernel_length)
    var = np.var(calc, axis = 0)
    kern /= np.sum(kern)
    chats = np.zeros_like(calc)
    length = calc.shape[0]
    for k in range(calc.shape[1]):
        c = calc[:, k]
        ctilde = c - np.convolve(kern, c, 'full')[0:length]
        cvar = np.convolve((ctilde**2), kern, 'full')[0:length]
        chat = ctilde / np.sqrt(cvar)
        chats[:, k] = chat
    if keep_var == True:
        chats /= np.sqrt(np.var(chats, axis = 0))
        chats *= np.sqrt(var)
    return chats

def getAUC(scores, labels, returnTPFP=False):
    scores = scores.flatten()
    labels = labels.flatten()
    labels = labels.astype('bool')

    order = np.flip(np.argsort(scores))
    labels_sort = labels[order]
    num1 = np.sum(labels)
    num0 = labels.size - num1
    
    tprs = np.zeros(scores.size + 1)
    fprs = scores.size*np.zeros(scores.size + 1)
    for i in range(1, scores.size):
        tprs[i] = tprs[i - 1]
        fprs[i] = fprs[i - 1]
        if labels_sort[i - 1]:
            tprs[i] += 1
        else:
            fprs[i] += 1
    tprs /= num1
    fprs /= num0
    tprs[-1] = 1
    fprs[-1] = 1

    auc = 0
    for i in range(1, tprs.size):
        auc += (fprs[i] - fprs[i-1])*(0.5*tprs[i] + 0.5*tprs[i])
    if returnTPFP:
        return auc, np.vstack((tprs, fprs)).T
    return auc

def armscore(data, o1, o2, c1, c2):
    for ind in [o1, o2, c1, c2]:
        ind = ind.flatten()
    oind = np.hstack((o1, o2))
    cind = np.hstack((c1, c2))
    x = np.hstack((data[oind], data[cind]))
    y = np.hstack((np.ones(oind.size), np.zeros(cind.size)))
    auc = getAUC(x, y)
    return 2*(auc - 0.5)

def armscorebehav(data, behav, bootstrap = False):
    o1 = behav['o1Ind'].flatten()
    o2 = behav['o2Ind'].flatten()
    c1 = behav['c1Ind'].flatten()
    c2 = behav['c2Ind'].flatten()
    oind = np.hstack((o1, o2))
    cind = np.hstack((c1, c2))
    x = np.hstack((data[oind], data[cind]))
    y = np.hstack((np.ones(oind.size), np.zeros(cind.size)))
    if bootstrap == True:
        y = y[np.random.permutation(y.size)]
    auc = getAUC(x, y)
    return 2*(auc - 0.5)

def rescale_epm_mpos(mpos, behav):
    #mpos = np.array(mpos)
    o1 = behav['o1Ind']
    o2 = behav['o2Ind']
    c1 = behav['c1Ind']
    c2 = behav['c2Ind']
    centerInd = behav['centerInd']
    
    centerx = (np.min(mpos[centerInd, 0]), np.max(mpos[centerInd, 0]))
    centery = (np.min(mpos[centerInd, 1]), np.max(mpos[centerInd, 1]))
    mpos[:, 0] = mpos[:, 0] - np.mean(centerx)
    mpos[:, 1] = mpos[:, 1] - np.mean(centery)
    
    xlim = (np.min(mpos[:, 0]), np.max(mpos[:, 0]))
    ylim = (np.min(mpos[:, 1]), np.max(mpos[:, 1]))
    
    mpos[mpos[:, 0] >= 0, 0] = mpos[mpos[:, 0] >= 0, 0] / xlim[1]
    mpos[mpos[:, 0] < 0, 0] = mpos[mpos[:, 0] < 0, 0] / np.abs(xlim[0])
    mpos[mpos[:, 1] >= 0, 1] = mpos[mpos[:, 1] >= 0, 1] / ylim[1]
    mpos[mpos[:, 1] < 0, 1] = mpos[mpos[:, 1] < 0, 1] / np.abs(ylim[0])
    
    return mpos

def epmthreat(mpos, behav):
    mpos = rescale_epm_mpos(mpos, behav)
    return np.abs(mpos[:, 0]) - np.abs(mpos[:, 1])

def minmaxscale(data_in, vmin = 0, vmax = 1, refmax = None):
    data = data_in.copy()
    dmin = np.min(data, axis = 0)
    dmax = np.max(data, axis = 0)
    data -= dmin
    if refmax != None:
        data /= (refmax - dmin)
    else:
        data /= np.max(data, axis = 0)
    
    data *= (vmax - vmin)
    data += vmin
    return data

def getinterleavedinds(stepsize, gap, L):
    i = 0
    toggle = True
    inds1 = np.zeros(0, dtype='int16')
    inds2 = np.zeros(0, dtype='int16')
    while i < L:
        inds = np.arange(i, np.minimum(i + stepsize - gap, L), dtype='int16')
        if toggle:
            inds1 = np.hstack((inds1, inds))
            toggle = False
        else:
            inds2 = np.hstack((inds2, inds))
            toggle = True
        i += stepsize
    return inds1, inds2

def get_crouchsniff(foldnames, mousenum, assay):
    fold = getfold(foldnames, mousenum, assay)
    try:
        d = sio.loadmat(fold + '\\Behavior_Rear_CrouchSniff_Groom_MS.mat')
        L = getcalc(foldnames, mousenum, assay).shape[0]
        for item in list(d):
            if 'Frame' in item:
                d[item] -= 3
                while np.min(d[item][0]) < 0:
                    d[item] = d[item][1:]
            elif 'Indices' in item:
                d[item] = d[item][3:L+3]
    except:
        return
    return d
##################
foldnames = getfoldnames('../')
mousenums = getmousenums('../')
print('mousenums: ', mousenums)
##################