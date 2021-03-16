import numpy as np
import torch
import torch.nn as nn
import copy
from sklearn.cross_decomposition import CCA
import pdb

def correlation(Z1, Z2):
    if Z1 is None or Z2 is None:
        return 0
    if Z1.shape[0] != Z2.shape[0]:
        raise ValueError('Data must be the same length')
    if Z1.shape[0] == 0:
        return 0
    Z1 = Z1 - torch.mean(Z1)
    Z2 = Z2 - torch.mean(Z2)
    
    cov = torch.dot(Z1[:, 0], Z2[:, 0]) / Z1.shape[0]
    std1 = torch.sqrt(torch.var(Z1))
    std2 = torch.sqrt(torch.var(Z2))
    return (cov/(std1*std2))

class Projector(nn.Module):
    def __init__(self, xdims, ydims, bias=True, devicename='cpu'):
        super(Projector, self).__init__()
        self.xdims = xdims
        self.ydims = ydims
        self.bias = bias
        self.device = torch.device(devicename)
        
        self.xWeights = nn.ParameterList(
            [nn.Parameter(2*torch.rand(dim, 1)/np.sqrt(dim)-1/np.sqrt(dim)) for dim in self.xdims]
        )
        self.yWeights = nn.ParameterList(
            [nn.Parameter(2*torch.rand(dim, 1)/np.sqrt(dim)-1/np.sqrt(dim)) for dim in self.ydims]
        )
        
        if bias:
            self.xBiases = nn.ParameterList(
                [nn.Parameter(2*torch.rand(1)/np.sqrt(dim)-1/np.sqrt(dim)) for dim in self.xdims]
            )
            self.yBiases = nn.ParameterList(
                [nn.Parameter(2*torch.rand(1)/np.sqrt(dim)-1/np.sqrt(dim)) for dim in self.ydims]
            )
        return
    
    def forward(self, X, Y, idx):
        # idx is a 2-tuple of the corresponding xModel and yModel
        if X is None:
            outx = None
        else:
            X = torch.from_numpy(X).float().to(self.device)
            outx = X@self.xWeights[idx[0]]
            if self.bias:
                outx += self.xBiases[idx[0]]
        if Y is None:
            outy = None
        else:
            Y = torch.from_numpy(Y).float().to(self.device)
            outy = Y@self.yWeights[idx[1]]
            if self.bias:
                outy += self.yBiases[idx[1]]
        return outx, outy

class CoCA:
    def __init__(self, n_components=1, lr=1e-2, scale=True, max_iter=10000, tol=1e-7, 
                 weight_decay=0, metric='sum', cca_init=True, devicename='cpu'):
        self.lr = lr
        self.scale = scale
        self.max_iter = max_iter
        self.tol = tol
        self.metric = metric
        self.weight_decay = weight_decay
        self.n_components = n_components
        self.cca_init = cca_init
        self.device = torch.device(devicename)
        self.devicename = devicename
    
    def fit(self, Xs, Ys, XYweights=None):
        # self.randomstate = np.random.get_state()
        Xs = copy.deepcopy(Xs)
        Ys = copy.deepcopy(Ys)

        if XYweights is None:
            XYweights = torch.ones(Xs.shape)  # use only for metric = 'sum' and 'sumsq'
        else:
            XYweights = torch.from_numpy(XYweights).float()
        xdims = []
        ydims = []
        for idx_x in range(Xs.shape[0]):
            for idx_y in range(Ys.shape[1]):
                X = Xs[idx_x, idx_y]
                if X is not None:
                    if X.shape[0] > 0:
                        xdims.append(X.shape[1])
                    break
            else:
                raise ValueError('All X rows must be non-empty')
        for idx_y in range(Ys.shape[1]):
            for idx_x in range(Xs.shape[0]):
                Y = Ys[idx_x, idx_y]
                if Y is not None:
                    if Y.shape[0] > 0:
                        ydims.append(Y.shape[1])
                    break
            else:
                raise ValueError('All Y columns must be non-empty')
        self.xdims = xdims
        self.ydims = ydims
        
        if self.scale:
            # note that mean and variance normalization is done independently of XYweights
            self.muX_ = np.array(
                [np.mean([np.mean(Xs[idx_x, idx_y], axis=0) for idx_y in range(Xs.shape[1]) \
                          if Xs[idx_x, idx_y] is not None], axis=0) for idx_x in range(Xs.shape[0])],
                dtype=list
            )
            self.muY_ = np.array(
                [np.mean([np.mean(Ys[idx_x, idx_y], axis=0) for idx_x in range(Ys.shape[0]) \
                          if Ys[idx_x, idx_y] is not None], axis=0) for idx_y in range(Ys.shape[1])],
                dtype=list
            )
            self.stdX_ = np.array(
                [np.sqrt(np.mean([np.var(Xs[idx_x, idx_y], axis=0) for idx_y in range(Xs.shape[1]) \
                                  if Xs[idx_x, idx_y] is not None], axis=0)) for idx_x in range(Xs.shape[0])],
                dtype=list
            )
            self.stdY_ = np.array(
                [np.sqrt(np.mean([np.var(Ys[idx_x, idx_y], axis=0) for idx_x in range(Ys.shape[0]) \
                                  if Ys[idx_x, idx_y] is not None], axis=0)) for idx_y in range(Ys.shape[1])],
                dtype=list
            )
        else:
            self.muX_ = np.zeros(Xs.shape[0])
            self.muY_ = np.zeros(Ys.shape[1])
            self.stdX_ = np.ones(Xs.shape[0])
            self.stdY_ = np.ones(Ys.shape[1])
        
        for idx_x in range(Xs.shape[0]):
            for idx_y in range(Ys.shape[1]):
                if Xs[idx_x, idx_y] is not None:
                    Xs[idx_x, idx_y] = (Xs[idx_x, idx_y] - self.muX_[idx_x]) / self.stdX_[idx_x]
                else:
                    XYweights[idx_x, idx_y] = 0 # ignore missing data
                if Ys[idx_x, idx_y] is not None:
                    Ys[idx_x, idx_y] = (Ys[idx_x, idx_y] - self.muY_[idx_y]) / self.stdY_[idx_y]
                else:
                    XYweights[idx_x, idx_y] = 0
        
        self.x_weights_ = np.array([np.zeros((dim, self.n_components)) for dim in xdims], dtype=list)
        self.y_weights_ = np.array([np.zeros((dim, self.n_components)) for dim in ydims], dtype=list)
        
        self.n_iter_ = np.zeros(self.n_components)
        for component in range(self.n_components):
            model = Projector(xdims, ydims, bias=False, devicename=self.devicename)
            
            if self.cca_init:
                weight_ccax = np.empty((Xs.shape[0], Ys.shape[1]), dtype=list)
                weight_ccay = np.empty((Xs.shape[0], Ys.shape[1]), dtype=list)
                for idx_x in range(Xs.shape[0]):
                    for idx_y in range(Ys.shape[1]):
                        if Xs[idx_x, idx_y] is not None and Ys[idx_x, idx_y] is not None:
                            cca = CCA(n_components=1, scale=self.scale)
                            cca.fit(Xs[idx_x, idx_y], Ys[idx_x, idx_y])
                            weight_ccax[idx_x, idx_y] = cca.x_weights_[:, 0]/cca.x_std_
                            weight_ccay[idx_x, idx_y] = cca.y_weights_[:, 0]/cca.y_std_
                
                with torch.no_grad():
                    for idx_x in range(Xs.shape[0]):
                        weights = np.array([w for w in weight_ccax[idx_x, :] if w is not None])
                        model.xWeights[idx_x].copy_(torch.FloatTensor(np.mean(weights, axis=0, keepdims=True).T))
                    for idx_y in range(Ys.shape[1]):
                        weights = np.array([w for w in weight_ccay[:, idx_y] if w is not None])
                        model.yWeights[idx_y].copy_(torch.FloatTensor(np.mean(weights, axis=0, keepdims=True).T))
            
            model.to(self.device)
            self.optimizer = torch.optim.Adam(model.parameters(), lr=self.lr, 
                                              weight_decay=self.weight_decay)

            last_loss = np.nan
            for idx_iter in range(self.max_iter):
                self.optimizer.zero_grad()
                iter_loss = 0
                corrs = torch.zeros((len(xdims), len(ydims)))
                for idx_x in range(len(xdims)):
                    for idx_y in range(len(ydims)):
                        X = Xs[idx_x, idx_y]
                        Y = Ys[idx_x, idx_y]
                        Zx, Zy = model(X, Y, (idx_x, idx_y))
                        corrs[idx_x, idx_y] = correlation(Zx, Zy)

                if self.metric == 'sum':
                    loss = -torch.sum(XYweights*corrs)
                elif self.metric == 'sumsq':
                    loss = -torch.sum(XYweights*torch.sign(corrs)*corrs**2)
                elif self.metric == 'min':
                    corrs += 1*(XYweights == 0)
                    # prevents treating missing data as having r=0
                    # ignores also XYweights=0
                    loss = -torch.min(corrs)
                elif self.metric == 'minsq':
                    # an alternative to 'min', but may perform worse
                    corrs += 1*(XYweights == 0)
                    loss = -torch.min(torch.sign(corrs)*corrs**2)
                elif self.metric == 'summinX':
                    corrs += 1*(XYweights == 0)
                    loss = -torch.sum(torch.min(corrs, axis=1)[0])
                elif self.metric == 'summinY':
                    corrs += 1*(XYweights == 0)
                    loss = -torch.sum(torch.min(corrs, axis=0)[0])
                elif self.metric == 'summinsqX':
                    corrs += 1*(XYweights == 0)
                    loss = -torch.sum(torch.min(torch.sign(corrs)*corrs**2, axis=1)[0])
                elif self.metric == 'summinsqY':
                    corrs += 1*(XYweights == 0)
                    loss = -torch.sum(torch.min(torch.sign(corrs)*corrs**2, axis=0)[0])
                else:
                    print('invalid metric')
                    return
                iter_loss = loss.item()
                loss.backward()
                self.optimizer.step()
                
                if np.abs(iter_loss - last_loss) < self.tol:
                    break
                last_loss = iter_loss
            self.n_iter_[component] = idx_iter + 1
            self.final_loss_ = last_loss
            
            for i, weight in enumerate(model.xWeights):
                v = weight.detach().cpu().numpy().flatten()
                dim = self.xdims[i]
                P = np.zeros((dim, dim))
                for j in range(component):
                    vj = self.x_weights_[i][:, j].reshape((-1, 1))
                    P += vj@vj.T
                u = v - P@v
                self.x_weights_[i][:, component] = u / np.linalg.norm(u)
            for i, weight in enumerate(model.yWeights):
                v = weight.detach().cpu().numpy().flatten()
                dim = self.ydims[i]
                P = np.zeros((dim, dim))
                for j in range(component):
                    vj = self.y_weights_[i][:, j].reshape((-1, 1))
                    P += vj@vj.T
                u = v - P@v
                self.y_weights_[i][:, component] = u / np.linalg.norm(u)
            
            for idx_x in range(len(xdims)):
                for idx_y in range(len(ydims)):
                    if Xs[idx_x, idx_y] is not None:
                        ux = self.x_weights_[idx_x][:, component].reshape((-1, 1)).astype('float64')
                        Xs[idx_x, idx_y] -= (Xs[idx_x, idx_y]@ux)@ux.T
                    if Ys[idx_x, idx_y] is not None:
                        uy = self.y_weights_[idx_y][:, component].reshape((-1, 1)).astype('float64')
                        #pdb.set_trace()
                        Ys[idx_x, idx_y] -= (Ys[idx_x, idx_y]@uy)@uy.T
        return
        
    def transform(self, Xs, Ys):
        xdims = self.xdims
        ydims = self.ydims
        Zxs = np.empty((len(xdims), len(ydims)), dtype=list)
        Zys = np.empty((len(xdims), len(ydims)), dtype=list)
        
        for idx_x in range(len(xdims)):
            for idx_y in range(len(ydims)):
                X = Xs[idx_x, idx_y]
                Y = Ys[idx_x, idx_y]
                if X is not None:
                    Zxs[idx_x, idx_y] = ((X - self.muX_[idx_x]) / self.stdX_[idx_x])@self.x_weights_[idx_x]
                if Y is not None:
                    Zys[idx_x, idx_y] = ((Y - self.muY_[idx_y]) / self.stdY_[idx_y])@self.y_weights_[idx_y]
        return Zxs, Zys
    
    def fit_transform(self, Xs, Ys, XYweights=None):
        self.fit(Xs, Ys, XYweights)
        return self.transform(Xs, Ys)
    
    def get_weights(self):
        return self.x_weights_, self.y_weights_