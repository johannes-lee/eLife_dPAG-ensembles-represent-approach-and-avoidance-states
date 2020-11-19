from helpers import *
import hdf5storage

mousenums = getmousenums('../')
os.chdir('../')

for mouse in mousenums:
    os.chdir(mouse)
    for filename in os.listdir():
        if filename.find('cellRegistered') > -1:
            CR = hdf5storage.loadmat(filename)['cell_registered_struct']['cell_to_index_map'][0].astype('int')
            #ind = np.arange(CR.shape[1])
            #ind[1] = 0
            #ind[0] = 1
            #
            #if CR.shape[1] == 7:
            #    CR = CR[:, [1, 0, 2, 3, 4, 5, 6]]
            #elif CR.shape[1] == 6:
            #    CR = CR[:, [1, 0, 2, 3, 4, 5]]
            crdict = {'CR': CR}
            print(mouse)
            sio.savemat('coreg.mat', crdict)
    os.chdir('../')