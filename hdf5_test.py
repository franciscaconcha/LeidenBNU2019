import numpy as np
import h5py


A = np.arange(10)

data = h5py.File('test.hdf5', 'w')
data.create_dataset('A', data=A)
data.close()


data = h5py.File('test.hdf5', 'r')
B = np.array(data['A'])
print B
