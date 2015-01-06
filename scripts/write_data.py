import numpy as np
import h5py

with h5py.File('corr.hdf5', 'r') as f:
    r = f['r'][:]
    corr = f['corr'][:]

# r.dtype = np.float32
# r.tofile('r.dat')
# corr.dtype = np.float32
# corr.tofile('corr.dat')

# r1 =np.fromfile('r.dat', dtype=np.float32)
# assert np.allclose(r, r1)
# print 'Done.'

f1=open('./corr.txt', 'w+')
for i in range(len(r)):
    print >>f1, '%.12f    %.12f' % (r[i], corr[i])