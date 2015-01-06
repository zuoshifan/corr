import numpy as np
import h5py

N = 768 # grid size in one dimension
num = 8 # contract number
# r0 = 1   # grid unit


def contract(den):
    xs, ys, zs = den.shape
    assert xs == ys == zs, 'Unequal size for all dimensions'
    assert xs % 2 == 0, 'Can not divide by 2'
    new_den = np.zeros((xs/2, xs/2, xs/2), dtype=den.dtype)
    for x in range(0, xs, 2):
        for y in range(0, xs, 2):
            for z in range(0, xs, 2):
                new_den[x/2, y/2, z/2] = (den[x, y, z] + den[x+1, y, z] + den[x, y+1, z] + den[x, y, z+1] + den[x+1, y+1, z] + den[x+1, y, z+1] + den[x, y+1, z+1] + den[x+1, y+1, z+1]) / 8.0

    return new_den


# read density data from file
den = np.fromfile('den_dm00.dat', dtype=np.float32)
den = den.reshape((N, N, N)).T # transpose because data is fortran array
den = den - 1.0 # original data is rho/rho_bar

r = []
corr = []
for count in range(num+1):
    print '%d of %d...' % (count, num+1)
    r.append(2**count)
    r.append(np.sqrt(2) * (2**count))
    r.append(np.sqrt(3) * (2**count))
    if count != 0:
        den = contract(den)
    Nden = den.shape[0]
    corr1 = 0.0
    corr2 = 0.0 # sqrt(2)
    corr3 = 0.0 # sqrt(3)
    for i in range(1, Nden-1):
        for j in range(1, Nden-1):
            for k in range(1, Nden-1):
                corr1 += den[i, j, k] * (den[i-1, j, k] + den[i+1, j, k] + den[i, j-1, k] + den[i, j+1, k] + den[i, j, k-1] + den[i, j, k+1]) / 6.0
                corr2 += den[i, j, k] * (den[i-1, j-1, k] + den[i-1, j+1, k] + den[i-1, j, k-1] + den[i-1, j, k+1] + den[i+1, j-1, k] + den[i+1, j+1, k] + den[i+1, j, k-1] + den[i+1, j, k+1] + den[i, j-1, k-1] + den[i, j+1, k-1] + den[i, j-1, k+1] + den[i, j+1, k+1]) / 12.0
                corr3 += den[i, j, k] * (den[i-1, j-1, k-1] + den[i-1, j-1, k+1] + den[i-1, j+1, k-1] + den[i-1, j+1, k+1] + den[i+1, j-1, k-1] + den[i+1, j-1, k+1] + den[i+1, j+1, k-1] + den[i+1, j+1, k+1]) / 8.0

    corr.append(corr1/(Nden-2)**3)
    corr.append(corr2/(Nden-2)**3)
    corr.append(corr3/(Nden-2)**3)


with h5py.File('corr.hdf5', 'w') as f:
    f.create_dataset('r', data=np.array(r))
    f.create_dataset('corr', data=np.array(corr))