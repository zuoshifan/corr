import numpy as np
import h5py


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


N = 768 # grid size in one dimension
num = 8 # contract number
# r0 = 1   # grid unit


# read density data from file
# dark mater data
dm_den = np.fromfile('den_dm00.dat', dtype=np.float32)
dm_den = dm_den.reshape((N, N, N)).T # transpose because data is fortran array
dm_den = dm_den - 1.0 # original data is rho/rho_bar

# neutrino data
nu_den = np.fromfile('den_nu00.dat', dtype=np.float32)
nu_den = nu_den.reshape((N, N, N)).T # transpose because data is fortran array
nu_den = nu_den - 1.0 # original data is rho/rho_bar

# ##----  test code ----------------------------------
# N = 12
# num = 2
# dm_den = np.arange(12**3).reshape(12, 12, 12).T
# nu_den = np.arange(12**3).reshape(12, 12, 12).T
# ##----  test code ----------------------------------

r = []
corr = []
for count in range(num+1):
    print '%d of %d...' % (count, num+1)
    r.append(2**count)
    r.append(np.sqrt(2) * (2**count))
    r.append(np.sqrt(3) * (2**count))
    if count != 0:
        dm_den = contract(dm_den)
        nu_den = contract(nu_den)
    Nden = dm_den.shape[0]
    corr1 = 0.0
    corr2 = 0.0 # sqrt(2)
    corr3 = 0.0 # sqrt(3)
    for i in range(1, Nden-1):
        for j in range(1, Nden-1):
            for k in range(1, Nden-1):
                corr1 += dm_den[i, j, k] * (nu_den[i-1, j, k] + nu_den[i+1, j, k] + nu_den[i, j-1, k] + nu_den[i, j+1, k] + nu_den[i, j, k-1] + nu_den[i, j, k+1]) / 6.0
                corr2 += dm_den[i, j, k] * (nu_den[i-1, j-1, k] + nu_den[i-1, j+1, k] + nu_den[i-1, j, k-1] + nu_den[i-1, j, k+1] + nu_den[i+1, j-1, k] + nu_den[i+1, j+1, k] + nu_den[i+1, j, k-1] + nu_den[i+1, j, k+1] + nu_den[i, j-1, k-1] + nu_den[i, j+1, k-1] + nu_den[i, j-1, k+1] + nu_den[i, j+1, k+1]) / 12.0
                corr3 += dm_den[i, j, k] * (nu_den[i-1, j-1, k-1] + nu_den[i-1, j-1, k+1] + nu_den[i-1, j+1, k-1] + nu_den[i-1, j+1, k+1] + nu_den[i+1, j-1, k-1] + nu_den[i+1, j-1, k+1] + nu_den[i+1, j+1, k-1] + nu_den[i+1, j+1, k+1]) / 8.0

    corr.append(corr1/(Nden-2)**3)
    corr.append(corr2/(Nden-2)**3)
    corr.append(corr3/(Nden-2)**3)


with h5py.File('corr.hdf5', 'w') as f:
    f.create_dataset('r', data=np.array(r))
    f.create_dataset('corr', data=np.array(corr))