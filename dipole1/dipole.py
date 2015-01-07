import numpy as np
import h5py
import correlate


N = 768 # grid size in one dimension
num = 8 # contract number
# r0 = 1   # grid unit


# read density data from file
# dark mater data
dm_den = np.fromfile('den_dm00.dat', dtype=np.float32)
dm_den = dm_den.reshape((N, N, N)).T # transpose because data is fortran array, now three dimensions are [x, y, z]
dm_den = dm_den - 1.0 # original data is rho/rho_bar

# dark mater velocity
dm_vx = np.fromfile('vx_dm.bin', dtype=np.float32)
dm_vx = dm_vx.reshape((N, N, N)).T
dm_vy = np.fromfile('vy_dm.bin', dtype=np.float32)
dm_vy = dm_vy.reshape((N, N, N)).T
dm_vz = np.fromfile('vz_dm.bin', dtype=np.float32)
dm_vz = dm_vz.reshape((N, N, N)).T

# neutrino data
nu_den = np.fromfile('den_nu00.dat', dtype=np.float32)
nu_den = nu_den.reshape((N, N, N)).T # transpose because data is fortran array, now three dimensions are [x, y, z]
nu_den = nu_den - 1.0 # original data is rho/rho_bar

# neutrino velocity
nu_vx = np.fromfile('vx_nu.bin', dtype=np.float32)
nu_vx = nu_vx.reshape((N, N, N)).T
nu_vy = np.fromfile('vy_nu.bin', dtype=np.float32)
nu_vy = nu_vy.reshape((N, N, N)).T
nu_vz = np.fromfile('vz_nu.bin', dtype=np.float32)
nu_vz = nu_vz.reshape((N, N, N)).T

# relative velocity between dark mater and neutrino
vx = nu_vx - dm_vx
vy = nu_vy - dm_vy
vz = nu_vz - dm_vz
# velocity normalization
vx, vy, vz = correlate.norm(vx, vy, vz)

# release no usefull arrays
del dm_vx
del dm_vy
del dm_vz
del nu_vx
del nu_vy
del nu_vz


# ##----  test code ----------------------------------
# N = 12
# num = 2
# dm_den = np.arange(12**3, dtype=np.float32).reshape(12, 12, 12).T
# nu_den = np.arange(12**3, dtype=np.float32).reshape(12, 12, 12).T
# vx = np.arange(12**3, dtype=np.float32).reshape(12, 12, 12).T
# vy = np.arange(12**3, dtype=np.float32).reshape(12, 12, 12).T
# vz = np.arange(12**3, dtype=np.float32).reshape(12, 12, 12).T
# vx, vy, vz = correlate.norm(vx, vy, vz)
# ##----  test code ----------------------------------


print 'Start to compute...'
r, corr = correlate.compute(dm_den, nu_den, vx, vy, vz, N, num)

with h5py.File('corr.hdf5', 'w') as f:
    f.create_dataset('r', data=np.array(r))
    f.create_dataset('corr', data=np.array(corr))
