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


def norm(vx, vy, vz):
    v = np.sqrt(vx**2 + vy**2 + vz**2)
    return vx/v, vy/v, vz/v


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
vx, vy, vz = norm(vx, vy, vz)

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
# dm_den = np.arange(12**3).reshape(12, 12, 12).T
# nu_den = np.arange(12**3).reshape(12, 12, 12).T
# vx = np.arange(12**3).reshape(12, 12, 12).T
# vy = np.arange(12**3).reshape(12, 12, 12).T
# vz = np.arange(12**3).reshape(12, 12, 12).T
# vx, vy, vz = norm(vx, vy, vz)
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
        vx = contract(vx)
        vy = contract(vy)
        vz = contract(vz)
        vx, vy, vz = norm(vx, vy, vz)
    Nden = dm_den.shape[0]
    corr1 = 0.0
    corr2 = 0.0 # sqrt(2)
    corr3 = 0.0 # sqrt(3)
    for x in range(1, Nden-1):
        for y in range(1, Nden-1):
            for z in range(1, Nden-1):
                corr1 += dm_den[x, y, z] * (nu_den[x-1, y, z] * (-1.0 * vx[x, y, z]) + nu_den[x+1, y, z] * (1.0 * vx[x, y, z]) + nu_den[x, y-1, z] * (-1.0 * vy[x, y, z]) + nu_den[x, y+1, z] * (1.0 * vy[x, y, z]) + nu_den[x, y, z-1] * (-1.0 * vz[x, y, z]) + nu_den[x, y, z+1] * (1.0 * vz[x, y, z])) / 6.0
                corr2 += dm_den[x, y, z] * (nu_den[x-1, y-1, z] * (-1.0 * vx[x, y, z] - 1.0 * vy[x, y, z]) + nu_den[x-1, y+1, z] * (-1.0 * vx[x, y, z] + 1.0 * vy[x, y, z]) + nu_den[x-1, y, z-1] * (-1.0 * vx[x, y, z] - 1.0 * vz[x, y, z]) + nu_den[x-1, y, z+1] * (-1.0 * vx[x, y, z] + 1.0 * vz[x, y, z]) + nu_den[x+1, y-1, z] * (1.0 * vx[x, y, z] - 1.0 * vy[x, y, z]) + nu_den[x+1, y+1, z] * (1.0 * vx[x, y, z] + 1.0 * vy[x, y, z]) + nu_den[x+1, y, z-1] * (1.0 * vx[x, y, z] - 1.0 * vz[x, y, z]) + nu_den[x+1, y, z+1] * (1.0 * vx[x, y, z] + 1.0 * vz[x, y, z]) + nu_den[x, y-1, z-1] * (-1.0 * vy[x, y, z] - 1.0 * vz[x, y, z] ) + nu_den[x, y+1, z-1] * (1.0 * vy[x, y, z] - 1.0 * vz[x, y, z]) + nu_den[x, y-1, z+1] * (-1.0 * vy[x, y, z] + 1.0 * vz[x, y, z]) + nu_den[x, y+1, z+1] * (1.0 * vy[x, y, z] + 1.0 * vz[x, y, z])) / (12.0 * np.sqrt(2))
                corr3 += dm_den[x, y, z] * (nu_den[x-1, y-1, z-1] * (-1.0 * vx[x, y, z] - 1.0 * vy[x, y, z] - 1.0 * vz[x, y, z]) + nu_den[x-1, y-1, z+1] * (-1.0 * vx[x, y, z] - 1.0 * vy[x, y, z] + 1.0 * vz[x, y, z]) + nu_den[x-1, y+1, z-1] * (-1.0 * vx[x, y, z] + 1.0 * vy[x, y, z] - 1.0 * vz[x, y, z]) + nu_den[x-1, y+1, z+1] * (-1.0 * vx[x, y, z] + 1.0 * vy[x, y, z] + 1.0 * vz[x, y, z]) + nu_den[x+1, y-1, z-1] * (1.0 * vx[x, y, z] - 1.0 * vy[x, y, z] - 1.0 * vz[x, y, z]) + nu_den[x+1, y-1, z+1] * (1.0 * vx[x, y, z] - 1.0 * vy[x, y, z] + 1.0 * vz[x, y, z]) + nu_den[x+1, y+1, z-1] * (1.0 * vx[x, y, z] + 1.0 * vy[x, y, z] - 1.0 * vz[x, y, z]) + nu_den[x+1, y+1, z+1] * (1.0 * vx[x, y, z] + 1.0 * vy[x, y, z] + 1.0 * vz[x, y, z])) / (8.0 * np.sqrt(3))

    corr.append(corr1/(Nden-2)**3)
    corr.append(corr2/(Nden-2)**3)
    corr.append(corr3/(Nden-2)**3)


with h5py.File('corr.hdf5', 'w') as f:
    f.create_dataset('r', data=np.array(r))
    f.create_dataset('corr', data=np.array(corr))
