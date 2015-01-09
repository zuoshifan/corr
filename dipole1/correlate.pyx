import sys
import numpy as np
cimport numpy as np


# We now need to fix a datatype for our arrays. I've used the variable
# DTYPE for this, which is assigned to the usual NumPy runtime
# type info object.
DTYPE = np.float32
# "ctypedef" assigns a corresponding compile-time type to DTYPE_t. For
# every type in the numpy module there's a corresponding compile-time
# type with a _t-suffix.
ctypedef np.float32_t DTYPE_t

cimport cython
@cython.boundscheck(False) # turn of bounds-checking for entire function
def contract(np.ndarray[DTYPE_t, ndim=3] den, unsigned int Nden):
    cdef unsigned int x, y, z
    cdef np.ndarray[DTYPE_t, ndim=3] new_den
    assert Nden % 2 == 0, 'Can not divide by 2'
    new_den = np.zeros((Nden/2, Nden/2, Nden/2), dtype=DTYPE)
    for x in range(0, Nden, 2):
        for y in range(0, Nden, 2):
            for z in range(0, Nden, 2):
                new_den[x/2, y/2, z/2] = (den[x, y, z] + den[x+1, y, z] + den[x, y+1, z] + den[x, y, z+1] + den[x+1, y+1, z] + den[x+1, y, z+1] + den[x, y+1, z+1] + den[x+1, y+1, z+1]) / 8.0

    return new_den


@cython.boundscheck(False) # turn of bounds-checking for entire function
def contractv(np.ndarray[DTYPE_t, ndim=4] den, unsigned int Nden):
    cdef unsigned int i, x, y, z
    cdef np.ndarray[DTYPE_t, ndim=4] new_den
    assert Nden % 2 == 0, 'Can not divide by 2'
    new_den = np.zeros((Nden/2, Nden/2, Nden/2, 3), dtype=DTYPE)
    for i in range(3):
        for x in range(0, Nden, 2):
            for y in range(0, Nden, 2):
                for z in range(0, Nden, 2):
                    new_den[x/2, y/2, z/2, i] = (den[x, y, z, i] + den[x+1, y, z, i] + den[x, y+1, z, i] + den[x, y, z+1, i] + den[x+1, y+1, z, i] + den[x+1, y, z+1, i] + den[x, y+1, z+1, i] + den[x+1, y+1, z+1, i]) / 8.0

    return new_den



@cython.boundscheck(False) # turn of bounds-checking for entire function
def norm(np.ndarray[DTYPE_t, ndim=1] v):
    v_mag = np.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    return v / v_mag


@cython.boundscheck(False) # turn of bounds-checking for entire function
cdef DTYPE_t dot(DTYPE_t x, DTYPE_t y, DTYPE_t z, np.ndarray[DTYPE_t, ndim=1] v):
    return x * v[0] + y * v[1] + z * v[2]



# @cython.boundscheck(False) # turn of bounds-checking for entire function
# def norm(np.ndarray[DTYPE_t, ndim=3] vx, np.ndarray[DTYPE_t, ndim=3] vy, np.ndarray[DTYPE_t, ndim=3] vz):
#     cdef np.ndarray[DTYPE_t, ndim=3] v
#     v = np.sqrt(vx**2 + vy**2 + vz**2)
#     return vx/v, vy/v, vz/v


@cython.boundscheck(False) # turn of bounds-checking for entire function
def compute(np.ndarray[DTYPE_t, ndim=3] dm_den, np.ndarray[DTYPE_t, ndim=3] nu_den, np.ndarray[DTYPE_t, ndim=4] v, unsigned int Nden, unsigned int num):
    cdef unsigned int count, x, y, z
    cdef float corr1, corr2, corr3
    r = []
    corr = []
    for count in range(num+1):
        print '%d of %d...' % (count, num+1)
        sys.stdout.flush()
        r.append(2**count)
        r.append(np.sqrt(2) * (2**count))
        r.append(np.sqrt(3) * (2**count))
        if count != 0:
            dm_den = contract(dm_den, Nden)
            nu_den = contract(nu_den, Nden)
            v = contractv(v, Nden)
            # vx = contract(vx, Nden)
            # vy = contract(vy, Nden)
            # vz = contract(vz, Nden)
            # vx, vy, vz = norm(vx, vy, vz)
            Nden /= 2
        # Nden = dm_den.shape[0]
        corr1 = 0.0
        corr2 = 0.0 # sqrt(2)
        corr3 = 0.0 # sqrt(3)
        for x in range(1, Nden-1):
            for y in range(1, Nden-1):
                for z in range(1, Nden-1):
                    corr1 += dm_den[x, y, z] * (nu_den[x-1, y, z] * dot(-1.0, 0.0, 0.0, norm(v[x, y, z])) + nu_den[x+1, y, z] * dot(1.0, 0.0, 0.0, norm(v[x, y, z])) + nu_den[x, y-1, z] * dot(0.0, -1.0, 0.0, norm(v[x, y, z])) + nu_den[x, y+1, z] * dot(0.0, 1.0, 0.0, norm(v[x, y, z])) + nu_den[x, y, z-1] * dot(0.0, 0.0, -1.0, norm(v[x, y, z])) + nu_den[x, y, z+1] * dot(0.0, 0.0, 1.0, norm(v[x, y, z]))) / 6.0
                    corr2 += dm_den[x, y, z] * (nu_den[x-1, y-1, z] * dot(-1.0, -1.0, 0.0, norm(v[x, y, z])) + nu_den[x-1, y+1, z] * dot(-1.0, 1.0, 0.0, norm(v[x, y, z])) + nu_den[x-1, y, z-1] * dot(-1.0, 0.0, -1.0, norm(v[x, y, z])) + nu_den[x-1, y, z+1] * dot(-1.0, 0.0, 1.0, norm(v[x, y, z])) + nu_den[x+1, y-1, z] * dot(1.0, -1.0, 0.0, norm(v[x, y, z])) + nu_den[x+1, y+1, z] * dot(1.0, 1.0, 0.0, norm(v[x, y, z])) + nu_den[x+1, y, z-1] * dot(1.0, 0.0, -1.0, norm(v[x, y, z])) + nu_den[x+1, y, z+1] * dot(1.0, 0.0, 1.0, norm(v[x, y, z])) + nu_den[x, y-1, z-1] * dot(0.0, -1.0, -1.0, norm(v[x, y, z])) + nu_den[x, y+1, z-1] * dot(0.0, 1.0, -1.0, norm(v[x, y, z])) + nu_den[x, y-1, z+1] * dot(0.0, -1.0, 1.0, norm(v[x, y, z])) + nu_den[x, y+1, z+1] * dot(0.0, 1.0, 1.0, norm(v[x, y, z]))) / (12.0 * np.sqrt(2))
                    corr3 += dm_den[x, y, z] * (nu_den[x-1, y-1, z-1] * dot(-1.0, -1.0, -1.0, norm(v[x, y, z])) + nu_den[x-1, y-1, z+1] * dot(-1.0, -1.0, 1.0, norm(v[x, y, z])) + nu_den[x-1, y+1, z-1] * dot(-1.0, 1.0, -1.0, norm(v[x, y, z])) + nu_den[x-1, y+1, z+1] * dot(-1.0, 1.0, 1.0, norm(v[x, y, z])) + nu_den[x+1, y-1, z-1] * dot(1.0, -1.0, -1.0, norm(v[x, y, z])) + nu_den[x+1, y-1, z+1] * dot(1.0, -1.0, 1.0, norm(v[x, y, z])) + nu_den[x+1, y+1, z-1] * dot(1.0, 1.0, -1.0, norm(v[x, y, z])) + nu_den[x+1, y+1, z+1] * dot(1.0, 1.0, 1.0, norm(v[x, y, z]))) / (8.0 * np.sqrt(3))

        corr.append(corr1/(Nden-2)**3)
        corr.append(corr2/(Nden-2)**3)
        corr.append(corr3/(Nden-2)**3)

    return r, corr


