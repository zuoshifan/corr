#!/usr/bin/env python

from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

if size != 8:
    raise Exception("Test needs %d processes." % size)

data_path = "/scratch/p/pen/emberson/cubep3m/neutrinos/cubep3m_6144_500.0_0.2nu"
out_path = "/scratch/p/pen/zhm"

Np = 8 # process grid in one dimension
n = 192 # data grid size in one file
N = 768 # data grid in one dimension

orig = [(0, 0, 0), (4, 0, 0), (0, 4, 0), (4, 4, 0), (0, 0, 4), (4, 0, 4), (0, 4, 4), (4, 4, 4)]

den = np.zeros((N, N, N), dtype=np.float32, order='C')

for z in range(4):
    for y in range(4):
        for x in range(4):
            num = (x + orig[rank][0]) + Np * (y + orig[rank][1]) + Np**2 * (z + orig[rank][2]) # file number
            den[n*x:n*(x+1), n*y:n*(y+1), n*z:n*(z+1)] = np.fromfile((data_path+'/node%d/0.000den%d_h.bin') % (num, num), dtype=np.float32).reshape(n, n, n).T # convert to C array

den = den.T # convert to fortran array
den.tofile((out_path + '/0.000den%d_h.bin') % rank)
