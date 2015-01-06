#!/usr/bin/env python

"""Visualize correlation function.

:Authors: Shifan Zuo
:Date: 2014-12-29
:email: sfzuo@bao.ac.cn
"""

import argparse


def visualize_corr(args):
    """Visualize correlation function in hdf5 files.

    Arguments
    ---------
    args : argparse namespace.
    """
    import numpy as np
    import h5py
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt

    # Read in correlation data
    with h5py.File(args.filename, 'r') as f:
        r = f['r'][:]
        corr = f['corr'][:]

    # Create output image file name
    if args.outfile:
        out_file = args.outfile
    else:
        out_file = 'corr.' + args.figfmt

    # Plot and save image
    plt.figure(figsize=(args.figlength,args.figwidth))
    plt.plot(r, corr, 'o')
    # plt.plot(r, np.log10(corr), 'o')
    # plt.plot(np.log10(r), np.log10(corr), 'o')
    plt.savefig(out_file)
    plt.close()


parser = argparse.ArgumentParser(description='Visualize a sky map in hdf5 files.')
parser.add_argument('filename', type=str, nargs='?', default='corr.hdf5', help='Input hdf5 correlation.')
parser.add_argument('-o', '--outfile', type=str, nargs='?', help='Name of the image file to save into.')
parser.add_argument('-f', '--figfmt', default='png', help='Output image format.')
# parser.add_argument('--min', type=float, help='The min value of the visualize range in the output image.')
# parser.add_argument('--max', type=float, help='The max value of the visualize range in the output image.')
parser.add_argument('-l', '--figlength', type=float, default=8, help='Output figure length.')
parser.add_argument('-w', '--figwidth', type=float, default=6, help='Output figure width.')
# parser.add_argument('-t', '--tight', action='store_true', help='Tight the figure marin space.')
parser.set_defaults(func=visualize_corr)

args = parser.parse_args()
args.func(args)
