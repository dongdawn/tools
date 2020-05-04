#!/usr/bin/env python
import os
import sys
import argparse

import numpy as np
import mdtraj as md
import pandas as pd

def cal_rmsf(trajf, topf, atom_sel="protein and name CA"):
    traj = md.load(trajf, top=topf)
    ref  = md.load(topf)
    atom_index = traj.top.select(atom_sel)
    traj = traj.superpose(ref, 0)
    traj = traj.atom_slice(atom_index)

    sumsquares = np.zeros((traj.top.n_atoms, 3))
    xyz_mean   = sumsquares.copy()
    n_frames   = traj.n_frames

    for i, frame in enumerate(traj):
        xyz = frame.xyz[0]
        sumsquares += (i / (i+1.0)) * (xyz - xyz_mean) ** 2
        xyz_mean = (i * xyz_mean + xyz) / (i + 1)

    rmsf = np.sqrt(sumsquares.sum(axis=1) / (n_frames))

    '''
    xyz = traj.xyz[:, :, :]
    xyz_mean = xyz.mean(axis=0)

    n_frame = traj.n_frames

    sumsquares = np.sum( n_frame /float(n_frame+1) * (xyz - xyz_mean)**2, axis=0)
    rmsf = np.sqrt(sumsquares.sum(axis=1)/(n_frame+1))
   '''
    return rmsf,atom_index

def parse_args():
    parser = argparse.ArgumentParser(description="Calculate the rmsf of trajectory file.")
    parser.add_argument("-x", dest="trajf", help="The xtc file to cal rmsf.", required=True)
    parser.add_argument("-t", dest="topf", help="The top file of the trjactory file.", required=True)
    parser.add_argument("-s", dest="sel", help="The selection to calculate rmsf. default: C-alpha", default="protein and name CA")
    parser.add_argument("-o", dest="outf", help="The result file path to save. default=rmsf-xtc.txt", default="rmsf-xtc.txt")

    args = parser.parse_args()

    return args.trajf, args.topf, args.sel, args.outf

def main():
    trajf, topf, sel, outf = parse_args()
    rmsf, atom_index = cal_rmsf(trajf, topf, atom_sel=sel)
    rs = np.array([atom_index, rmsf]).T
    print rs.shape
    print rs[0].shape, rs[1].shape
    np.savetxt(outf, rs, fmt="%d %.4f")

if __name__ == "__main__":
    main()
