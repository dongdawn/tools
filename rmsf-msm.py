#!/usr/bin/env python
import os
import sys
import argparse

import numpy as np
import mdtraj as md
import pandas as pd

from rmsf_xtc import cal_rmsf
def parse_args():
    parser = argparse.ArgumentParser(description="Calculate the rmsf for the msm ensemble.")
    parser.add_argument("-m", dest="msmf", help="The msm pickle file")
    parser.add_argument("-s", dest="samf", help="The msm samples floder.")
    parser.add_argument("--sel", dest="sel", help="Select atom to cal rmsf. default: protein and name CA", default="protein and name CA")
    parser.add_argument("-t", dest="topf", help="The top file for the sample file.")
    parser.add_argument("-o", dest="outf", help="The result file to save.")

    args = parser.parse_args()
    return args.msmf, args.samf, args.topf, args.outf, args.sel

def main():
    msmf, samf, topf, outf, sel = parse_args()
    M = pd.read_pickle(msmf)
    msm_rmsf = None

    for c_i in xrange(M.n_states_):
        xtcf = os.path.join(samf, "%d.xtc"%c_i)
        popu = M.populations_[c_i]

        rmsf_xtc, atom_index = cal_rmsf(xtcf, topf, atom_sel=sel)
        if msm_rmsf is None:
            msm_rmsf = rmsf_xtc * popu
        else:
            msm_rmsf += rmsf_xtc * popu

    rs = np.array([atom_index, msm_rmsf]).T
    np.savetxt(outf, rs, fmt="%d %.4f")

if __name__ == "__main__":
    main()
