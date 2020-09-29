import os
os.chdir('/home/dongdong/SCR/chignolin/anton/CLN025-0-protein/')
import mdtraj as md
pdb = md.load_pdb('crystal.pdb')
topology = pdb.topology
CA=topology.select('name CA')
for i in range(54):
    traj = md.load_dcd('CLN025-0-protein-%03d.dcd' %i, top='CLN025-0-protein.pdb')
    rmsds = md.rmsd(traj, pdb, 0,CA)
    np.savetxt('rmsd%03d.dat' %i, rmsds)
    print(np.min(rmsds))
