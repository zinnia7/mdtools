#!/usr/bin/python3
# LAMMPS units should be "real"
# position data should be unwrapped coordinates (xu, yu, zu)

import sys
from ase.io import read
from ase.calculators.lammps import convert
import numpy as np

# original msd method
def msd_original(r):
  shifts = np.arange(len(r))
  msds = np.zeros((shifts.size,3))
  for i, shift in enumerate(shifts):
    diffs = r[:-shift if shift else None] - r[shift:]
    sqdist = np.square(diffs)
    msds[i,0] = sqdist[:,0].mean()
    msds[i,1] = sqdist[:,1].mean()
    msds[i,2] = sqdist[:,2].mean()
  return msds

# msd by FFT
def msd_fft(r):
  n = len(r)
  s1 = np.zeros((n,3))
  s2 = np.zeros((n,3))
  msds = np.zeros((n,3))
  for i in range(3):
    fk = np.fft.fft(r[:,i], n=2*n)
    power = fk * fk.conjugate()
    res = np.fft.ifft(power)[:n].real
    s2[:,i] = res/(n*np.ones(n)-np.arange(0,n))
    x2 = r[:,i]**2
    s1[0,i] = np.average(x2)*2.0
    for m in range(1,n):
      s1[m,i] = np.average(x2[m:] + x2[:-m])
    msds[:,i] = s1[:,i] - 2*s2[:,i]
  return msds

if len(sys.argv) != 4:
  print("# usage: msd.py lammps_trajectory_file target_atom_type save_timestep(fs)")
  exit()
else:
  lammpstrj = sys.argv[1]
  target_type = int(sys.argv[2])
  save_timestep = float(sys.argv[3])

# read lammpstrj as ase format
trjs = read(lammpstrj, format="lammps-dump-text", index=":", units="real")

# summary of lammpstrj
natoms = trjs[0].get_global_number_of_atoms()
ntypes = len(set(trjs[0].get_chemical_symbols()))
ntrjs = len(trjs)

print("# number of atoms: {}".format(natoms))
print("# number of atom types: {}".format(ntypes))
print("# number of snapshots: {}".format(ntrjs))
print("# target atom type: {}".format(target_type))
print("# time interval between snapshots: {}".format(save_timestep))

msds = np.zeros((ntrjs,3))
ntarget = 0
for i in range(0,natoms):
  rus = np.empty((0,3))
  if trjs[0].get_atomic_numbers()[i] == target_type:
    ntarget += 1
    for j in range(0, ntrjs):
      rus = np.append(rus, trjs[j].positions[i].reshape(1,-1),axis=0)
    rus = convert(rus, "distance", "ASE", "real")
    msds += msd_fft(rus)
#    msds += msd_original(rus)
msds /= ntarget

print("# number of target atoms: {}".format(ntarget))

f = open("msd.dat","w")
f.write("# time(fs) msdx msdy msdz msd(A^2)\n")
for i in range(ntrjs):
  f.write("{0:.1f} {1:.5f} {2:.5f} {3:.5f} {4:.5f}\n".format(i*save_timestep,msds[i,0],msds[i,1],msds[i,2],msds[i,0]+msds[i,1]+msds[i,2]))
f.close()

print("# msd.dat file has been saved")
