#!/usr/bin/python3
# unit of velocity should be A/fs (real unit)

import sys
from ase.io import read
from ase.calculators.lammps import convert
import numpy as np

# vac by FFT
def vac_fft(v):
  n = len(v)
  vac = np.zeros((n,3))
  for i in range(3):
    fk = np.fft.fft(v[:,i], n=2*n)
    power = fk * fk.conjugate()
    vac[:,i] = np.fft.ifft(power)[:n].real/n
  return vac

if len(sys.argv) != 4:
  print("# usage: vac.py lammps_trajectory_file target_atom_type save_timestep(fs)")
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

vac = np.zeros((ntrjs,3))
ntarget = 0

for i in range(0,natoms):
  vs = np.empty((0,3))
  if trjs[0].get_atomic_numbers()[i] == target_type:
    ntarget += 1
    for j in range(0, ntrjs):
      vs = np.append(vs, trjs[j].get_velocities()[i].reshape(1,-1),axis=0)
    vs = convert(vs, "velocity", "ASE", "real")
    vac += vac_fft(vs)
vac /= ntarget

print("# number of target atoms: {}".format(ntarget))

f = open("vac.dat","w")
f.write("# time(fs) vacx vacy vacz vac(A^2/fs^2) ivac(A^2/fs)\n")
ivac = 0 # for integral
for i in range(ntrjs):
  vac_tot = vac[i,0] + vac[i,1] + vac[i,2]
  ivac += vac_tot
  f.write("{:.1f} {:.5e} {:.5e} {:.5e} {:.5e} {:.5e}\n".format(i*save_timestep,vac[i,0],vac[i,1],vac[i,2],vac_tot,ivac*save_timestep))
f.close()

print("# vac.dat file has been saved")
