#!/usr/bin/python3

import sys
import numpy as np
import pandas as pd

filename = sys.argv[1]
name = filename.split(".")[0]

with open(filename) as file:
  lines = [[word for word in line.split()] for line in file]

tag = ""
masses = []
atoms = []
bonds = []
angles = []
dihedrals = []
impropers = []
pair_coeffs = []
bond_coeffs = []
angle_coeffs = []
dihedral_coeffs = []
improper_coeffs = []

for n, line in enumerate(lines):
  if len(line) == 2 and line[1] == "atoms":
    natoms = int(line[0])
  if len(line) == 2 and line[1] == "bonds":
    nbonds = int(line[0])
  if len(line) == 2 and line[1] == "angles":
    nangles = int(line[0])
  if len(line) == 2 and line[1] == "dihedrals":
    ndihedrals = int(line[0])
  if len(line) == 2 and line[1] == "impropers":
    nimpropers = int(line[0])
  if len(line) == 3 and line[1] == "atom":
    ntatoms = int(line[0])
  if len(line) == 3 and line[1] == "bond":
    ntbonds = int(line[0])
  if len(line) == 3 and line[1] == "angle":
    ntangles = int(line[0])
  if len(line) == 3 and line[1] == "dihedral":
    ntdihedrals = int(line[0])
  if len(line) == 3 and line[1] == "improper":
    ntimpropers = int(line[0])

  if len(line) != 0 and (line[0] == "Masses" or line[0] == "Pair" or line[0] == "Atoms" \
     or line[0] == "Bond" or line[0] == "Angle" or line[0] == "Dihedral" or line[0] == "Improper" \
     or line[0] == "Bonds" or line[0] == "Angles" or line[0] == "Dihedrals" or line[0] == "Impropers"):
    tag = line[0]
  elif len(line) !=0 and tag != "":
    if tag == "Masses":
      masses.append(line)
    if tag == "Pair":
      pair_coeffs.append(line)
    if tag == "Atoms":
      atoms.append(line)
    if tag == "Bonds":
      bonds.append(line)
    if tag == "Angles":
      angles.append(line)
    if tag == "Dihedrals":
      dihedrals.append(line)
    if tag == "Impropers":
      impropers.append(line)
    if tag == "Bond":
      bond_coeffs.append(line)
    if tag == "Angle":
      angle_coeffs.append(line)
    if tag == "Dihedral":
      dihedral_coeffs.append(line)
    if tag == "Improper":
      improper_coeffs.append(line)

paircoeffs = np.array(pair_coeffs, dtype="float64")
u_pc, u_pid, u_prv = np.unique(paircoeffs[:,1:3], axis=0, return_index=True, return_inverse=True)
newpcs = paircoeffs[np.sort(u_pid)]
newpids = []
for i in u_prv:
  for ar in paircoeffs:
    if (ar[1:3] == u_pc[i]).all():
      newpids.append(str(int(ar[0])))
      break

bondcoeffs = np.array(bond_coeffs, dtype="float64")
u_bc, u_bid, u_brv = np.unique(bondcoeffs[:,1:3], axis=0, return_index=True, return_inverse=True)
newbcs = bondcoeffs[np.sort(u_bid)]
newbids = []
for i in u_brv:
  for ar in bondcoeffs:
    if (ar[1:3] == u_bc[i]).all():
      newbids.append(str(int(ar[0])))
      break

anglecoeffs = np.array(angle_coeffs, dtype="float64")
u_ac, u_aid, u_arv = np.unique(anglecoeffs[:,1:3], axis=0, return_index=True, return_inverse=True)
newacs = anglecoeffs[np.sort(u_aid)]
newaids = []
for i in u_arv:
  for ar in anglecoeffs:
    if (ar[1:3] == u_ac[i]).all():
      newaids.append(str(int(ar[0])))
      break

dihedralcoeffs = np.array(dihedral_coeffs, dtype="float64")
u_dc, u_did, u_drv = np.unique(dihedralcoeffs[:,1:5], axis=0, return_index=True, return_inverse=True)
newdcs = dihedralcoeffs[np.sort(u_did)]
newdids = []
for i in u_drv:
  for ar in dihedralcoeffs:
    if (ar[1:5] == u_dc[i]).all():
      newdids.append(str(int(ar[0])))
      break

impropercoeffs = np.array(improper_coeffs, dtype="float64")
u_ic, u_iid, u_irv = np.unique(impropercoeffs[:,1:5], axis=0, return_index=True, return_inverse=True)
newics = impropercoeffs[np.sort(u_iid)]
newiids = []
for i in u_irv:
  for ar in impropercoeffs:
    if (ar[1:5] == u_ic[i]).all():
      newiids.append(str(int(ar[0])))
      break


f = open(name + ".lt", "w")
f.write(name + " {\n\n")

f.write('write("Data Atoms") {\n')
for i in range(natoms):
  f.write('{0:15s}{1:10s}{2:15s}{3:10.4f}{4:15.10f}{5:15.10f}{6:15.10f}\n'.format("$atom:A"+atoms[i][0],"$mol:M"+atoms[i][1],"@atom:A"+newpids[i],float(atoms[i][3]),float(atoms[i][4]),float(atoms[i][5]),float(atoms[i][6])))
f.write('}\n\n')

if ntbonds != 0:
  f.write('write("Data Bonds") {\n')
  for i in range(nbonds):
    f.write('{0:20s}{1:20s}{2:15s}{3:15s}\n'.format("$bond:B"+bonds[i][0],"@bond:B"+newbids[i],"$atom:A"+bonds[i][2],"$atom:A"+bonds[i][3]))
  f.write('}\n\n')

if ntangles != 0:
  f.write('write("Data Angles") {\n')
  for i in range(nangles):
    f.write('{0:20s}{1:20s}{2:15s}{3:15s}{4:15s}\n'.format("$angle:A"+angles[i][0],"@angle:A"+newaids[i],"$atom:A"+angles[i][2],"$atom:A"+angles[i][3],"$atom:A"+angles[i][4]))
  f.write('}\n\n')

if ntdihedrals != 0:
  f.write('write("Data Dihedrals") {\n')
  for i in range(ndihedrals):
    f.write('{0:20s}{1:20s}{2:15s}{3:15s}{4:15s}{5:15s}\n'.format("$dihedral:D"+dihedrals[i][0],"@dihedral:D"+newdids[i],"$atom:A"+dihedrals[i][2],"$atom:A"+dihedrals[i][3],"$atom:A"+dihedrals[i][4],"$atom:A"+dihedrals[i][5]))
  f.write('}\n\n')

if ntimpropers != 0:
  f.write('write("Data Impropers") {\n')
  for i in range(nimpropers):
    f.write('{0:20s}{1:20s}{2:15s}{3:15s}{4:15s}{5:15s}\n'.format("$improper:I"+impropers[i][0],"@improper:I"+newiids[i],"$atom:A"+impropers[i][2],"$atom:A"+impropers[i][3],"$atom:A"+impropers[i][4],"$atom:A"+impropers[i][5]))
  f.write('}\n\n')

f.write('write_once("Data Masses") {\n')
for i in range(len(u_pc)):
  f.write('{0:15s}{1:10.4f}\n'.format("@atom:A"+masses[int(newpcs[i,0])-1][0],float(masses[int(newpcs[i,0])-1][1])))
f.write('}\n\n')

f.write('write_once("In Settings") {\n')
for i in range(len(newpcs)):
  f.write('{0:15s}{1:15s}{1:15s}{2:15.5f}{3:15.5f}\n'.format("pair_coeff","@atom:A"+str(int(newpcs[i,0])),newpcs[i,1],newpcs[i,2]))
f.write('\n')
for i in range(len(newbcs)):
  f.write('{0:15s}{1:15s}{2:15.5f}{3:15.5f}\n'.format("bond_coeff","@bond:B"+str(int(newbcs[i,0])),newbcs[i,1],newbcs[i,2]))
f.write('\n')
for i in range(len(newacs)):
  f.write('{0:15s}{1:15s}{2:15.5f}{3:15.5f}\n'.format("angle_coeff","@angle:A"+str(int(newacs[i,0])),newacs[i,1],newacs[i,2]))
f.write('\n')
for i in range(len(newdcs)):
  f.write('{0:15s}{1:15s}{2:10.4f}{3:10.4f}{4:10.4f}{5:10.4f}\n'.format("dihedral_coeff","@dihedral:D"+str(int(newdcs[i,0])),newdcs[i,1],newdcs[i,2],newdcs[i,3],newdcs[i,4]))
f.write('\n')
for i in range(len(newics)):
  f.write('{0:15s}{1:15s}{2:10.4f}{3:5d}{4:5d}\n'.format("improper_coeff","@improper:I"+str(int(newics[i,0])),newics[i,1],int(newics[i,2]),int(newics[i,3])))
f.write('}\n\n')

f.write('write_once("In Init") {\n')
f.write('units           real\n')
f.write('atom_style      full\n')
f.write('pair_style      lj/charmm/coul/long 9.0 10.0\n')
if ntbonds != 0:
  f.write('bond_style      harmonic\n')
if ntangles != 0:
  f.write('angle_style     harmonic\n')
if ntdihedrals != 0:
  f.write('dihedral_style  opls\n')
if ntimpropers != 0:
  f.write('improper_style  cvff\n')
f.write('kspace_style    pppm 0.0001\n')
f.write('pair_modify     mix geometric\n')
f.write('special_bonds   lj/coul 0.0 0.0 0.5\n')
f.write('}\n\n')

f.write('}')
f.close()
