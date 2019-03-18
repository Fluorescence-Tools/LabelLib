#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import LabelLib as ll
import numpy as np
import mdtraj as md
import warnings

# Prints inter-dye distance distribution histograms for each frame of a trajectory.
# Labels are: donor at residue 2, acceptor at residue 27 
def main(args):
	traj = md.load_pdb('http://www.rcsb.org/pdb/files/2EQQ.pdb')
	
	#strip sidechains of residues 2 and 27
	idxs=traj.topology.select('name CB CA N C or not (resSeq 2 or resSeq 27)')
	traj=traj.atom_slice(idxs)
	
	bins=np.linspace(0.0,80.0,41)
	bin_centers=bins[:-1]+(bins[1]-bins[0])/2.0
	
	#get the index of donor and acceptor attachment atoms
	donorAttIdx=traj.topology.select('name CB and resSeq 2')[0]
	acceptorAttIdx=traj.topology.select('name CB and resSeq 27')[0]
	
	#print header
	print('Frame#\t'+'A\t'.join('{:.1f}'.format(e) for e in bin_centers))
	
	for frIdx in range(traj.n_frames):
		gridDonor=genAV(traj, frIdx, donorAttIdx, 20.0, 2.0, 3.5, 0.9, 0.3, 3.5+3.0)
		gridAcceptor=genAV(traj, frIdx, acceptorAttIdx, 22.0, 2.0, 3.5, 0.9, 0.4, 3.5+3.0)
		
		donorVolSize=np.count_nonzero(np.array(gridDonor.grid) > 1.0)
		accVolSize=np.count_nonzero(np.array(gridAcceptor.grid) > 1.0)
		if donorVolSize==0 or accVolSize==0:	
			continue
		
		distances=ll.sampleDistanceDistInv(gridDonor,gridAcceptor,1000000)
		freq, bin_edges = np.histogram(distances,bins=bins)
		print('{}\t'.format(frIdx)+'\t'.join(str(e) for e in freq))
			

# Converts mdtraj frame to xyzR array for LabelLib
def xyzr(traj, frameId):
	vdw=np.empty(traj.n_atoms)
	for i,atom in enumerate(traj.topology.atoms):
		vdw[i]=atom.element.radius
	# For better performance `vdw` should be cached
		
	return np.column_stack([traj.xyz[frameId],vdw]).astype(np.float32).T*10.0

#cvRadius does not include the dyeRadius, so one would typically set cvRadius=dyeR+const
def genAV(traj, frameId, attachmentIdx, length, width, radius, resolution, cvFrac=None, cvRadius=0.0):
	
	atoms=xyzr(traj,frameId)
	atoms[3,attachmentIdx]=0.0 #attachment atom should not be an obstacle
	source=atoms[:3,attachmentIdx]
	av=ll.dyeDensityAV1(atoms, source, length, width, radius, resolution)
	if cvFrac is None:
		return av
	#ACV
	labels=np.full([1,atoms.shape[1]], 2.0) #density >= 2 means "contact volume"
	surfaceAtoms=np.vstack([atoms,labels])
	surfaceAtoms[3]+=cvRadius
	acv = ll.addWeights(av,surfaceAtoms)
	
	#reweight
	acvarr=np.array(acv.grid)
	volCV = np.count_nonzero(acvarr > 1.0)
	volFree = np.count_nonzero(acvarr == 1.0)
	newWeight = 1.0
	if volCV==0.0:
		warnings.warn('Contact volume is empty, please check the settings.')
		return acv
	newWeight = volFree*cvFrac/(volCV*(1.0-cvFrac))
	acvarr[acvarr>1.0]=newWeight
	acv.grid=list(acvarr)
	return acv

def savePqr(fileName, grid):
  with open(fileName, "w") as out:
      area = grid.shape[0] * grid.shape[1]
      nx, ny, nz = grid.shape
      ox, oy, oz = grid.originXYZ
      dx = grid.discStep
      g = np.array(grid.grid).reshape((nx, ny, nz),order='F')
      
      iat = 0
      for iz in range(nz):
        for iy in range(ny):
          for ix in range(nx):
            val = g[ix, iy, iz]
            if val <= 0.0:
              continue
            
            iat += 1
            resi = int(iat / 10)
            
            x = ix * dx + ox
            y = iy * dx + oy
            z = iz * dx + oz
            
            sz = 'ATOM{0: 7}   AV  AV{1: 6}{2:12.1f}{3:8.1f}{4:8.1f}{5:8.2f}{6:7.3f}\n'
            sz = sz.format(iat, resi, x, y, z, val, dx * 0.5)
            out.write(sz)

def savePqrFromAtoms(fileName, atoms):
  sz = 'ATOM{0: 7}    X   X{1: 6}{2:12.1f}{3:8.1f}{4:8.1f}{5:8.2f}{6:7.3f}\n'
  with open(fileName, "w") as out:
    for i,at in enumerate(atoms.T):
      out.write(sz.format(i, i, at[0], at[1], at[2], 1.0, at[3]))

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
