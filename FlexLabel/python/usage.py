import LabelLib as ll
import numpy as np

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
    
  
atoms=np.array([
[0.0, -4.0, 22.0, 1.5],
[9.0, 0.0, 0.0, 3.0],
[9.0, 8.0, 0.0, 3.0],
[0.0, -10.5, 0.0, 1.5],
[5.0, 0.0, 0.0, 2.0],
[0.0, -4.0, -10.5, 1.3],
[0.0, -4.0, -11.5, 1.0],
[0.0, -4.0, -12.5, 2.5],
[0.0, -4.0, -13.5, 1.7],
[0.0, -4.0, -14.5, 1.8],
[0.0, -4.0, 95.0, 1.5]]).astype(float).T
source=np.array([0.0, -4.0, 0.0]).astype(float)

av1 = ll.dyeDensityAV1(atoms, source, 20.0, 2.0, 3.5, 0.9)
minLengthGrid = ll.minLinkerLength(atoms, source, 20.0, 2.0, 3.5, 0.9)

print('Saving AVs...')
savePqrFromAtoms('atoms.pqr', atoms)
savePqr('AV1.pqr', av1)
savePqr('minLinkerLength.pqr', minLengthGrid)

Emean = ll.meanEfficiency(av1,av1,52.0,100000)
print('Emean = {}'.format(Emean))

#Contact volume (re)weigting
labels=np.full([1,surfaceAtoms.shape[1]],10.123) # density close to surfaceAtoms will be 10.123 units higher
surfaceAtoms=np.vstack([atoms,labels])
surfaceAtoms[3]+=2.34 #contact radius is larger than vdW radius
acv = ll.addWeights(av1,surfaceAtoms)
savePqr('ACV.pqr', acv)
print('done.')

# Grid3D initialization
g = ll.Grid3D(av1.shape, av1.originXYZ, av1.discStep)
