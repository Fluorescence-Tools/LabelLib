import LabelLib as ll
import numpy as np

def savePqr(fileName, grid):
  out = open(fileName, "w")

  area=grid.shape[0]*grid.shape[1]
  for idx in range(len(grid.grid)):
    val = grid.grid[idx]
    if val<=0.0:
      continue
    
    k = int(idx / area)
    tmp = idx - k * area
    j = int(tmp / grid.shape[0])
    i = int(tmp % grid.shape[0])
    
    x = i * grid.discStep + grid.originXYZ[0]
    y = j * grid.discStep + grid.originXYZ[1]
    z = k * grid.discStep + grid.originXYZ[2]
    
    sz = 'ATOM{0: 7}   AV  AV{1: 6}{2:12.1f}{3:8.1f}{4:8.1f}{5:8.2f}{6:7.3f}\n'
    sz = sz.format(idx, idx, x, y, z, val, grid.discStep * 0.5)
    out.write(sz)

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
source=np.array([[0.0, -4.0, 0.0],]).astype(float).T

av1 = ll.dyeDensityAV1(atoms, source, 20.0, 2.0, 3.5, 0.9)
minLengthGrid = ll.minLinkerLength(atoms, source, 20.0, 2.0, 3.5, 0.9)

print('Saving AVs. This can take ~10 min...')
#savePqr('AV1.pqr', av1)
savePqr('minLinkerLength.pqr', minLengthGrid)
