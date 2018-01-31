from pymol import cmd
import chempy
import numpy as np
import LabelLib as ll

# Usage example: 
# genAV('chain A', 'chain A and resi 123 and name CB')
def genAV(obstacles, attachment, linker_length=20.0, linker_diameter=2.0, dye_radius=3.5, disc_step=0.9, name=None, state=1, stripAttSC=True):
  
  source = np.array(cmd.get_model(attachment, state).get_coord_list())
  if (source.shape[0]!=1):
    print('attachment selection must contain exactly one atom, selected: {}'.format(source.shape[0]))
    return
  source=source.reshape(3)
  
  srcAt = cmd.get_model(attachment, state).atom[0]
  srcModelName = cmd.get_names('objects',0,attachment)[0]
  
  obstacles = '(' + obstacles + ') and not (' + attachment + ')'
  if stripAttSC:
    obstacles += ' and not (' + srcModelName
    if len(srcAt.chain)>0:
      obstacles += ' and chain ' + srcAt.chain
    obstacles+=' and resi '+srcAt.resi+' and sidechain'+')'
  
  atoms=cmd.get_model(obstacles, state).atom
  nAtoms=len(atoms)
  xyzRT=np.zeros((nAtoms,4))
  for i,at in enumerate(atoms):
    xyzRT[i]=[at.coord[0],at.coord[1],at.coord[2],at.vdw]
  
  av1=ll.dyeDensityAV1(xyzRT.T,source,linker_length, linker_diameter, dye_radius, disc_step)
  m=avToModel(av1)
  if name is None:
    name = srcModelName + '_'
    if len(srcAt.chain)>0:
      name += srcAt.chain + '-'
    name +=  srcAt.resi + '-' + srcAt.name
  cmd.load_model(m, name)

def makeAtom(index, xyz, vdw):
    atom = chempy.Atom()
    atom.index = index
    atom.name = 'AV'
    atom.symbol = 'AV'
    atom.resn = 'AV'
    atom.chain = 'A'
    atom.resi = 1
    atom.resi_number = 1
    atom.coord = xyz
    atom.vdw=vdw
    atom.hetatm = False
    return atom

def avToModel(av):
  m = chempy.models.Indexed()
  
  area = av.shape[0] * av.shape[1]
  nx, ny, nz = av.shape
  ox, oy, oz = av.originXYZ
  dx = av.discStep
  g = np.array(av.grid).reshape((nx, ny, nz),order='F')
  
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
        
        m.add_atom(makeAtom(iat,[x,y,z],dx * 0.5))
  
  m.update_index()
  if iat==0:
    print('Empty AV. Is attachment position buried?')
  return m


cmd.extend("genAV", genAV)
 
