from pymol import cmd
import chempy
import numpy as np
import LabelLib as ll

# Usage example:
# fetch 1BNA, async=0
# remove solvent 
# genAV('1BNA', '/1BNA/B/B/19/C5', allowed_sphere_radius=1.5)
## genAV('1BNA', '/1BNA/B/B/19/C5', linker_length=22.0, linker_diameter=3.0, dye_radius=4.0, disc_step=0.7, allowed_sphere_radius=2.0)
def genAV(obstacles, attachment, linker_length=20.0, linker_diameter=2.0, dye_radius=3.5, disc_step=0.9, name=None, state=1, stripsc=True, allowed_sphere_radius=0.0, smoothSurf=True):
  
  source = np.array(cmd.get_model(attachment, state).get_coord_list())
  if (source.shape[0]!=1):
    print('attachment selection must contain exactly one atom, selected: {}'.format(source.shape[0]))
    return
  source=source.reshape(3)
  
  srcAt = cmd.get_model(attachment, state).atom[0]
  srcModelName = cmd.get_names('objects',0,attachment)[0]
  
  obstacles = '(' + obstacles + ') and not (' + attachment + ')'
  if stripsc and isAA(srcAt.resn):
    obstacles += ' and not (' + srcModelName
    if len(srcAt.chain)>0:
      obstacles += ' and chain ' + srcAt.chain
    obstacles+=' and resi '+srcAt.resi+' and sidechain'+')'
  if allowed_sphere_radius > 0.0:
    obstacles+=' and not (({}) around {})'.format(attachment, allowed_sphere_radius)
  
  xyzRT=np.zeros((1,4))
  nAtoms=cmd.count_atoms(obstacles)
  if nAtoms>0:
    atoms=cmd.get_model(obstacles, state).atom
    nAtoms=len(atoms)
    xyzRT=np.zeros((nAtoms,4))
    for i,at in enumerate(atoms):
      xyzRT[i]=[at.coord[0],at.coord[1],at.coord[2],at.vdw]
  
  av1=ll.dyeDensityAV1(xyzRT.T,source,linker_length, linker_diameter, dye_radius, disc_step)
  m=avToModel(av1)
  if len(m.atom)==0:
    print('Failed: Empty AV. Is attachment position buried?')
    return
  if name is None:
    name = srcModelName + '_'
    if len(srcAt.chain)>0:
      name += srcAt.chain + '-'
    name +=  srcAt.resi + '-' + srcAt.name
  cmd.load_model(m, name)
  
  if smoothSurf:
    surfName=name+'_surf'
    mapName=name+'_map'
    gRes=cmd.get('gaussian_resolution')
    cmd.set('gaussian_resolution',3.0)
    cmd.map_new(mapName,'gaussian', 1.0, name, 6)
    cmd.isosurface(surfName,mapName,0.9)
    cmd.set('gaussian_resolution',gRes)
    cmd.disable(name)

def isAA(resn):
  names=['ALA', 'ARG', 'ASN', 'ASP', 'ASX',
	 'CYS', 'GLU', 'GLN', 'GLX', 'GLY',
	 'HIS', 'ILE', 'LEU', 'LYS', 'MET',
	 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
  if resn in names:
    return True
  return False

def makeAtom(index, xyz, vdw, name='AV'):
    atom = chempy.Atom()
    atom.index = index
    atom.name = name
    atom.symbol = 'He'
    atom.resn = 'AV'
    atom.chain = 'A'
    atom.resi = 1
    atom.resi_number = 1
    atom.coord = xyz
    atom.vdw=vdw
    atom.hetatm = False
    atom.b = 100
    atom.q = 1
    return atom

def avToModel(av):
  m = chempy.models.Indexed()
  
  area = av.shape[0] * av.shape[1]
  nx, ny, nz = av.shape
  ox, oy, oz = av.originXYZ
  dx = av.discStep
  g = np.array(av.grid).reshape((nx, ny, nz),order='F')
  
  MP=np.array([0.0,0.0,0.0])
  vol=0.0
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
        
        MP+=np.array([x,y,z])*val
        vol+=val
        
        m.add_atom(makeAtom(iat,[x,y,z],dx * 0.5))
  
  MP=MP/vol
  if iat>0:
    m.add_atom(makeAtom(iat+1,list(MP),2.0,'AVmp'))
  
  m.update_index()
  return m


cmd.extend("genAV", genAV)
 
