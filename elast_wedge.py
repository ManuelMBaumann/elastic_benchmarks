from __future__ import print_function, division

import matplotlib
matplotlib.use('agg')

from nutils import *
import scipy
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.io import mmwrite

@log.title
def makeplots( domain, geom, Lx, Lz, value, name, title, ndigits=0, index=None, clim=None, lineOn=False, imgtype=None,):
  points, colors = domain.elem_eval( [ geom, value ], ischeme='bezier3', separate=True )

  with plot.PyPlot( name, ndigits=ndigits, figsize=(5,6), index=index, imgtype=imgtype ) as plt:
    plt.mesh( points, colors, triangulate='bezier', edgecolors='none' )
    plt.title(title)
    plt.xlabel('x [m]')
    plt.ylabel('z [m]')

    plt.xticks( [0, Lx/2.0, Lx], ['0', '300', '600'] )
    plt.yticks( [-0, -0.4*Lz, -0.8*Lz, -Lz], ['0', '400', '800', '1000'] )

    if clim is not None:
      plt.clim(*clim)
    plt.colorbar()

    if lineOn:
      # Only for wedge problem in 2D
      plt.plot( [0, 600],[-400, -500],'k' )
      plt.plot( [0, 600],[-800, -600],'k' )

def makevtk(domain, geom, rho, lam, mu, cp, cs, sol, freq, vec_basis, name):
  Nom = len(freq)
  vtk_geom, vtk_rho, vtk_lam, vtk_mu, vtk_cp, vtk_cs = domain.simplex.elem_eval( [ geom, rho, lam, mu, cp, cs ], ischeme='vtk', separate=True )
  with plot.VTKFile( name ) as vtk:
      vtk.unstructuredgrid( vtk_geom )
      vtk.pointdataarray( 'rho', vtk_rho )
      vtk.pointdataarray( 'lambda', vtk_lam )
      vtk.pointdataarray( 'mu', vtk_mu )
      vtk.pointdataarray( 'cp', vtk_cp )
      vtk.pointdataarray( 'cs', vtk_cs )
      for k in range(0,Nom):
          disp = vec_basis.dot( sol[k,:] ).real
          vtk_disp = domain.simplex.elem_eval( disp, ischeme='vtk', separate=True )
          vtk.pointdataarray( 'disp_f'+str(freq[k]), vtk_disp )

def makespyplot( matrix, name, imgtype=None ):
  if not scipy.sparse.isspmatrix( matrix ):
      matrix = matrix.toscipy()
  with plot.PyPlot( name, ndigits=0, imgtype=imgtype ) as plt:
    plt.spy( matrix, markersize=0.8, color='black')
    plt.title( name+', nnz = '+str(matrix.nnz) )

def point_eval(func, domain, geom, point):
  domain = domain[tuple(slice(0, p) if p > 0 else slice(None) for p in point)]
  for p in point:
      domain = domain.boundary['right' if p > 0 else 'left']
  return numpy.asarray(domain.integrate( func, geometry=geom, ischeme='gauss2' ).toscipy().todense())

def elast_mat(rho, cp, cs, lam, mu, ndims, nx, ny, nz, vec_basis, domain, geom):
  # define PDE
  stress = lambda u: lam*u.div(geom)[:,_,_]*function.eye(ndims) + 2.0*mu*u.symgrad(geom)
  elasticity = function.outer( stress(vec_basis), vec_basis.grad(geom) ).sum([2,3])

  w_mass = lambda u: rho*u
  mass = function.outer( w_mass(vec_basis), vec_basis ).sum(-1)

  # define BC
  n = geom.normal()
  t = np.eye(ndims)
  t = t-(t*n[_,:]).sum(1)
  B_bc = cp*n[:,_]*n[_,:]+cs*(t[:,:,_]*t[:,_,:]).sum(0)

  bc_fun = lambda u: rho*(B_bc*u[:,_,:]).sum(-1)
  sommerfeld = function.outer( bc_fun(vec_basis), vec_basis ).sum(-1)

  if ndims == 2:
      sommerfeld_boundary = 'left,right,bottom'
      source_position = nx//2, nz
  else:
      sommerfeld_boundary = 'left,right,bottom,front,back'
      source_position = nx//2, ny//2, nz

  # Build matrices
  K, M = domain.integrate( [elasticity, mass], geometry=geom, ischeme='gauss2' )
  C = domain.boundary[sommerfeld_boundary].integrate( sommerfeld, geometry=geom, ischeme='gauss2' )

  # Build RHS
  rhs = point_eval(vec_basis, domain, geom, source_position)[:,-1]

  return K, C, M, rhs


def main( ndims=2,          # problem dimension (2,3) 
          dx=10.0,          # grid size in x-direction 
          dy=10.0,          # grid size in y-direction          
          dz=10.0,          # grid size in z-direction  
          freq=[32.0],      # frequency in Hz 
          plots=True,       # plot of parameters and num. solution
          spy=True,         # provide spy plot of matrices
          storing=False,    # storing=True saves matrices in matrix market format (no solve) 
          degree=1 ):       # degree of FEM splines

  # domain size
  Lx = 600.0
  Ly = 600.0
  Lz = 1000.0

  # problem parameters
  freq = np.array(freq)  
  om   = 2.0*np.pi*freq
  Nom  = len(om)

  # define physical params
  rho0 = 1800.0
  rho1 = 2100.0
  rho2 = 1950.0

  cp0  = 2000.0
  cp1  = 3000.0
  cp2  = 2300.0

  cs0  = 800.0
  cs1  = 1600.0
  cs2  = 1100.0

  # define Cartesian grid
  nx = int(np.round(Lx/dx))+1
  nz = int(np.round(Lz/dz))+1
  verts_x = np.linspace( 0, Lx, nx )
  verts_z = np.linspace( -Lz, 0, nz )

  if ndims == 2:
      ny = 1
      dy = 0.
      verts = [verts_x, verts_z]
  elif ndims == 3:
      ny = int(np.round(Ly/dy))+1
      verts_y = np.linspace( 0, Ly, ny )
      verts = [verts_x, verts_y, verts_z]

  domain, geom = mesh.rectilinear(verts)
  vec_basis    = domain.splinefunc( degree=degree ).vector( ndims )

  # define wedge problem
  rho = function.select(
      [function.greater(geom[-1]+0.4*Lz+geom[0]/6, 0), function.greater(geom[-1]+0.8*Lz-geom[0]/3, 0)],
      [rho0, rho1], rho2)
  cp  = function.select(
      [function.greater(geom[-1]+0.4*Lz+geom[0]/6, 0), function.greater(geom[-1]+0.8*Lz-geom[0]/3, 0)],
      [cp0, cp1], cp2)
  cs  = function.select(
      [function.greater(geom[-1]+0.4*Lz+geom[0]/6, 0), function.greater(geom[-1]+0.8*Lz-geom[0]/3, 0)],
      [cs0, cs1], cs2)
  
  mu   = cs**2 * rho
  lam  = rho * (cp**2 - 2.0*cs**2) 
    
  # problem summary
  print( '----     WEDGE PROBLEM     ----' )
  ppw = 20.0  
  print( 'problem size   : ' + str(nx-1+degree)+' x '+str(ny-1+degree)+' x '+str(nz-1+degree) )
  print( '# dofs         : ' + str(len(vec_basis)) )
  print( 'max. frequency : ' + str( min(cs0,cs1,cs2,cp0,cp1,cp2)/(ppw*max(dx,dy,dz)) ) )
  print( '-------------------------------\n' )

  # Create discretization matrices using nutils
  K, C, M, rhs = elast_mat(rho, cp, cs, lam, mu, ndims, nx, ny, nz, vec_basis, domain, geom,)

  if storing:
    mmwrite('matrix_io/K.mtx', K.toscipy())
    mmwrite('matrix_io/C.mtx', C.toscipy())
    mmwrite('matrix_io/M.mtx', M.toscipy())
  else:
    print('Use pythons sparse linear solver solver...')
    sol = np.zeros((Nom, len(vec_basis)), dtype=complex)
    for k in range(0,Nom):
      matrix = K + 1j*om[k]*C - om[k]**2*M
      sol[k,:] = scipy.sparse.linalg.spsolve( matrix.toscipy().tocsc(), rhs )
      if spy:
          makespyplot( matrix, 'spy_plot' )

    if plots:  
      if(ndims ==2):
          makeplots( domain, geom, Lx, Lz, rho, 'rho', 'rho [kg/m**3]' )
          #makeplots( domain, geom, Lx, Lz, cp, 'cp', 'c_p [m/s]' )
          #makeplots( domain, geom, Lx, Lz, cs, 'cs', 'c_s [m/s]' )
          
          for k in range(0,Nom):
              disp     = vec_basis.dot( sol[k,:] )          # FEM summation
              disp_x   = disp[0].real                       # Plot Re(u_x)  
              disp_z   = disp[-1].real                      # Plot Re(u_z)
              makeplots( domain, geom, Lx, Lz, disp_x, 'disp_x'+str(k), 'u_x at {} Hz'.format(freq[k]), lineOn=True )
              makeplots( domain, geom, Lx, Lz, disp_z, 'disp_z'+str(k), 'u_z at {} Hz'.format(freq[k]), lineOn=True )
              makevtk(domain, geom, rho, lam, mu, cp, cs, sol, freq, vec_basis, 'wedge2d')

      elif(ndims==3):
        makevtk(domain, geom, rho, lam, mu, cp, cs, sol.T, freq, vec_basis, 'wedge3d')

util.run( main )