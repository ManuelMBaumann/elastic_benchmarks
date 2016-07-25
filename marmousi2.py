from __future__ import print_function, division

import matplotlib
matplotlib.use('agg')

from nutils import *
import scipy
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from obspy.segy.core import readSEGY
from mpl_toolkits.axes_grid1 import make_axes_locatable


@log.title
def makeplots( domain, geom, verts_x, verts_z, value, name, title, ndigits=0, index=None, imgtype=None):
  points, colors = domain.elem_eval( [ geom, value ], ischeme='bezier3', separate=True )

  with plot.PyPlot( name, ndigits=ndigits, index=index, imgtype=imgtype  ) as plt:

      plt.mesh( points, colors, triangulate='bezier', edgecolors='none' )
     
      plt.title(title)
          
      plt.xlabel('x [m]', fontsize=10)
      plt.xticks([0, max(verts_x)/2.0, max(verts_x)], ['0', '2000', '4000'], fontsize=10)
      plt.ylabel('z [m]', fontsize=10)
      plt.yticks([max(verts_z), min(verts_z)], ['0', '1850'], fontsize=10)
      plt.ylim
          
      ax = plt.gca()
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05)
      cb = plt.colorbar(cax=cax)
      cb.ax.tick_params(labelsize=10)   
          
          
def makespyplot( matrix, name, imgtype=None ):
  if not scipy.sparse.isspmatrix( matrix ):
      matrix = matrix.toscipy()

  with plot.PyPlot( name, ndigits=0, imgtype=imgtype ) as plt:
    plt.spy( matrix, markersize=0.8, color='black')
    plt.title( name+', nnz = '+str(matrix.nnz) )


def makevtk(domain, geom, rho, lam, mu, cp, cs, sol, freq, vec_basis, name):
  Nom = sol.shape[0]
  vtk_geom, vtk_rho, vtk_lam, vtk_mu, vtk_cp, vtk_cs = domain.simplex.elem_eval( [ geom, rho, lam, mu, cp, cs ], ischeme='vtk', separate=True )
  with plot.VTKFile( name ) as vtk:
      vtk.unstructuredgrid( vtk_geom )
      vtk.pointdataarray( 'rho', vtk_rho )
      vtk.pointdataarray( 'lambda', vtk_lam )
      vtk.pointdataarray( 'mu', vtk_mu )
      vtk.pointdataarray( 'cp', vtk_cp )
      vtk.pointdataarray( 'cs', vtk_cs )
      for i in range(0,Nom):
          disp = vec_basis.dot( sol[i,:] ).real
          vtk_disp = domain.simplex.elem_eval( disp, ischeme='vtk', separate=True )
          vtk.pointdataarray( 'disp_'+str(int(freq[i])), vtk_disp )
          

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

  # build matrices
  K, M = domain.integrate( [elasticity, mass], geometry=geom, ischeme='gauss2' )
  C = domain.boundary['left,right,bottom'].integrate( sommerfeld, geometry=geom, ischeme='gauss2' )
  
  # build RHS
  source_position = nx//3, nz
  rhs = point_eval(vec_basis, domain, geom, source_position)[:,-1]
  #rhs = rhs+point_eval(vec_basis, domain, geom, source_position)[:,0]
  
  return K, C, M, rhs    

def shrink_marmousi(segy1, segy2, segy3, n_course, x0, xe, ze):
  # original problem dimensions
  Lx = 17000.0
  Lz = 3500.0
  nx = len(segy1.traces)
  nz = len(segy1.traces[0].data)
  dx = Lx/(nx-1)
  dz = Lz/(nz-1) 

  # cast segy data into arrays
  rho_coeffs = np.zeros((nz, nx))
  cp_coeffs  = np.zeros((nz, nx))
  cs_coeffs  = np.zeros((nz, nx))

  for i, tr in enumerate(segy1.traces):
    rho_coeffs[:,i] = tr.data
  for i, tr in enumerate(segy2.traces):
    cp_coeffs[:,i] = tr.data
  for i, tr in enumerate(segy3.traces):
    cs_coeffs[:,i] = tr.data
  rho_coeffs *= 1000.0

  # Shrink computational domain
  z0_ind = 370 # delete water layer
  ze_ind = 0
  x0_ind = 0
  xe_ind = 0
  
  z0 = z0_ind*dz 
  ze = ze + z0
    
  for i in range(z0_ind,nz):
      if (i*dz>ze):
          ze_ind = i
          break

  for i in range(0,nx):
      if (i*dx>x0):
          x0_ind = i
          break
      
  for i in range(x0_ind,nx):
      if (i*dx>xe):
          xe_ind = i
          break
      
  rho_coeffs = rho_coeffs[z0_ind:ze_ind,x0_ind:xe_ind]
  cp_coeffs  =  cp_coeffs[z0_ind:ze_ind,x0_ind:xe_ind]
  cs_coeffs  =  cs_coeffs[z0_ind:ze_ind,x0_ind:xe_ind] 
  
  rho_coeffs = rho_coeffs[::-n_course, ::n_course]
  cp_coeffs  =  cp_coeffs[::-n_course, ::n_course]
  cs_coeffs  =  cs_coeffs[::-n_course, ::n_course]
    
  nx = rho_coeffs.shape[1]
  nz = rho_coeffs.shape[0]
  Lx = xe-x0
  Lz = ze-z0
  
  # Coursen the data
  dx *= n_course
  dz *= n_course
  verts_x = np.linspace( 0.0, Lx+dz, nx )
  verts_z = np.linspace( -Lz-z0, -z0+dz, nz )
  
  rho_coeffs = rho_coeffs.T.ravel()
  cp_coeffs  = cp_coeffs.T.ravel()
  cs_coeffs  = cs_coeffs.T.ravel()

  return verts_x, verts_z, rho_coeffs, cp_coeffs, cs_coeffs
    

def main( n_course=16,      # coursening of original problem
          freq=[4.0],       # frequency in Hz 
          plots=True,       # plot of parameters and num. solution
          spy=True,         # provide spy plot of matrices
          storing=False,    # storing=True saves matrices in matrix market format (no solve) 
          degree=1 ):       # degree of FEM splines

  
  ndims  = 2
  freq   = np.array(freq)
  om     = 2.0*np.pi*freq
  Nom    = len(om)

  print("Reading Marmousi-II data...\n")
  segy1 = readSEGY("data/MODEL_DENSITY_1.25m.segy")
  segy2 = readSEGY("data/MODEL_P-WAVE_VELOCITY_1.25m.segy")
  segy3 = readSEGY("data/MODEL_S-WAVE_VELOCITY_1.25m.segy")

  # shrink domain to [0,4000] x [0,1850]
  [verts_x, verts_z, rho_coeffs, cp_coeffs, cs_coeffs] = shrink_marmousi(segy1, segy2, segy3, n_course, 7500, 11500, 1850 )

  mu_coeffs  = cs_coeffs**2 * rho_coeffs
  lam_coeffs = rho_coeffs * (cp_coeffs**2 - 2.0*cs_coeffs**2) 

  nx = len(verts_x)
  nz = len(verts_z)
  dx = verts_x[1]-verts_x[0]
  dz = verts_z[1]-verts_z[0]

  # define Cartesian grid
  domain, geom = mesh.rectilinear( [verts_x, verts_z] )     
  vec_basis    = domain.splinefunc( degree=degree ).vector( ndims )
  scal_basis   = domain.splinefunc( degree=1 )

  # problem summary
  print( '--- MARMOUSI-II PROBLEM  ---' )
  print( 'problem size   : '+str(nx-1+degree)+' x '+str(nz-1+degree) )
  print( '# dofs         : '+str(len(vec_basis)) )
  print( 'grid size      : '+str(round(dx,1))+' x '+str(round(dz,1)) )
  ppw = 20.0  
  print( 'max. frequency : '+str( round(min(np.amin(cp_coeffs),np.amin(cs_coeffs))/(ppw*max(dx,dz)),1) ) )
  print( '----------------------------\n' )
  
  # Create discretization matrices using nutils
  mu  = scal_basis.dot(mu_coeffs)
  lam = scal_basis.dot(lam_coeffs)
  rho = scal_basis.dot(rho_coeffs)
  cp  = scal_basis.dot(cp_coeffs)
  cs  = scal_basis.dot(cs_coeffs)
  
  K, C, M, rhs = elast_mat(rho, cp, cs, lam, mu, ndims, nx, 1, nz, vec_basis, domain, geom)

  if storing:
    mmwrite('matrix_io/K.mtx', K.toscipy())
    mmwrite('matrix_io/C.mtx', C.toscipy())
    mmwrite('matrix_io/M.mtx', M.toscipy())
  else:
    print('Use pythons sparse linear solver solver...')
    sol = np.zeros((Nom, len(vec_basis)), dtype=complex)
    for k in range(0,Nom):
      matrix = K + 1j*om[k]*C - om[k]**2*M  
      if spy:
        makespyplot( matrix, 'spy_plot' )
      sol[k,:] = scipy.sparse.linalg.spsolve( matrix.toscipy().tocsc(), rhs )
          
    if plots:
      makeplots( domain, geom, verts_x, verts_z, rho, 'rho', 'rho [kg/m**3]')
      #makeplots( domain, geom, verts_x, verts_z, cp, 'cp', 'c_p [m/s]')
      #makeplots( domain, geom, verts_x, verts_z, cs, 'cs', 'c_s [m/s]')
      
      for k in range(0,Nom):
          disp     = vec_basis.dot( sol[k,:] )     # FEM summation
          disp_x   = disp[0].real                  # Plot Re(u_x)  
          disp_z   = disp[-1].real                 # Plot Re(u_z)
          makeplots( domain, geom, verts_x, verts_z, disp_x, 'disp_x'+str(k), 'u_x at {} Hz'.format(freq[k]) )
          makeplots( domain, geom, verts_x, verts_z, disp_z, 'disp_z'+str(k), 'u_z at {} Hz'.format(freq[k]) )
          makevtk(domain, geom, rho, lam, mu, cp, cs, sol, freq, vec_basis, 'marmousi2')
                  
util.run( main )