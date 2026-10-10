#!/usr/bin/env python3
"""Compile an elastic volume into a small pressure-driven deformation engine.

Eight-node hexahedra, three displacement components, fixed lateral backing and
anterior/posterior attachments. Selective integration of bulk strain avoids
volumetric locking near incompressibility; deviatoric strain is fully integrated.
Only body setup uses an eigenproblem. No eigenfrequency is a commanded pitch.
"""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
from scipy.sparse import coo_matrix,diags
from scipy.sparse.linalg import eigsh,spsolve

ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--config',type=Path,default=ROOT/'data/tissue_body.json')
p.add_argument('--out',type=Path,default=ROOT/'out/compiled-tissue-body')
p.add_argument('--header',type=Path,default=ROOT/'include/ducklips/tissue_body_data.hpp')
args=p.parse_args();c=json.loads(args.config.read_text());args.out.mkdir(parents=True,exist_ok=True)
nx,nz=c['length_nodes'],c['height_nodes'];nb,nc=c['body_depth_cells'],c['cover_depth_cells'];ny=nb+nc+1
xs=np.linspace(0,c['length'],nx);zs=np.linspace(0,c['height'],nz)
ys=np.r_[np.linspace(0,c['body_depth'],nb+1),c['body_depth']+np.linspace(0,c['cover_depth'],nc+1)[1:]]
xyz=np.array([[x,y,z] for z in zs for y in ys for x in xs])
def node(x,y,z):return (z*ny+y)*nx+x
signs=np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],[-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]])
def shape(point,lengths):
    f=1+signs*point
    N=np.prod(f,axis=1)/8
    grad=np.column_stack([signs[:,j]*f[:,(j+1)%3]*f[:,(j+2)%3]/(4*lengths[j]) for j in range(3)])
    B=np.zeros((6,24))
    for a,(dx,dy,dz) in enumerate(grad):
        j=3*a;B[0,j]=dx;B[1,j+1]=dy;B[2,j+2]=dz
        B[3,j]=dy;B[3,j+1]=dx;B[4,j+1]=dz;B[4,j+2]=dy;B[5,j]=dz;B[5,j+2]=dx
    return N,grad,B

rows=[];cols=[];ks=[];ds=[];mass=np.zeros(3*len(xyz))
q=1/np.sqrt(3);points=np.array([[a,b,d] for a in [-q,q] for b in [-q,q] for d in [-q,q]])
dev=np.zeros((6,6));dev[:3,:3]=2*(np.eye(3)-np.ones((3,3))/3);dev[3:,3:]=np.eye(3)
for iz in range(nz-1):
 for iy in range(ny-1):
  for ix in range(nx-1):
    nodes=np.array([node(ix+(s[0]+1)//2,iy+(s[1]+1)//2,iz+(s[2]+1)//2) for s in signs])
    dofs=(3*nodes[:,None]+np.arange(3)).ravel()
    lengths=np.array([xs[ix+1]-xs[ix],ys[iy+1]-ys[iy],zs[iz+1]-zs[iz]])
    volume=np.prod(lengths);jac=volume/8
    cover=iy>=nb
    G=c['cover_shear' if cover else 'body_shear'];eta=c['cover_viscosity' if cover else 'body_viscosity']
    stress=c['longitudinal_stress' if cover else 'body_longitudinal_stress']
    bulk=2*G*(1+c['poisson_ratio'])/(3*(1-2*c['poisson_ratio']))
    K=np.zeros((24,24));D=np.zeros((24,24))
    for point in points:
        N,grad,B=shape(point,lengths)
        K+=jac*(G*(B.T@dev@B)+stress*np.kron(np.outer(grad[:,0],grad[:,0]),np.eye(3)))
        D+=jac*eta*(B.T@dev@B)
        mass[dofs]+=np.repeat(c['density']*N*jac,3)
    _,_,B=shape(np.zeros(3),lengths);div=B[:3].sum(axis=0)
    K+=bulk*volume*np.outer(div,div)
    rows.extend(np.repeat(dofs,24));cols.extend(np.tile(dofs,24));ks.extend(K.ravel());ds.extend(D.ravel())
n=len(mass);K=coo_matrix((ks,(rows,cols)),shape=(n,n)).tocsr();D=coo_matrix((ds,(rows,cols)),shape=(n,n)).tocsr()
fixed=np.repeat((xyz[:,0]==0)|(xyz[:,0]==c['length'])|(xyz[:,1]==0),3)
free=np.where(~fixed)[0];inverse_mass=1/np.sqrt(mass[free]);scale=diags(inverse_mass)
Kr=scale@K[free][:,free]@scale
nm=c['retained_modes'];eigenvalues,eigenvectors=eigsh(Kr,k=nm,sigma=0,which='LM',v0=np.linspace(1,2,len(free)),tol=1e-11)
order=np.argsort(eigenvalues);eigenvalues=eigenvalues[order];eigenvectors=eigenvectors[:,order]
assert np.all(eigenvalues>0)
phi=np.zeros((n,nm));phi[free]=inverse_mass[:,None]*eigenvectors
for i in range(nm):
    if phi[np.argmax(abs(phi[:,i])),i]<0:phi[:,i]*=-1
C=phi.T@D@phi;C=.5*(C+C.T)
residual=K[free]@phi-eigenvalues*(mass[free,None]*phi[free])
relative_residual=float(np.linalg.norm(residual)/np.linalg.norm(K[free]@phi))
orthogonality=float(np.max(abs(phi.T@(mass[:,None]*phi)-np.eye(nm))))
assert relative_residual<1e-7 and orthogonality<1e-8 and np.linalg.eigvalsh(C).min()>-1e-8
surface=np.array([node(x,ny-1,z) for z in range(nz) for x in range(nx)])
inside=np.array([node(x,nb//2,z) for z in range(nz) for x in range(nx)])
surface_basis=phi[3*surface+1];inside_basis=phi[3*inside+1]
area=np.outer(np.r_[.5,np.ones(nz-2),.5],np.r_[.5,np.ones(nx-2),.5]).ravel()*(xs[1]-xs[0])*(zs[1]-zs[0])
# Static posturing influence: prescribe unit normal medial displacement,
# allow interior/vertical/longitudinal displacement to equilibrate, then
# project the resulting volume deformation. Runtime posture changes only
# the equilibrium coordinates, preserving physical position and momentum.
posture=np.zeros(n);along=xyz[surface,0]/c['length'];posture[3*surface+1]=4*along*(1-along)
static_fixed=fixed.copy();static_fixed[3*surface+1]=True
static_free=np.where(~static_fixed)[0]
posture[static_free]=spsolve(K[static_free][:,static_free],-(K[static_free]@posture))
posture_coordinates=phi.T@(mass*posture)
posture_error=float(np.linalg.norm(surface_basis@posture_coordinates-posture[3*surface+1])/np.linalg.norm(posture[3*surface+1]))
np.savez_compressed(args.out/'body.npz',xyz=xyz,modes=phi.reshape(-1,3,nm),eigenvalues=eigenvalues,damping=C,mass=mass,surface=surface,inside=inside)
receipt={'config':c,'config_sha256':hashlib.sha256(args.config.read_bytes()).hexdigest(),
 'compiler_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'free_degrees':len(free),
 'retained_modes':nm,'frequencies_hz':(np.sqrt(eigenvalues)/(2*np.pi)).tolist(),
 'relative_eigen_residual':relative_residual,'mass_orthogonality_error':orthogonality,'posture_surface_relative_error':posture_error,
 'method':'linear isotropic volume elasticity; nearly incompressible selective bulk integration; initial longitudinal stress; projected shear viscosity',
 'limits':['linear strain about reference posture','symmetric folds','no moving fluid volume or acoustic feedback yet']}
(args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
def array(name,a):
    flat=np.asarray(a).ravel()
    lines=[', '.join(format(v,'.17g') for v in flat[i:i+8]) for i in range(0,len(flat),8)]
    return f'inline constexpr std::array<double,{len(flat)}> {name}{{{{\n    '+',\n    '.join(lines)+'\n}};\n'
header='// Generated by tools/compile_tissue_body.py from data/tissue_body.json.\n#pragma once\n#include <array>\n#include <cstddef>\nnamespace ducklips::compiled_tissue {\n'
header+=f'inline constexpr std::size_t nx={nx},nz={nz},modes={nm};\n'
for name in ['length','height','cover_depth','body_depth','density','cover_shear','body_shear','cover_viscosity','body_viscosity','longitudinal_stress','body_longitudinal_stress']:
    header+=f'inline constexpr double {name}={c[name]:.17g};\n'
header+=array('stiffness',eigenvalues)+array('damping',C)+array('surface_basis',surface_basis)+array('inside_basis',inside_basis)+array('surface_area',area)+array('posture',posture_coordinates)+'}\n'
args.header.write_text(header)
print(json.dumps(receipt,indent=2))
