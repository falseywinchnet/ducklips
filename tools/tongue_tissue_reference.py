#!/usr/bin/env python3
"""Constrained sagittal tissue experiment; NOT a speech renderer or calibrated anatomy.

Use a native fixture from tongue_geometry_probe.cpp. A finite dorsal patch is
prescribed toward the palate while lower tissue is anchored. Minimize plane-strain
compressible neo-Hookean energy, constrain total slice area, and prohibit wall
penetration. Element Jacobians expose local compression and inversion. Constant
out-of-plane thickness is assumed, not inferred from a human tongue.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.spatial import Delaunay, ConvexHull


def mesh(n, centre, radius, wall, width):
    u, v = np.meshgrid(np.linspace(-1, 1, n), np.linspace(-1, 1, n))
    uv = np.c_[u.ravel(), v.ravel()]
    disk = uv * np.sqrt(1 - uv[:, ::-1]**2 / 2)
    # Put a material boundary point on the native vertical target ray.
    sine = (wall[0] - centre[0]) / radius
    if abs(sine) >= 1:
        raise ValueError('Fixture target ray misses the resting circle')
    cosine = np.sqrt(1-sine*sine)
    rotation = np.array([[cosine, sine], [-sine, cosine]])
    rest = centre + radius * disk @ rotation.T
    # Pin the physical patch endpoints across refinements. Otherwise changing
    # mesh resolution silently changes the imposed contact region.
    for x in [wall[0]-width/2, wall[0], wall[0]+width/2]:
        if abs(x-centre[0]) >= radius:
            raise ValueError('Dorsal patch extends beyond the resting body')
        point = np.array([x, centre[1]+np.sqrt(radius**2-(x-centre[0])**2)])
        if np.min(np.linalg.norm(rest-point,axis=1)) > 1e-8:
            rest = np.vstack([rest, point])
    triangles = Delaunay(rest).simplices
    boundary = ConvexHull(rest).vertices
    return rest, triangles, boundary


def ceiling(x):
    # Native wall: quarter-circle bend, then piecewise linear palate.
    knots = np.array([2., 4.5, 6., 7.2, 8.4, 9.4])
    heights = 9 + np.array([.9, 1.1, .9, .55, .5, .6])
    h = np.interp(x, knots, heights)
    index = np.clip(np.searchsorted(knots, x, side='right')-1, 0, len(knots)-2)
    slope = np.diff(heights)[index]/np.diff(knots)[index]
    slope[(x >= knots[-1])] = 0
    bend = x < 2
    root = np.sqrt(np.maximum(2.9**2-(x[bend]-2)**2, 0))
    h[bend] = 7+root
    slope[bend] = -(x[bend]-2)/np.maximum(root, 1e-6)
    return h, slope


class RejectedGeometry(RuntimeError):
    def __init__(self, record):
        super().__init__(json.dumps(record))
        self.record = record


class Tissue:
    def __init__(self, fixture, n, width, poisson, incompressible=False):
        self.fixture = fixture
        self.centre = np.array(fixture['centre_cm'])
        self.wall = np.array(fixture['wall_target_cm'])
        self.rest, self.triangles, self.boundary = mesh(n, self.centre, fixture['radius_cm'], self.wall, width)
        p, tri = self.rest, self.triangles
        d = np.stack([p[tri[:, 1]]-p[tri[:, 0]], p[tri[:, 2]]-p[tri[:, 0]]], axis=-1)
        self.inverse = np.linalg.inv(d)
        self.area = np.linalg.det(d)/2
        assert self.area.min() > 0
        self.total_area = self.area.sum()
        self.lam = 2*poisson/(1-2*poisson)  # shear modulus normalized to 1
        # Explicit fixture assumptions, not measured anatomical insertions.
        self.anchor = np.flatnonzero(p[:, 1] <= self.centre[1]-.5*fixture['radius_cm'])
        b = self.boundary
        self.patch = b[(abs(p[b, 0]-self.wall[0]) <= width/2+1e-9) & (p[b, 1] > self.centre[1])]
        if not len(self.patch):
            raise ValueError('No material nodes in dorsal patch')
        self.fixed = np.union1d(self.anchor, self.patch)
        self.free = np.setdiff1d(np.arange(len(p)), self.fixed)
        self.prescribed = p.copy()
        self.incompressible = incompressible
        self.active_elements = np.any(np.isin(tri, self.free), axis=1)

    def unpack(self, z):
        p = self.prescribed.copy()
        p[self.free] = z.reshape(-1, 2)
        return p

    def state(self, z):
        p = self.unpack(z)
        tri = self.triangles
        d = np.stack([p[tri[:, 1]]-p[tri[:, 0]], p[tri[:, 2]]-p[tri[:, 0]]], axis=-1)
        f = d @ self.inverse
        jacobian = np.linalg.det(f)
        return p, f, jacobian

    def assemble(self, stress):
        h = self.area[:, None, None] * (stress @ self.inverse.transpose(0, 2, 1))
        grad = np.zeros_like(self.rest)
        np.add.at(grad, self.triangles[:, 1], h[:, :, 0])
        np.add.at(grad, self.triangles[:, 2], h[:, :, 1])
        np.add.at(grad, self.triangles[:, 0], -h.sum(axis=2))
        return grad[self.free].ravel()

    @staticmethod
    def cofactor(f):
        c = np.empty_like(f)
        c[:, 0, 0] = f[:, 1, 1]
        c[:, 1, 1] = f[:, 0, 0]
        c[:, 0, 1] = -f[:, 1, 0]
        c[:, 1, 0] = -f[:, 0, 1]
        return c

    def energy(self, z):
        p, f, jacobian = self.state(z)
        # Smooth extension outside the admissible region allows line search
        # to return from inverted trial steps. Such steps are never accepted
        # as a valid final geometry; minimum J is checked independently.
        safe = np.maximum(jacobian, .1)
        below = np.minimum(jacobian-.1, 0)/.1
        log = np.log(safe)+below-.5*below**2
        derivative = np.where(jacobian >= .1, 1/safe, (1-below)/.1)
        w = .5*(np.sum(f*f, axis=(1, 2))-2)-log+.5*self.lam*log**2
        stress = f + ((self.lam*log-1)*derivative)[:, None, None]*self.cofactor(f)
        return float(self.area@w), self.assemble(stress)

    def volume(self, z):
        p, f, jacobian = self.state(z)
        return float(self.area@jacobian/self.total_area-1)

    def volume_gradient(self, z):
        p, f, jacobian = self.state(z)
        return self.assemble(self.cofactor(f))/self.total_area

    def element_volume(self, z):
        return self.state(z)[2][self.active_elements]-1

    def element_volume_gradient(self, z):
        p, f, jacobian = self.state(z)
        h = self.cofactor(f) @ self.inverse.transpose(0, 2, 1)
        grad = np.zeros((len(self.triangles), len(p), 2))
        ids = np.arange(len(self.triangles))
        grad[ids, self.triangles[:, 1]] = h[:, :, 0]
        grad[ids, self.triangles[:, 2]] = h[:, :, 1]
        grad[ids, self.triangles[:, 0]] = -h.sum(axis=2)
        return grad[self.active_elements][:, self.free].reshape(np.count_nonzero(self.active_elements), -1)

    def contact(self, z):
        p = self.unpack(z)[self.free]
        height, slope = ceiling(p[:, 0])
        return np.r_[p[:, 0]+.9, height-p[:, 1]]

    def contact_gradient(self, z):
        p = self.unpack(z)[self.free]
        height, slope = ceiling(p[:, 0])
        n = len(p)
        result = np.zeros((2*n, 2*n))
        result[np.arange(n), 2*np.arange(n)] = 1
        result[n+np.arange(n), 2*np.arange(n)] = slope
        result[n+np.arange(n), 2*np.arange(n)+1] = -1
        return result

    def solve(self, previous, activation):
        amount = activation**2*(3-2*activation)
        patch = self.patch
        self.prescribed = self.rest.copy()
        target_y = ceiling(self.rest[patch, 0])[0]-self.fixture['gap_cm']
        self.prescribed[patch, 1] += amount*(target_y-self.rest[patch, 1])
        result = minimize(self.energy, previous[self.free].ravel(), jac=True, method='SLSQP',
                          constraints=[{'type':'eq', 'fun':self.element_volume if self.incompressible else self.volume,
                                        'jac':self.element_volume_gradient if self.incompressible else self.volume_gradient},
                                       {'type':'ineq', 'fun':self.contact, 'jac':self.contact_gradient}],
                          options={'maxiter':1200, 'ftol':1e-10})
        p, f, jacobian = self.state(result.x)
        h = ceiling(p[:, 0])[0]
        record = {'activation':float(activation), 'converged':bool(result.success),
                  'iterations':int(result.nit), 'message':result.message,
                  'normalized_elastic_energy':float(result.fun),
                  'slice_area_relative_error':self.volume(result.x),
                  'element_jacobian_min':float(jacobian.min()),
                  'element_jacobian_max':float(jacobian.max()),
                  'wall_penetration_cm':float(max(0, np.max(p[:, 1]-h), np.max(-.9-p[:, 0]))),
                  'patch_clearance_cm':(h[patch]-p[patch, 1]).tolist(),
                  'anchor_max_displacement_cm':float(np.linalg.norm(p[self.anchor]-self.rest[self.anchor],axis=1).max())}
        if (not result.success or abs(record['slice_area_relative_error'])>1e-8
                or jacobian.min()<=0 or record['wall_penetration_cm']>1e-7
                or (self.incompressible and np.max(abs(jacobian-1))>1e-8)):
            raise RejectedGeometry(record)
        constraints = [self.element_volume_gradient(result.x) if self.incompressible
                       else self.volume_gradient(result.x)[None, :]]
        active = self.contact(result.x)<1e-7
        if active.any():
            constraints.append(self.contact_gradient(result.x)[active])
        matrix = np.concatenate(constraints)
        gradient = self.energy(result.x)[1]
        multiplier = np.linalg.lstsq(matrix.T, -gradient, rcond=1e-10)[0]
        record['projected_stationarity_inf'] = float(np.max(abs(gradient+matrix.T@multiplier)))
        if record['projected_stationarity_inf']>1e-3:
            raise RejectedGeometry(record)
        return p, record

    def verify_derivatives(self):
        rng = np.random.default_rng(481)
        z = self.rest[self.free].ravel()+rng.normal(0, .003, 2*len(self.free))
        direction = rng.normal(size=len(z)); direction /= np.linalg.norm(direction)
        step = 1e-6
        energy = abs((self.energy(z+step*direction)[0]-self.energy(z-step*direction)[0])/(2*step)-self.energy(z)[1]@direction)
        volume = abs((self.volume(z+step*direction)-self.volume(z-step*direction))/(2*step)-self.volume_gradient(z)@direction)
        contact = np.max(abs((self.contact(z+step*direction)-self.contact(z-step*direction))/(2*step)-self.contact_gradient(z)@direction))
        local = np.max(abs((self.element_volume(z+step*direction)-self.element_volume(z-step*direction))/(2*step)-self.element_volume_gradient(z)@direction))
        if max(energy, volume, contact, local)>1e-6:
            raise RuntimeError(f'derivative mismatch: {energy}, {volume}, {contact}')
        return {'energy_absolute_error':float(energy),'area_absolute_error':float(volume),'contact_absolute_error':float(contact),'local_area_absolute_error':float(local)}


def ray_polygon(wall, inward, vertices):
    start = vertices
    edge = np.roll(vertices, -1, axis=0)-vertices
    offset = start-wall
    # A pressed body already containing the wall origin seals the ray. The
    # positive exit intersection is not an air gap.
    end = np.roll(vertices, -1, axis=0)
    crosses = (start[:, 1] > wall[1]) != (end[:, 1] > wall[1])
    crossing_x = np.full(len(vertices), -np.inf)
    crossing_x[crosses] = start[crosses, 0] + (wall[1]-start[crosses, 1])*edge[crosses, 0]/edge[crosses, 1]
    if np.count_nonzero(crosses & (wall[0] < crossing_x)) % 2:
        return 0.0
    denominator = inward[0]*edge[:, 1]-inward[1]*edge[:, 0]
    valid = abs(denominator)>1e-12
    t = np.full(len(edge), np.inf); u = t.copy()
    t[valid] = (offset[valid, 0]*edge[valid, 1]-offset[valid, 1]*edge[valid, 0])/denominator[valid]
    u[valid] = (offset[valid, 0]*inward[1]-offset[valid, 1]*inward[0])/denominator[valid]
    hit = valid & (t >= -1e-8) & (u >= 0) & (u <= 1)
    return float(np.maximum(t[hit], 0).min()) if hit.any() else np.inf


def body_airway(fixture, vertices):
    # Deliberately isolate the BODY. No fitted residuals, lips, tip or blade.
    sections = np.array(fixture['sections'])
    arc = sections[:, 0]
    gaps = np.array([min(row[5], ray_polygon(row[1:3], -row[3:5], vertices)) for row in sections])
    blend = np.clip((arc-7)/np.pi, 0, 1)
    area = (1-blend)*1.6*gaps**1.35 + blend*1.9*gaps**1.45
    return float(np.trapezoid(area, arc)), gaps


def run(args):
    fixture = json.loads(args.fixture.read_text())
    if not np.allclose(fixture['normal'], [0, 1]):
        raise ValueError('This reference currently supports oral, upward closures only')
    if fixture['gap_cm'] != 0:
        raise ValueError('This experiment tests complete closure, not calibrated finite-area grooves')
    square = np.array([[-1., -1.], [1., -1.], [1., 1.], [-1., 1.]])
    assert ray_polygon(np.array([0., 0.]), np.array([0., -1.]), square) == 0
    assert abs(ray_polygon(np.array([0., 3.]), np.array([0., -1.]), square)-2) < 1e-12
    assert np.isinf(ray_polygon(np.array([3., 3.]), np.array([0., -1.]), square))
    sections = np.array(fixture['sections'])
    upper = sections[sections[:, 0] >= 7]
    assert np.max(abs(ceiling(upper[:, 1])[0]-upper[:, 2])) < 1e-8
    args.out.mkdir(parents=True, exist_ok=True)
    receipt = {'fixture_sha256':hashlib.sha256(args.fixture.read_bytes()).hexdigest(),
               'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               'model':'quasistatic plane-strain neo-Hookean triangle mesh, selected area constraints, unilateral outer wall contact',
               'assumptions':{'poisson_ratio':args.poisson, 'shear_modulus':'normalized; no force calibration',
                              'elementwise_incompressible':args.incompressible, 'continuation_steps':args.steps, 'patch_width_cm':args.width,'anchor':'nodes below centre minus half radius',
                              'out_of_plane_thickness':'constant, unmeasured'},
               'scope':'geometry reference only; not integrated speech and not validated human tissue', 'meshes':[]}
    plot_data = []
    for n in args.mesh:
        tissue = Tissue(fixture, n, args.width, args.poisson, args.incompressible)
        entry = {'grid_size':n,'nodes':len(tissue.rest),'triangles':len(tissue.triangles),
                 'patch_nodes':len(tissue.patch),'actual_patch_span_cm':float(np.ptp(tissue.rest[tissue.patch,0])),'rest_slice_area_cm2':float(tissue.total_area),
                 'derivatives':tissue.verify_derivatives(),'steps':[]}
        p = tissue.rest.copy()
        for activation in np.linspace(0, 1, args.steps+1):
            try:
                p, record = tissue.solve(p, activation)
            except RejectedGeometry as failure:
                entry['rejected_trial'] = failure.record
                receipt['meshes'].append(entry)
                receipt['numerical_status'] = 'rejected; invalid trial was not used as a geometry'
                np.savez(args.out/f'last-valid-{n}.npz',rest=tissue.rest,last_valid=p,triangles=tissue.triangles,boundary=tissue.boundary)
                (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
                print(json.dumps({'grid_size':n,'rejected_trial':failure.record}),flush=True)
                return 1
            volume, gaps = body_airway(fixture, p[tissue.boundary])
            record['body_only_airway_volume_ml'] = volume
            entry['steps'].append(record)
        shift = np.array(fixture['poses'][-1]['centre_cm'])-fixture['centre_cm']
        rigid_volume, rigid_gaps = body_airway(fixture, tissue.rest[tissue.boundary]+shift)
        entry['rigid_body_only_airway_volume_ml'] = rigid_volume
        rigid = tissue.rest[tissue.boundary]+shift
        entry['rigid_wall_penetration_cm'] = float(np.max(rigid[:,1]-ceiling(rigid[:,0])[0]))
        entry['rest_to_closure_airway_loss_ml'] = entry['steps'][0]['body_only_airway_volume_ml']-entry['steps'][-1]['body_only_airway_volume_ml']
        entry['rigid_airway_loss_ml'] = entry['steps'][0]['body_only_airway_volume_ml']-rigid_volume
        np.savez(args.out/f'tissue-{n}.npz',rest=tissue.rest,closed=p,triangles=tissue.triangles,boundary=tissue.boundary,anchor=tissue.anchor,patch=tissue.patch)
        receipt['meshes'].append(entry)
        plot_data.append((tissue,p))
        print(json.dumps({'grid_size':n,'closure':entry['steps'][-1],'airway_loss_ml':entry['rest_to_closure_airway_loss_ml'],'rigid_airway_loss_ml':entry['rigid_airway_loss_ml']}),flush=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1,len(plot_data),figsize=(5*len(plot_data),5.8),squeeze=False)
    sections=np.array(fixture['sections'])
    for ax,(tissue,p) in zip(axes[0],plot_data):
        boundary=np.r_[tissue.boundary,tissue.boundary[0]]
        ax.plot(sections[:,1],sections[:,2],color='black',label='Outer wall')
        ax.plot(tissue.rest[boundary,0],tissue.rest[boundary,1],color='gray',linestyle='--',label='Rest body')
        ax.triplot(p[:,0],p[:,1],tissue.triangles,color='#b45129',lw=.5)
        ax.plot(p[boundary,0],p[boundary,1],color='#b45129',label='Deformed body')
        shift=np.array(fixture['poses'][-1]['centre_cm'])-fixture['centre_cm']
        rigid=tissue.rest[boundary]+shift
        ax.plot(rigid[:,0],rigid[:,1],color='#397bb5',label='Current rigid gesture')
        ax.scatter(p[tissue.anchor,0],p[tissue.anchor,1],s=12,color='black',label='Fixed lower tissue')
        ax.scatter(p[tissue.patch,0],p[tissue.patch,1],s=16,color='#cc208e',label='Driven dorsal patch')
        ax.set(xlabel='Anterior position (cm)',ylabel='Height (cm)',title=f'{len(tissue.rest)} nodes',xlim=(-1.2,9.6),ylim=(1.7,10.5),aspect='equal')
        ax.legend(fontsize=8,loc='lower right')
    fig.suptitle('Quebec closure: constrained sagittal tissue experiment — not a speech render')
    fig.tight_layout();fig.savefig(args.out/'tissue-closure.png',dpi=160);plt.close(fig)
    receipt['numerical_status'] = 'all solves passed; anatomical validity is unproven'
    (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return 0


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fixture',type=Path)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--mesh',type=int,nargs='+',default=[7,11,15])
    parser.add_argument('--width',type=float,default=.9)
    parser.add_argument('--poisson',type=float,default=.49)
    parser.add_argument('--steps',type=int,default=10)
    parser.add_argument('--incompressible',action='store_true',help='Conserve every triangle area, not just their sum')
    args=parser.parse_args()
    if any(n<5 or n%2==0 for n in args.mesh) or not 0<args.poisson<.5 or args.width<=0 or args.steps<1:
        parser.error('Require odd mesh sizes >=5, 0<poisson<.5, positive patch width')
    raise SystemExit(run(args))
