#!/usr/bin/env python3
"""Energy-normalized dynamic realization of the frozen surface scattering mesh.

Research tool: no geometry, paths, losses, or radiators are removed. Orthogonal
projection of this contractive scattering realization preserves passivity.
"""
import argparse
import json
import time
from pathlib import Path
import numpy as np
from scipy.sparse import block_diag, coo_matrix, csr_matrix, diags, eye
from scipy.sparse.linalg import spsolve
from surface_ports import FrozenPorts


def first_order(b, a):
    """Balanced scalar passive realization; certify the small scattering map."""
    pole, direct = -a[1], b[0]
    residue = b[1] + direct*pole
    root = np.sqrt(abs(residue))
    return (np.array([[pole]]), np.array([[root]]),
            np.array([[np.copysign(root, residue)]]), direct)


def cascade(one, two):
    a, b, c, d = one
    e, f, g, h = two
    return (np.block([[a, np.zeros((len(a), len(e)))], [f@c, e]]),
            np.vstack([b, f*d]), np.hstack([h*c, g]), h*d)


def delay_line(length):
    if length < 1:
        raise ValueError('Surface paths must have at least one integer sample of delay')
    a = np.zeros((length, length))
    a[1:, :-1] = np.eye(length-1)
    b = np.zeros((length, 1)); b[0] = 1
    c = np.zeros((1, length)); c[0, -1] = 1
    return a, b, c, 0.0


def path_filter(path):
    b, a = np.array(path[4]), np.array(path[5])
    if abs(a[2]) + abs(a[1]) < 1e-15:
        # Actual linear fractional delay (convex two-tap interpolation).
        return first_order(b, a)
    poles = np.roots(a)
    if np.max(abs(poles.imag)) > 1e-9:
        raise ValueError('Expected real all-pass and shelf poles')
    poles = poles.real
    # An all-pass pole has a reciprocal zero: evaluate the numerator at that
    # zero without dividing by the pole (also handles an all-pass pole at 0).
    errors = [abs(b[0] + b[1]*p + b[2]*p*p) for p in poles]
    which = int(np.argmin(errors))
    allpass = -poles[which]
    shelf_pole = poles[1-which]
    shelf_b1 = b[2]
    shelf_b0 = b[1] - allpass*shelf_b1
    if not np.allclose([allpass*shelf_b0, shelf_b0+allpass*shelf_b1, shelf_b1], b, atol=1e-12, rtol=1e-10):
        raise ValueError('Cannot separate exported path filters without changing them')
    root = np.sqrt(1-allpass*allpass)
    ap = (np.array([[-allpass]]), np.array([[root]]), np.array([[root]]), allpass)
    shelf = first_order([shelf_b0,shelf_b1,0], [1,-shelf_pole,0])
    return cascade(ap, shelf)


def wall_filter(port, fs):
    """Recover normalized mass/resistance/stiffness from the original biquad.

    The wall is a mass/spring/damper driven through its acoustic wave port.
    Its bilinear state map is contractive in the mechanical energy metric.
    """
    n, d = np.array(port[3]), np.array(port[4])
    h = 0.5/fs
    def continuous(poly):
        return np.array([poly.sum(), 2*h*(poly[0]-poly[2]), h*h*(poly[0]-poly[1]+poly[2])])
    numerator, denominator = continuous(d+n), continuous(d-n)
    mass, resistance, stiffness = numerator[[2,1,0]] / denominator[1]
    if min(mass, resistance, stiffness) <= 0:
        raise ValueError('Nonpositive wall energy parameter')
    omega = np.sqrt(stiffness/mass)
    ac = np.array([[0,omega],[-omega,-(resistance+1)/mass]])
    bc = np.array([[0],[np.sqrt(2/mass)]])
    cc = -bc.T
    inverse = np.linalg.inv(np.eye(2)-ac/(2*fs))
    ad = inverse@(np.eye(2)+ac/(2*fs))
    bd = inverse@bc/np.sqrt(fs)
    cd = cc@inverse/np.sqrt(fs)
    dd = float((1 + cc@inverse@bc/(2*fs)).item())
    return ad, bd, cd, dd


def certify(device):
    a,b,c,d = device
    step = np.block([[a,b],[c,np.array([[d]])]])
    gain = np.linalg.svd(step,compute_uv=False)[0]
    if gain > 1+2e-10 or not np.isfinite(gain):
        raise ValueError(f'Active component: energy-normalized step singular value {gain}')
    return float(gain)


class PassiveSurface:
    def __init__(self, snapshot):
        self.snapshot = snapshot
        self.fs = snapshot['sample_rate']
        self.devices = []
        self.paths = len(snapshot['paths'])
        from_node, to_node, admittance, reverse = [], [], [], []
        self.component_gain = 0.0
        self.component_filter_error = 0.0
        for index, path in enumerate(snapshot['paths']):
            device = cascade(delay_line(path[3]), path_filter(path))
            self.devices.append(device)
            from_node.append(path[0]); to_node.append(path[1]); admittance.append(path[2]); reverse.append(index^1)
        for port in snapshot['ports']:
            device = first_order(port[3],port[4]) if port[2] else wall_filter(port,self.fs)
            reverse.append(len(self.devices))
            self.devices.append(device)
            from_node.append(port[0]); to_node.append(port[0]); admittance.append(port[1])
        frequencies = [48,500,3000,8000,12000]
        for index, device in enumerate(self.devices):
            self.component_gain = max(self.component_gain,certify(device))
            a,b,c,d = device
            source = snapshot['paths'][index] if index < self.paths else snapshot['ports'][index-self.paths]
            numerator, denominator = (source[4],source[5]) if index < self.paths else (source[3],source[4])
            for frequency in frequencies:
                z = np.exp(2j*np.pi*frequency/self.fs); q=1/z
                expected = FrozenPorts.transfer(np.array([numerator]),np.array([denominator]),q)[0]
                if index < self.paths: expected *= q**source[3]
                actual = (c@np.linalg.solve(z*np.eye(len(a))-a,b)).item()+d
                self.component_filter_error = max(self.component_filter_error,float(abs(expected-actual)))
        if self.component_filter_error > 1e-8:
            raise ValueError(f'Component transfer mismatch: {self.component_filter_error}')
        count = len(self.devices)
        self.from_node, self.to_node = np.array(from_node), np.array(to_node)
        self.y = np.array(admittance); self.root_y = np.sqrt(self.y)
        self.direct = np.array([device[3] for device in self.devices])
        self.offset = np.r_[0,np.cumsum([len(device[0]) for device in self.devices])]
        self.states = int(self.offset[-1])
        tubes=snapshot.get('piece_tubes',[0]*snapshot['pieces'])
        names=[tuple(sorted((tubes[a],tubes[b]))) for a,b in zip(from_node,to_node)]
        groups={name:i for i,name in enumerate(sorted(set(names)))}
        self.state_groups=np.repeat([groups[name] for name in names],np.diff(self.offset))
        self.external_nodes = np.array(snapshot['pressure_pieces'],dtype=int)
        if len(np.unique(self.external_nodes)) != len(self.external_nodes):
            raise ValueError('External observation/source ports must be distinct')
        node_sum = np.bincount(self.to_node,weights=self.y,minlength=snapshot['pieces'])
        # Reference impedances define wave coordinates only. Closing these
        # external ports with prescribed flows recovers the original network.
        self.external_y = node_sum[self.external_nodes].copy()
        self.external_root = np.sqrt(self.external_y)
        denominator = node_sum-np.bincount(self.to_node,weights=2*self.y*self.direct/(1+self.direct),minlength=len(node_sum))
        denominator[self.external_nodes] += self.external_y
        if np.any(denominator <= 0): raise ValueError('Disconnected junction')
        channel = np.arange(count)
        a0=block_diag([csr_matrix(d[0]) for d in self.devices],format='csr')
        self.b0=block_diag([csr_matrix(d[1]) for d in self.devices],format='csr')
        self.c0=block_diag([csr_matrix(d[2]) for d in self.devices],format='csr')
        gather=coo_matrix((2*self.root_y/(1+self.direct)/denominator[self.to_node],(self.to_node,channel)),shape=(len(node_sum),count)).tocsr()
        self.px=gather@self.c0
        self.pu=coo_matrix((2*self.external_root/denominator[self.external_nodes],(self.external_nodes,np.arange(len(self.external_nodes)))),shape=(len(node_sum),len(self.external_nodes))).tocsr()
        emit=coo_matrix((self.root_y/(1+self.direct),(channel,self.from_node)),shape=(count,len(node_sum))).tocsr()
        subtract=coo_matrix((1/(1+self.direct),(channel,reverse)),shape=(count,count)).tocsr()
        self.A=(a0+self.b0@(emit@self.px-subtract@self.c0)).tocsr();self.A.eliminate_zeros()
        self.B=(self.b0@emit@self.pu).toarray()
        self.C=(diags(self.external_root)@self.px[self.external_nodes]).toarray()
        self.D=(diags(self.external_root)@self.pu[self.external_nodes]).toarray()-np.eye(len(self.external_nodes))
        # Radiation-flow readout (before its discrete differentiator).
        weights=np.zeros(count); pressure_weights=np.zeros(len(node_sum))
        for j,port in enumerate(snapshot['ports']):
            if port[2]:
                i=self.paths+j
                pressure_weights[port[0]] += self.y[i]*(1-self.direct[i])/(1+self.direct[i])
                weights[i] = -2*self.root_y[i]/(1+self.direct[i])
        self.O=np.asarray(pressure_weights@self.px + weights@self.c0).reshape(1,-1)
        self.OU=np.asarray(pressure_weights@self.pu).reshape(1,-1)

    def project(self, basis):
        if not np.allclose(basis.T@basis,np.eye(basis.shape[1]),atol=1e-10):
            raise ValueError('Basis must be orthonormal in stored-energy coordinates')
        return FlowModel(basis.T@(self.A@basis),basis.T@self.B,self.C@basis,self.D,
                         self.O@basis,self.OU,self.external_root,self.fs)

    def full_flow(self):
        return FlowModel(self.A,self.B,self.C,self.D,self.O,self.OU,self.external_root,self.fs)


class FlowModel:
    def __init__(self,a,b,c,d,o,ou,root_y,fs):
        self.fs=fs
        feedback=np.linalg.solve(np.eye(len(d))-d,np.eye(len(d)))
        self.A=a+csr_matrix(b)@csr_matrix(feedback)@csr_matrix(c) if hasattr(a,'tocsr') else a+b@feedback@c
        self.B=b@feedback@np.diag(1/root_y)
        self.C=2*np.diag(1/root_y)@feedback@c
        self.D=2*np.diag(1/root_y)@feedback@np.diag(1/root_y)-np.diag(1/(root_y*root_y))
        self.O=o+ou@feedback@c
        self.OU=ou@feedback@np.diag(1/root_y)

    def responses(self,frequency,damping=40):
        """All port drives in one factorization; columns remain independent."""
        z=np.exp((damping+2j*np.pi*frequency)/self.fs)
        if hasattr(self.A,'tocsr'):
            state=spsolve(z*eye(self.A.shape[0],format='csc')-self.A.tocsc(),self.B.astype(complex))
        else:
            state=np.linalg.solve(z*np.eye(len(self.A))-self.A,self.B)
        return np.vstack([self.C@state+self.D,self.fs*(1-1/z)*(self.O@state+self.OU)])

    def response(self,frequency,source=0,damping=40):
        z=np.exp((damping+2j*np.pi*frequency)/self.fs)
        if hasattr(self.A,'tocsr'):
            state=spsolve(z*eye(self.A.shape[0],format='csc')-self.A.tocsc(),self.B[:,source])
        else:
            state=np.linalg.solve(z*np.eye(len(self.A))-self.A,self.B[:,source])
        pressures=self.C@state+self.D[:,source]
        flow=(self.O@state+self.OU[:,source]).item()
        return np.append(pressures,self.fs*(1-1/z)*flow)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('snapshot',type=Path)
    p.add_argument('--out',type=Path,required=True)
    args=p.parse_args()
    snapshot=json.loads(args.snapshot.read_text())
    started=time.perf_counter();model=PassiveSurface(snapshot);prepare=time.perf_counter()-started
    print('Prepared',model.states,'energy states;',model.A.nnz,'nonzeros;',prepare,'seconds',flush=True)
    full=model.full_flow();reference=FrozenPorts(snapshot)
    rng=np.random.default_rng(9182)
    energies=[]
    for unused in range(20):
        x=rng.normal(size=model.states);u=rng.normal(size=len(model.external_y))
        next_x=model.A@x+model.B@u;out=model.C@x+model.D@u
        energies.append(float((next_x@next_x+out@out)/(x@x+u@u)))
    errors=[]
    for f in [48,100,500,1000,2000,4000,8000,12000]:
        actual=full.responses(f)
        expected=np.column_stack([reference.response(f,40,source_piece=int(piece))[0] for piece in model.external_nodes])
        error=float(np.max(abs(actual-expected)/np.maximum(abs(expected),1e-12)))
        errors.append(error); print(f,'Hz relative error',error,flush=True)
    report={'fixture':snapshot['fixture'],'input_ports_tested':len(model.external_nodes),'states':model.states,'matrix_nonzeros':model.A.nnz,
            'prepare_seconds':prepare,'maximum_component_step_gain':model.component_gain,
            'maximum_component_transfer_error':model.component_filter_error,
            'maximum_random_energy_ratio':max(energies),'max_relative_complex_error':max(errors)}
    args.out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
    if max(errors)>1e-5 or max(energies)>1+1e-10:raise SystemExit('Energy realization verification failed')

if __name__=='__main__':main()
