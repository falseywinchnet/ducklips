#!/usr/bin/env python3
"""Project frozen acoustic energy states; measure accuracy before accepting order."""
import argparse
import json
import hashlib
import time
from pathlib import Path
import numpy as np
from scipy.linalg import qr
from scipy.sparse import eye
from scipy.sparse.linalg import splu
from passive_surface import PassiveSurface
from surface_ports import FrozenPorts


def reduce(model, frequencies, rank_tolerance=1e-10, anatomical=False):
    if anatomical and len(model.snapshot.get('piece_tubes',[]))!=model.snapshot['pieces']:
        raise ValueError('Anatomical projection requires a fresh snapshot with piece_tubes metadata')
    flow = model.full_flow()
    columns=[]
    identity=eye(model.states,format='csc')
    for frequency in frequencies:
        z=np.exp((40+2j*np.pi*frequency)/model.fs)
        states=splu(z*identity-flow.A.tocsc()).solve(flow.B.astype(complex))
        for vector in np.hstack([states.real,states.imag]).T:
            norm=np.linalg.norm(vector)
            if norm>1e-20:columns.append(vector/norm)
    snapshots=np.column_stack(columns)
    groups=np.unique(model.state_groups) if anatomical else [None]
    blocks=[]
    for group in groups:
        rows=np.arange(model.states) if group is None else np.flatnonzero(model.state_groups==group)
        q,r,unused=qr(snapshots[rows],mode='economic',pivoting=True)
        diagonal=abs(np.diag(r))
        rank=int(np.sum(diagonal>diagonal[0]*rank_tolerance)) if diagonal[0]>0 else 0
        if rank:
            block=np.zeros((model.states,rank));block[rows]=q[:,:rank];blocks.append(block)
    basis=np.hstack(blocks)
    reduced=model.project(basis)
    ar=basis.T@(model.A@basis);br=basis.T@model.B;cr=model.C@basis
    step=np.block([[ar,br],[cr,model.D]])
    gain=np.linalg.svd(step,compute_uv=False)[0]
    radius=np.max(abs(np.linalg.eigvals(reduced.A)))
    return reduced,basis,float(gain),float(radius)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('snapshot',type=Path)
    parser.add_argument('--anchors',type=int,default=12)
    parser.add_argument('--validation-points',type=int,default=961)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--save-basis',action='store_true')
    parser.add_argument('--rank-tolerance',type=float,default=1e-10)
    parser.add_argument('--anatomical',action='store_true',help='Keep each cavity and connecting channel in a separate projection block')
    args=parser.parse_args()
    if args.validation_points<2 or args.anchors<2 or not np.isfinite(args.rank_tolerance) or not 0<args.rank_tolerance<1:
        parser.error('anchor and validation counts must be at least 2 and rank tolerance must be finite in (0,1)')
    snapshot=json.loads(args.snapshot.read_text())
    full=PassiveSurface(snapshot)
    frequencies=np.linspace(100,12000,args.anchors)
    started=time.perf_counter()
    reduced,basis,gain,radius=reduce(full,frequencies,args.rank_tolerance,args.anatomical)
    elapsed=time.perf_counter()-started
    print('Projected',full.states,'->',basis.shape[1],'states in',elapsed,'s',flush=True)
    reference=FrozenPorts(snapshot)
    test_hz=np.linspace(50,12000,args.validation_points)
    expected=np.array([[reference.response(f,0,source_piece=piece)[0]
                        for piece in snapshot['pressure_pieces']] for f in test_hz])
    actual=np.array([reduced.responses(f,0).T for f in test_hz])
    peak=np.max(abs(expected),axis=0)
    relative=abs(actual-expected)/np.maximum(abs(expected),np.maximum(peak*0.01,1e-30))
    channels=[]
    for source in range(actual.shape[1]):
        for i in range(actual.shape[2]):
            ex=expected[:,source,i];ac=actual[:,source,i]
            mask=abs(ex)>peak[source,i]*0.01
            channels.append({'input':source,'channel':i,
                'normalized_rms_error':float(np.linalg.norm(ac-ex)/max(np.linalg.norm(ex),1e-30)),
                'max_relative_complex_error_with_1pct_floor':float(relative[:,source,i].max()),
                'max_magnitude_error_db_above_1pct_peak':float(np.max(abs(20*np.log10(np.maximum(abs(ac[mask]),1e-30)/abs(ex[mask]))))) if mask.any() else 0.0,
                'max_phase_error_rad_above_1pct_peak':float(np.max(abs(np.angle(ac[mask]/ex[mask])))) if mask.any() else 0.0})
    report={'fixture':snapshot['fixture'],'snapshot_sha256':hashlib.sha256(args.snapshot.read_bytes()).hexdigest(),'full_states':full.states,'reduced_states':basis.shape[1],
            'anchor_frequencies_hz':frequencies.tolist(),'rank_tolerance':args.rank_tolerance,
            'anatomical_blocks':args.anatomical,'validation_frequency_count':len(test_hz),
            'validation_frequency_range_hz':[float(test_hz.min()),float(test_hz.max())],'preparation_seconds':elapsed,
            'scattering_step_largest_singular_value':gain,'flow_step_spectral_radius':radius,
            'channels':channels,'accepted_as_speech_renderer':False}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    extra={'basis':basis} if args.save_basis else {}
    np.savez_compressed(args.out.with_suffix('.npz'),A=reduced.A,B=reduced.B,C=reduced.C,D=reduced.D,
                        O=reduced.O,OU=reduced.OU,fs=full.fs,hz=test_hz,reference=expected,reduced=actual,**extra)
    print(json.dumps(report,indent=2),flush=True)
    if gain>1+1e-9 or radius>1+1e-9:raise SystemExit('Passivity/stability gate failed')

if __name__=='__main__':main()
