#!/usr/bin/env python3
"""Diagonalize a verified frozen reduced operator for a linear-cost runtime."""
import argparse
import json
import hashlib
from pathlib import Path
import numpy as np
from scipy.linalg import eig


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('model',type=Path);p.add_argument('--out',type=Path,required=True)
    args=p.parse_args()
    with np.load(args.model) as data:
        a,b,c,d,o,ou=[data[k] for k in ['A','B','C','D','O','OU']]
        fs=float(data['fs'])
    if not all(np.isfinite(value).all() for value in [a,b,c,d,o,ou]) or not np.isfinite(fs) or fs<=24000:
        raise ValueError('Expected finite coefficients and sample rate above 24 kHz')
    outputs=np.vstack([c,o]);direct=np.vstack([d,ou])
    poles,vectors=eig(a)
    keep=np.flatnonzero(poles.imag>=0)
    selected=poles[keep]
    # Solve for coordinates in a real eigenbasis. A complex solve of the full
    # conjugate eigenbasis leaves spurious imaginary roundoff on real poles.
    # x = 2 Re(v*z) for a pair, and x = v*r for a real mode.
    columns=[]
    for index in keep:
        if poles[index].imag>0:
            columns.extend([2*vectors[:,index].real,-2*vectors[:,index].imag])
        else:columns.append(vectors[:,index].real)
    real_basis=np.column_stack(columns)
    coordinates=np.linalg.solve(real_basis,b)
    drive=np.zeros((len(keep),b.shape[1]),dtype=complex)
    observe=np.zeros((outputs.shape[0],len(keep)),dtype=complex)
    cursor=0
    for i,index in enumerate(keep):
        drive[i]=coordinates[cursor];cursor+=1
        observe[:,i]=outputs@vectors[:,index]
        if poles[index].imag>0:
            drive[i]+=1j*coordinates[cursor];cursor+=1
            observe[:,i]*=2
    if cursor!=len(a):raise ValueError('Unpaired eigenvalues in a real model')
    expected_responses=[]
    actual_responses=[]
    worst=0.0
    worst_location=None
    for f in np.linspace(50,12000,241):
        z=np.exp(2j*np.pi*f/fs)
        expected=outputs@np.linalg.solve(z*np.eye(len(a))-a,b)+direct
        actual=direct.astype(complex).copy()
        for i,pole in enumerate(selected):
            if pole.imag>0:
                actual += np.outer(observe[:,i]/2,drive[i])/(z-pole)
                actual += np.outer(observe[:,i].conj()/2,drive[i].conj())/(z-pole.conjugate())
            else:actual += np.outer(observe[:,i],drive[i])/(z-pole)
        scale=np.maximum(abs(expected),1e-30)
        ratio=abs(actual-expected)/scale
        if float(ratio.max())>worst:
            worst=float(ratio.max())
            channel,source=np.unravel_index(np.argmax(ratio),ratio.shape)
            worst_location={'frequency_hz':float(f),'output':int(channel),'input':int(source),
                            'absolute_error':float(abs(actual-expected)[channel,source])}
        expected_responses.append(expected);actual_responses.append(actual)
    expected_responses=np.array(expected_responses);actual_responses=np.array(actual_responses)
    # Outputs mix pressure and flow units. Use each transfer's own frequency
    # range as the floor; never a louder channel or a different physical unit.
    peaks=abs(expected_responses).max(axis=0)
    scaled=abs(actual_responses-expected_responses)/np.maximum(abs(expected_responses),np.maximum(peaks*.01,1e-30))
    channel_error=scaled.max(axis=0)
    bounded_error=float(channel_error.max())
    if not np.isfinite(bounded_error) or bounded_error>1e-5:raise ValueError(f'Modal-coordinate roundoff too large: {bounded_error}')
    args.out.parent.mkdir(parents=True,exist_ok=True)
    with args.out.open('w') as out:
        out.write(f'DUCKLIPS_FROZEN_MODES_V1 {fs:.17g} {len(keep)} {b.shape[1]} {outputs.shape[0]}\n')
        np.savetxt(out,direct,fmt='%.17g')
        for index,pole in enumerate(selected):
            values=[pole.real,pole.imag]
            for value in drive[index]:values.extend([value.real,value.imag])
            for value in observe[:,index]:values.extend([value.real,value.imag])
            out.write(' '.join(f'{v:.17g}' for v in values)+'\n')
    report={'model_sha256':hashlib.sha256(args.model.read_bytes()).hexdigest(),'real_state_order':len(a),'stored_modal_states':len(keep),
            'sample_rate':fs,'eigenvector_condition_number':float(np.linalg.cond(vectors)),
            'real_coordinate_condition_number':float(np.linalg.cond(real_basis)),
            'max_modal_coordinate_relative_transfer_error_unfloored':worst,
            'worst_unfloored_location':worst_location,
            'max_modal_error_with_1pct_per_transfer_floor':bounded_error,
            'per_transfer_modal_error_with_1pct_floor':channel_error.tolist(),
            'maximum_pole_radius':float(np.max(abs(selected)))}
    Path(str(args.out)+'.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
