#!/usr/bin/env python3
"""Verify the native conservative material map on the previous Quebec fixture.
This is a kinematic geometry reference, not a muscle model or speech renderer.
"""
import argparse
import hashlib
import io
import json
import subprocess
from pathlib import Path
import numpy as np
from tongue_tissue_reference import ceiling, body_airway


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def native(args, points, command, steps, anchor, free, centre):
    output = subprocess.run([str(args.probe.resolve()), str((args.out/'profile.csv').resolve()),
                             str(anchor), str(free), str(centre), str(command), str(steps)],
                            input=''.join(f'{x:.17g},{y:.17g}\n' for x,y in points),
                            text=True, capture_output=True, check=True).stdout
    result = np.loadtxt(io.StringIO(output), delimiter=',', skiprows=1, ndmin=2)
    if len(result) != len(points) or not np.isfinite(result).all() or not np.all(result[:,8] == 1):
        raise RuntimeError('Native material mapping failed')
    return result[:,2:4], result[:,4:8].reshape(-1,2,2)


def polygon_area(points):
    return float(.5*np.sum(points[:,0]*np.roll(points[:,1],-1)-points[:,1]*np.roll(points[:,0],-1)))


def main(args):
    args.out.mkdir(parents=True,exist_ok=True)
    fixture=json.loads(args.fixture.read_text())
    centre=np.array(fixture['centre_cm']);radius=fixture['radius_cm'];target=fixture['wall_target_cm'][0]
    half=.45;shoulder=1.5
    if not np.allclose(fixture['normal'],[0,1]) or fixture['gap_cm']!=0:
        raise ValueError('This fixture tests the upward zero-gap oral closure')
    endpoints=np.array([target-half,target+half])
    anchor=centre[1]-radius/2
    free=float(np.min(centre[1]+np.sqrt(radius**2-(endpoints-centre[0])**2)))
    x=np.unique(np.r_[np.linspace(target-half-shoulder,target+half+shoulder,257),endpoints,target,2.,4.5])
    x=x[(x>=target-half-shoulder)&(x<=target+half+shoulder)]
    if np.any(abs(x-centre[0])>=radius):
        raise ValueError('Lift support extends beyond the resting body')
    radial=np.sqrt(radius**2-(x-centre[0])**2)
    top=centre[1]+radial;roof,roof_slope=ceiling(x)
    q=np.clip((abs(x-target)-half)/shoulder,0,1)
    weight=1-q**3*(10+q*(-15+6*q))
    weight_slope=-30*q*q*(1-q)**2*np.sign(x-target)/shoulder
    lift=(roof-top)*weight
    slope=(roof_slope+(x-centre[0])/radial)*weight+(roof-top)*weight_slope
    lift[[0,-1]]=0;slope[[0,-1]]=0
    profile=np.c_[x,lift,slope]
    np.savetxt(args.out/'profile.csv',profile,delimiter=',',fmt='%.17g')
    n=8192
    angle=np.arange(n)*2*np.pi/n
    boundary=centre+radius*np.c_[np.cos(angle),np.sin(angle)]
    rng=np.random.default_rng(891)
    theta=rng.uniform(0,2*np.pi,4096);r=radius*np.sqrt(rng.uniform(0,1,4096))
    interior=centre+r[:,None]*np.c_[np.cos(theta),np.sin(theta)]
    patch_x=np.linspace(target-half,target+half,129)
    patch=np.c_[patch_x,centre[1]+np.sqrt(radius**2-(patch_x-centre[0])**2)]
    points=np.vstack([boundary,interior,patch])
    records=[];frames=[]
    for command in np.linspace(0,1,11):
        mapped,gradient=native(args,points,float(command),64,anchor,free,target)
        determinant=np.linalg.det(gradient)
        singular=np.linalg.svd(gradient,compute_uv=False)
        fixed=points[:,1]<=anchor
        clearance=ceiling(mapped[:,0])[0]-mapped[:,1]
        patch_clearance=clearance[-len(patch):]
        record={'command':float(command),'local_volume_error':float(np.max(abs(determinant-1))),
                'maximum_principal_stretch':float(singular.max()),
                'anchor_displacement_cm':float(np.linalg.norm(mapped[fixed]-points[fixed],axis=1).max()),
                'wall_penetration_cm':float(max(0,-clearance.min(),np.max(-.9-mapped[:,0]))),
                'patch_clearance_cm_range':[float(patch_clearance.min()),float(patch_clearance.max())],
                'boundary_polygon_area_cm2':polygon_area(mapped[:n])}
        if record['local_volume_error']>1e-10 or record['anchor_displacement_cm']!=0 or record['wall_penetration_cm']>1e-7:
            raise RuntimeError(json.dumps(record))
        records.append(record);frames.append(mapped)
    closed=frames[-1]
    if max(abs(np.array(records[-1]['patch_clearance_cm_range'])))>1e-7:
        raise RuntimeError('Dorsal patch did not reach the roof')
    back,unused=native(args,closed,-1,64,anchor,free,target)
    inverse_error=float(np.max(np.linalg.norm(back-points,axis=1)))
    check=interior[:256];step=1e-5
    difference_points=np.vstack([check+[step,0],check-[step,0],check+[0,step],check-[0,step]])
    difference,unused=native(args,difference_points,1,64,anchor,free,target)
    nominal,analytic=native(args,check,1,64,anchor,free,target)
    finite=np.stack([(difference[:256]-difference[256:512])/(2*step),
                     (difference[512:768]-difference[768:])/(2*step)],axis=-1)
    gradient_error=float(np.max(abs(finite-analytic)))
    fd_volume_error=float(np.max(abs(np.linalg.det(finite)-1)))
    if inverse_error>1e-9 or gradient_error>2e-6 or fd_volume_error>2e-6:
        raise RuntimeError('Independent material map derivative or inverse check failed')
    coarse,unused=native(args,points,1,16,anchor,free,target)
    medium,unused=native(args,points,1,32,anchor,free,target)
    fine,unused=native(args,points,1,128,anchor,free,target)
    refinement={'16_to_32_max_cm':float(np.linalg.norm(coarse-medium,axis=1).max()),
                '32_to_64_max_cm':float(np.linalg.norm(medium-closed,axis=1).max()),
                '64_to_128_max_cm':float(np.linalg.norm(closed-fine,axis=1).max())}
    if not refinement['64_to_128_max_cm']<refinement['32_to_64_max_cm']<refinement['16_to_32_max_cm']:
        raise RuntimeError('Material trajectories did not converge under step refinement')
    contour=[]
    for stride in [16,4,1]:
        rest_area=polygon_area(boundary[::stride]);closed_area=polygon_area(closed[:n:stride])
        contour.append({'vertices':n//stride,'rest_cm2':rest_area,'closed_cm2':closed_area,
                        'relative_change':closed_area/rest_area-1})
    volume_rest,gap=body_airway(fixture,boundary)
    volume_closed,gap=body_airway(fixture,closed[:n])
    shift=np.array(fixture['poses'][-1]['centre_cm'])-centre
    volume_rigid,gap=body_airway(fixture,boundary+shift)
    np.savez(args.out/'material.npz',rest=points,closed=closed,gradient=gradient,boundary_count=n,half=frames[5])
    report={'baseline_commit':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
            'source_hashes':{str(p):sha(p) for p in ['include/ducklips/tongue_flow.hpp','src/tongue_flow.cpp',
                'tools/tongue_flow_probe.cpp','tests/tongue_flow_tests.cpp','tools/tongue_flow_reference.py','tools/tongue_tissue_reference.py']},
            'probe_sha256':sha(args.probe),'fixture_sha256':sha(args.fixture),'profile_sha256':sha(args.out/'profile.csv'),
            'fixture':{k:fixture[k] for k in ['centre_cm','radius_cm','place_cm','gap_cm','wall_target_cm']},
            'profile':{'knots':len(x),'patch_width_cm':2*half,'shoulder_cm':shoulder,'anchor_y_cm':anchor,'free_y_cm':free},
            'samples_per_command':len(points),'records':records,
            'independent_checks':{'inverse_recovery_max_cm':inverse_error,'gradient_max_absolute_error':gradient_error,
                                  'finite_difference_volume_max_error':fd_volume_error},
            'step_refinement':refinement,'contour_refinement':contour,
            'body_only_airway_ml':{'rest':volume_rest,'deformed':volume_closed,'rigid':volume_rigid},
            'scope':'Native kinematic material-map component; no calibrated elasticity, transverse anatomy, speech integration, pronunciation or performance claim',
            'material_sha256':sha(args.out/'material.npz')}
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(8,7))
    sections=np.array(fixture['sections']);ax.plot(sections[:,1],sections[:,2],color='black',label='Native outer wall')
    for p,label,color,style in [(boundary,'Rest','#888888','--'),(frames[5][:n],'Half command','#ae8f58','-'),
                                 (closed[:n],'Conservative map: full command','#ab492c','-'),
                                 (boundary+shift,'Current rigid gesture','#447ba9','--')]:
        loop=np.vstack([p,p[0]]);ax.plot(loop[:,0],loop[:,1],color=color,linestyle=style,label=label,lw=1.6)
    ax.scatter(closed[n:n+len(interior):8,0],closed[n:n+len(interior):8,1],s=3,color='#ab492c',alpha=.4)
    ax.axhline(anchor,color='#777777',linestyle=':',label='Fixed tissue below this height')
    ax.set(xlabel='Anterior position (cm)',ylabel='Height (cm)',aspect='equal',xlim=(-1.2,8.3),ylim=(2.2,10.5),
           title='Native conservative tongue map\nKinematic closure, not a speech render')
    ax.legend(loc='lower right',fontsize=8);fig.tight_layout();fig.savefig(args.out/'material-closure.png',dpi=160);plt.close(fig)
    (args.out/'receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'closure':records[-1],'independent_checks':report['independent_checks'],
                      'step_refinement':refinement,'contour_refinement':contour,
                      'body_only_airway_ml':report['body_only_airway_ml']},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fixture',type=Path)
    parser.add_argument('--probe',type=Path,default=Path('build/ducklips_tongue_flow_probe'))
    parser.add_argument('--out',type=Path,required=True)
    main(parser.parse_args())
