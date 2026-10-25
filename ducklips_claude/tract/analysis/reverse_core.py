"""The reversible core: the vocal tube as a physical two-port that runs both
ways (the owner, 2026-10-19: "a physical model that, when excited by sound,
replicates the mouth shapes ... end to end reversibility").

Forward, a replica of Tube::step_air for a rigid, lossless tube: travelling
pressure waves hop one node per step (kStepRate 96 kHz, nodes c / 96000 =
0.3646 cm apart), each node a junction of its two links (link area the mean
of the two nodes'), the glottis end driven by the flow entering the tube,
the lips loaded by the engine's radiation (R parallel L for the last node's
area). Backward, the same tube driven by the sound: layer peeling, node by
node from the lips to the glottis (peel). (Inverting the transfer by FFT
failed: the nearly lossless tube rings far past the signal's end, and the
cut-off ringing smeared the inverse over every band: a synthetic round trip
came back at 4-36 dB.)

usage: reverse_core.py AREAS RAW   (sustain's AREAS= and RAW= files; RAW
holds float32 pairs per step: the flow entering the vocal tube, the flow
radiated at the lips)"""
import os
import sys

import numpy
from scipy.signal import lfilter
from numba import njit

STEP_RATE = 96000.0
AIR_DENSITY = 1.14e-3      # g/cm^3
SOUND_SPEED = 35000.0      # cm/s
RHO_C = AIR_DENSITY * SOUND_SPEED
H = 1.0 / STEP_RATE
NODE_LENGTH = SOUND_SPEED / STEP_RATE


def viscous_per_cm(area):
    """The engine's viscous resistance per cm of a link (Tube's
    viscous_resistance_per_cm)."""
    safe = numpy.maximum(area, 1e-4)
    return numpy.maximum(0.03 * safe ** -2.6, 0.0065 / safe * numpy.sqrt(2250.0))


def radiation(area, end_extra=0.0):
    """The lips' load for an opening of `area` (cm^2): resistance and
    inertance (Tube::step_air, scales 1), the inertance with the air of the
    tube's fraction of a node past its last whole node (`end_extra` cm)."""
    resistance = 128.0 * RHO_C / (9.0 * numpy.pi ** 2 * area)
    inertance = 8.0 * AIR_DENSITY / (3.0 * numpy.pi * numpy.sqrt(numpy.pi * area)) + AIR_DENSITY * end_extra / area
    return resistance, inertance


def forward(areas, flow_in, kinetic=True, end_extra=0.0, viscous=0.0, hop=0.0):
    """The flow radiated at the lips for the flow entering the tube, step by
    step, as the engine's rigid tube does it: lossless hops, and on each link
    the kinetic loss of air expanding out of a narrowing (Borda-Carnot, the
    engine's series_element with no viscous resistance)."""
    areas = numpy.asarray(areas, float)
    n = len(areas)
    link = 0.5 * (areas[:-1] + areas[1:])
    y = link / RHO_C                       # each link's admittance
    f = numpy.zeros(n - 1)                 # waves on link i heading right (arriving at node i + 1)
    b = numpy.zeros(n - 1)                 # waves on link i heading left (arriving at node i)
    resistance, inertance = radiation(areas[-1], end_extra)
    gain = 0.5 * H / inertance
    inertance_flow = 0.0
    previous = 0.0
    out = numpy.zeros(len(flow_in))
    p = numpy.zeros(n)
    resistance_link = viscous * viscous_per_cm(link) * NODE_LENGTH
    f_smooth = numpy.zeros(n - 1)
    b_smooth = numpy.zeros(n - 1)
    for step, u in enumerate(flow_in):
        # Junction pressures from the waves arriving this step.
        p[0] = u / y[0] + 2.0 * b[0]
        p[1:-1] = (2.0 * y[:-1] * f[:-1] + 2.0 * y[1:] * b[1:]) / (y[:-1] + y[1:])
        admittance = y[-1] + 1.0 / resistance + gain
        p[-1] = (2.0 * y[-1] * f[-1] - inertance_flow - gain * previous) / admittance
        inertance_flow += gain * (p[-1] + previous)
        previous = p[-1]
        out[step] = p[-1] / resistance + inertance_flow
        # The waves leaving each junction arrive at its neighbours next step.
        # Each hop passes the tube's one-pole low-pass (wave_smoothing).
        f_smooth = (1.0 - hop) * (p[:-1] - b) + hop * f_smooth
        b_smooth = (1.0 - hop) * (p[1:] - f) + hop * b_smooth
        f, b = scatter(f_smooth, b_smooth, link, areas, resistance_link, kinetic)
    return out


def scatter(forward_waves, backward_waves, link, areas, resistance, kinetic=True):
    """Each link's series element: 2 (a - b) = (2 Z + R) U + K U|U|, R the
    viscous resistance, K the Borda-Carnot loss into the wider side the air
    flows to (Tube's series_element)."""
    z = RHO_C / link
    drive = 2.0 * (forward_waves - backward_waves)
    rightward = drive >= 0.0
    narrow = numpy.where(rightward, areas[:-1], areas[1:])
    wide = numpy.where(rightward, areas[1:], areas[:-1])
    difference = 1.0 / numpy.maximum(narrow, 1e-4) - 1.0 / wide
    k = numpy.where(wide > narrow, 0.5 * AIR_DENSITY * difference ** 2, 0.0) if kinetic else 0.0
    linear = 2.0 * z + resistance
    magnitude = 2.0 * numpy.abs(drive) / (linear + numpy.sqrt(linear ** 2 + 4.0 * k * numpy.abs(drive)))
    flow = numpy.where(rightward, magnitude, -magnitude)
    return backward_waves + z * flow, forward_waves - z * flow


def delayed(signal, steps=1):
    """The signal `steps` steps later (zeros before it began)."""
    out = numpy.zeros_like(signal)
    out[steps:] = signal[:-steps]
    return out


def advanced(signal, steps=1):
    """The signal `steps` steps earlier (zeros past its end)."""
    out = numpy.zeros_like(signal)
    out[:-steps] = signal[steps:]
    return out


def peel(areas, flow_out, kinetic=True, end_extra=0.0, viscous=0.0, hop=0.0):
    """The tube driven backward by the sound: the flow radiated at the lips
    carried node by node to the flow entering at the glottis, by the same
    junction, link and load equations as forward, solved the other way.
    Each node inward costs one step: the flow entering at step t is known
    from the radiated flow up to step t + (nodes - 1).

    At the lips the load gives the pressure step by step (radiated = p / R
    + inertance flow); the junction there gives the wave that arrived. On
    each link the arriving (post-link) wave and the wave that left the far
    node fix the link's flow, U = (f' - b) / Z, and from it the wave that
    entered the link, a = b + Z U + K U|U| / 2 (the kinetic loss made
    explicit); that wave left the near node one step earlier, which gives
    the near node's pressure; its junction gives the next link's arriving
    wave."""
    areas = numpy.asarray(areas, float)
    n = len(areas)
    link = 0.5 * (areas[:-1] + areas[1:])
    y = link / RHO_C
    z = RHO_C / link
    count = len(flow_out)
    # The lips: pressure from the load, step by step.
    resistance, inertance = radiation(areas[-1], end_extra)
    gain = 0.5 * H / inertance
    pressure = numpy.zeros(count)
    inertance_flow = numpy.zeros(count)
    flow_l = 0.0
    previous = 0.0
    for step in range(count):
        pressure[step] = (flow_out[step] - flow_l - gain * previous) / (1.0 / resistance + gain)
        flow_l += gain * (pressure[step] + previous)
        inertance_flow[step] = flow_l
        previous = pressure[step]
    admittance = y[-1] + 1.0 / resistance + gain
    arriving = (pressure * admittance + delayed(inertance_flow) + gain * delayed(pressure)) / (2.0 * y[-1])
    resistance_link = viscous * viscous_per_cm(link) * NODE_LENGTH
    for j in range(n - 2, -1, -1):
        # Link j: `arriving` is its rightward wave reaching node j + 1; the
        # leftward wave entering it left node j + 1 one step earlier, through
        # the hop's low-pass.
        leftward_in = lfilter([1.0 - hop], [1.0, -hop], delayed(pressure) - delayed(arriving))
        flow = (arriving - leftward_in) / z[j]
        rightward_in = leftward_in + (z[j] + 0.5 * resistance_link[j]) * flow
        if kinetic:
            narrow = numpy.where(flow >= 0.0, areas[j], areas[j + 1])
            wide = numpy.where(flow >= 0.0, areas[j + 1], areas[j])
            difference = 1.0 / numpy.maximum(narrow, 1e-4) - 1.0 / wide
            k = numpy.where(wide > narrow, 0.5 * AIR_DENSITY * difference ** 2, 0.0)
            rightward_in = rightward_in + 0.5 * k * flow * numpy.abs(flow)
        leftward_out = rightward_in - z[j] * flow
        # The rightward wave entering link j left node j one step earlier,
        # before the hop's low-pass (undone): node j's pressure is it plus
        # the leftward wave that node received.
        rightward_raw = (rightward_in - hop * delayed(rightward_in)) / (1.0 - hop)
        pressure = advanced(rightward_raw) + leftward_out
        if j > 0:
            arriving = (pressure * (y[j - 1] + y[j]) - 2.0 * y[j] * leftward_out) / (2.0 * y[j - 1])
        else:
            return y[0] * (pressure - 2.0 * leftward_out)
    return None


# ---- The tube with its walls (each node's lining local: no axial tension) ----
# walls: per node surface, cover mass, body mass, hoop, coupling and support
# stiffness, cover and body damping, rest radius (sustain's WALLDUMP); a node
# whose cover mass is 0 is rigid.


@njit(cache=True)
def _lining_springs(walls, q, v, bq, bv):
    """Tube::step_lining for local linings: each cover on its hoop and
    coupling springs (implicit midpoint), then each body on coupling and
    support (backward Euler), as one step."""
    h = 1.0 / 96000.0
    quarter = 0.25 * h * h
    for i in range(walls.shape[0]):
        mc = walls[i, 1]
        if mc <= 0.0:
            continue
        mb = walls[i, 2]
        hoop = walls[i, 3]
        coupling = walls[i, 4]
        support = walls[i, 5]
        cc = walls[i, 6]
        bc = walls[i, 7]
        diag = mc + quarter * (hoop + coupling)
        rhs = mc * (2.0 * q[i] + h * v[i]) + 0.5 * h * h * coupling * bq[i]
        new_q = rhs / diag - q[i]
        v[i] = 2.0 * (new_q - q[i]) / h - v[i]
        q[i] = new_q
        stiffness = coupling + support
        pull = coupling * q[i] - stiffness * bq[i] + cc * v[i]
        bv[i] = (mb * bv[i] + h * pull) / (mb + h * (bc + cc) + h * h * stiffness)
        bq[i] += h * bv[i]


@njit(cache=True)
def _junction(i, walls, known, admittance, v, bv):
    """A node's junction with its lining (Tube::step_air): `known` is the
    waves' twice-sum (or its flow-source equivalent), `admittance` the
    node's; the cover velocity and the pressure found together. Returns the
    pressure (v updated)."""
    h = 1.0 / 96000.0
    mc = walls[i, 1]
    if mc <= 0.0:
        return known / admittance
    surface = walls[i, 0]
    cc = walls[i, 6]
    resistance = surface * surface / admittance + cc
    drive = surface * known / admittance + cc * bv[i]
    following = ((mc - 0.5 * h * resistance) * v[i] + h * drive) / (mc + 0.5 * h * resistance)
    middle = 0.5 * (v[i] + following)
    v[i] = following
    return (known - surface * middle) / admittance


@njit(cache=True)
def _wall_middle(i, walls, pressure, v, bv):
    """The same junction solved the other way: with the node's pressure
    known, the cover velocity's step and its mean over the step (v updated)."""
    h = 1.0 / 96000.0
    mc = walls[i, 1]
    if mc <= 0.0:
        return 0.0
    surface = walls[i, 0]
    cc = walls[i, 6]
    following = (mc * v[i] + h * surface * pressure - h * cc * (0.5 * v[i] - bv[i])) / (mc + 0.5 * h * cc)
    middle = 0.5 * (v[i] + following)
    v[i] = following
    return middle


@njit(cache=True)
def _series(forward_wave, backward_wave, z, resistance, kinetic, first_area, second_area, source=0.0):
    """One link's series element; returns (flow, to_second, to_first)."""
    drive = 2.0 * (forward_wave - backward_wave) + source
    if drive >= 0.0:
        narrow, wide = first_area, second_area
    else:
        narrow, wide = second_area, first_area
    k = 0.0
    if kinetic and wide > narrow:
        difference = 1.0 / max(narrow, 1e-4) - 1.0 / wide
        k = 0.5 * 1.14e-3 * difference * difference
    linear = 2.0 * z + resistance
    magnitude = 2.0 * abs(drive) / (linear + numpy.sqrt(linear * linear + 4.0 * k * abs(drive)))
    flow = magnitude if drive >= 0.0 else -magnitude
    return flow, backward_wave + z * flow, forward_wave - z * flow


@njit(cache=True)
def _forward_tube(areas, walls, flow_in, resistance_link, hop, kinetic, end_extra):
    n = areas.shape[0]
    rho_c = 1.14e-3 * 35000.0
    h = 1.0 / 96000.0
    link = 0.5 * (areas[:-1] + areas[1:])
    y = link / rho_c
    z = rho_c / link
    f = numpy.zeros(n - 1)
    b = numpy.zeros(n - 1)
    f_smooth = numpy.zeros(n - 1)
    b_smooth = numpy.zeros(n - 1)
    q = numpy.zeros(n)
    v = numpy.zeros(n)
    bq = numpy.zeros(n)
    bv = numpy.zeros(n)
    p = numpy.zeros(n)
    resistance = 128.0 * rho_c / (9.0 * numpy.pi ** 2 * areas[-1])
    inertance = 8.0 * 1.14e-3 / (3.0 * numpy.pi * numpy.sqrt(numpy.pi * areas[-1])) + 1.14e-3 * end_extra / areas[-1]
    gain = 0.5 * h / inertance
    inertance_flow = 0.0
    previous = 0.0
    out = numpy.zeros(flow_in.shape[0])
    for step in range(flow_in.shape[0]):
        _lining_springs(walls, q, v, bq, bv)
        # Node 0, driven by the flow entering: p y0 = U + 2 y0 b0 - s middle.
        p[0] = _junction(0, walls, flow_in[step] + 2.0 * y[0] * b[0], y[0], v, bv)
        for i in range(1, n - 1):
            p[i] = _junction(i, walls, 2.0 * y[i - 1] * f[i - 1] + 2.0 * y[i] * b[i], y[i - 1] + y[i], v, bv)
        admittance = y[n - 2] + 1.0 / resistance + gain
        p[n - 1] = _junction(n - 1, walls, 2.0 * y[n - 2] * f[n - 2] - inertance_flow - gain * previous, admittance, v, bv)
        inertance_flow += gain * (p[n - 1] + previous)
        previous = p[n - 1]
        out[step] = p[n - 1] / resistance + inertance_flow
        for i in range(n - 1):
            f_smooth[i] = (1.0 - hop) * (p[i] - b[i]) + hop * f_smooth[i]
            b_smooth[i] = (1.0 - hop) * (p[i + 1] - f[i]) + hop * b_smooth[i]
        for i in range(n - 1):
            flow, f[i], b[i] = _series(f_smooth[i], b_smooth[i], z[i], resistance_link[i], kinetic, areas[i], areas[i + 1])
    return out


@njit(cache=True)
def _peel_tube(areas, walls, flow_out, resistance_link, hop, kinetic, end_extra):
    n = areas.shape[0]
    count = flow_out.shape[0]
    rho_c = 1.14e-3 * 35000.0
    h = 1.0 / 96000.0
    link = 0.5 * (areas[:-1] + areas[1:])
    y = link / rho_c
    z = rho_c / link
    q = numpy.zeros(n)
    v = numpy.zeros(n)
    bq = numpy.zeros(n)
    bv = numpy.zeros(n)
    # The lips: the pressure from the load, the lining with it, the wave that
    # arrived.
    resistance = 128.0 * rho_c / (9.0 * numpy.pi ** 2 * areas[-1])
    inertance = 8.0 * 1.14e-3 / (3.0 * numpy.pi * numpy.sqrt(numpy.pi * areas[-1])) + 1.14e-3 * end_extra / areas[-1]
    gain = 0.5 * h / inertance
    admittance = y[n - 2] + 1.0 / resistance + gain
    pressure = numpy.zeros(count)
    arriving = numpy.zeros(count)
    flow_l = 0.0
    previous = 0.0
    for t in range(count):
        pressure[t] = (flow_out[t] - flow_l - gain * previous) / (1.0 / resistance + gain)
        _lining_one(n - 1, walls, q, v, bq, bv)
        middle = _wall_middle(n - 1, walls, pressure[t], v, bv)
        known = pressure[t] * admittance + walls[n - 1, 0] * middle
        arriving[t] = (known + flow_l + gain * previous) / (2.0 * y[n - 2])
        flow_l += gain * (pressure[t] + previous)
        previous = pressure[t]
    leftward_in = numpy.zeros(count)
    rightward_in = numpy.zeros(count)
    leftward_out = numpy.zeros(count)
    for j in range(n - 2, -1, -1):
        # Link j: the leftward wave entering it left node j + 1 a step ago,
        # through the hop's low-pass; with the arriving rightward wave it
        # fixes the link's flow and so the wave that entered from node j.
        smooth = 0.0
        for t in range(count):
            raw = (pressure[t - 1] - arriving[t - 1]) if t > 0 else 0.0
            smooth = (1.0 - hop) * raw + hop * smooth
            leftward_in[t] = smooth
            flow = (arriving[t] - smooth) / z[j]
            entered = smooth + (z[j] + 0.5 * resistance_link[j]) * flow
            if kinetic:
                if flow >= 0.0:
                    narrow, wide = areas[j], areas[j + 1]
                else:
                    narrow, wide = areas[j + 1], areas[j]
                if wide > narrow:
                    difference = 1.0 / max(narrow, 1e-4) - 1.0 / wide
                    entered += 0.5 * (0.5 * 1.14e-3 * difference * difference) * flow * abs(flow)
            rightward_in[t] = entered
            leftward_out[t] = entered - z[j] * flow
        # Node j: the rightward wave entering link j left it a step earlier,
        # before the hop's low-pass (undone); its pressure, its lining, and
        # the wave that arrived from the left (or, at node 0, the flow in).
        for t in range(count):
            following = rightward_in[t + 1] if t + 1 < count else 0.0
            current = rightward_in[t]
            raw = (following - hop * current) / (1.0 - hop)
            pressure[t] = raw + leftward_out[t]
        for t in range(count):
            _lining_one(j, walls, q, v, bq, bv)
            middle = _wall_middle(j, walls, pressure[t], v, bv)
            if j > 0:
                known = pressure[t] * (y[j - 1] + y[j]) + walls[j, 0] * middle
                arriving[t] = (known - 2.0 * y[j] * leftward_out[t]) / (2.0 * y[j - 1])
            else:
                arriving[t] = pressure[t] * y[0] + walls[0, 0] * middle - 2.0 * y[0] * leftward_out[t]
    return arriving, pressure


@njit(cache=True)
def _lining_one(i, walls, q, v, bq, bv):
    """_lining_springs for one node."""
    h = 1.0 / 96000.0
    quarter = 0.25 * h * h
    mc = walls[i, 1]
    if mc <= 0.0:
        return
    mb = walls[i, 2]
    hoop = walls[i, 3]
    coupling = walls[i, 4]
    support = walls[i, 5]
    cc = walls[i, 6]
    bc = walls[i, 7]
    diag = mc + quarter * (hoop + coupling)
    rhs = mc * (2.0 * q[i] + h * v[i]) + 0.5 * h * h * coupling * bq[i]
    new_q = rhs / diag - q[i]
    v[i] = 2.0 * (new_q - q[i]) / h - v[i]
    q[i] = new_q
    stiffness = coupling + support
    pull = coupling * q[i] - stiffness * bq[i] + cc * v[i]
    bv[i] = (mb * bv[i] + h * pull) / (mb + h * (bc + cc) + h * h * stiffness)
    bq[i] += h * bv[i]


# The engine's wall materials (Tube: soft_wall, firm_wall, lip_wall): cover
# density and thickness, hoop tension, cover loss, coupling stiffness, body
# density and depth, support stiffness and loss.
MATERIALS = {
    "soft": (1.04, 0.10, 3000.0, 0.15, 4.0e8, 1.04, 1.5, 1.0e5, 1.0),
    "firm": (1.04, 0.10, 6000.0, 0.15, 1.0e9, 1.04, 1.5, 4.0e5, 1.0),
    "lip": (1.04, 0.15, 2000.0, 0.15, 3.0e8, 1.04, 1.2, 8.0e4, 1.0),
}


def lining_constants(areas):
    """Each node's lining as the engine compiles it (Tube::compile_node, the
    vowel layout's materials: firm within 2 cm of the glottis, lip in the
    last cm, soft between; support factor 1): the WALLDUMP columns."""
    areas = numpy.asarray(areas, float)
    count = len(areas)
    total = count * NODE_LENGTH
    out = numpy.zeros((count, 9))
    for node in range(count):
        x = (node + 0.5) * NODE_LENGTH
        name = "firm" if x < 2.0 else ("lip" if x > total - 1.0 else "soft")
        density, thickness, hoop_tension, cover_loss, coupling, body_density, depth, support, support_loss = MATERIALS[name]
        radius = numpy.sqrt(max(areas[node], 1e-4) / numpy.pi)
        surface = 2.0 * numpy.pi * radius * NODE_LENGTH
        cover_mass = density * thickness * surface
        body_mass = body_density * depth * surface
        hoop_k = hoop_tension * surface / (radius * radius)
        coupling_k = coupling * surface
        support_k = support * surface
        reduced = cover_mass * body_mass / (cover_mass + body_mass)
        cover_c = 2.0 * cover_loss * numpy.sqrt(reduced * (hoop_k + coupling_k))
        body_c = 2.0 * support_loss * numpy.sqrt(body_mass * support_k)
        out[node] = (surface, cover_mass, body_mass, hoop_k, coupling_k, support_k, cover_c, body_c, radius)
    return out


def forward_tube(areas, flow_in, walls=None, viscous=0.0, hop=0.0, kinetic=True, end_extra=0.0):
    """The flow radiated for the flow entering, with the walls' linings."""
    areas = numpy.asarray(areas, float)
    if walls is None:
        walls = numpy.zeros((len(areas), 9))
    link = 0.5 * (areas[:-1] + areas[1:])
    return _forward_tube(areas, numpy.asarray(walls, float), numpy.asarray(flow_in, float),
                         viscous * viscous_per_cm(link) * NODE_LENGTH, float(hop), bool(kinetic), float(end_extra))


def peel_tube(areas, flow_out, walls=None, viscous=0.0, hop=0.0, kinetic=True, end_extra=0.0):
    """The tube with its walls driven backward by the sound (as peel)."""
    areas = numpy.asarray(areas, float)
    if walls is None:
        walls = numpy.zeros((len(areas), 9))
    link = 0.5 * (areas[:-1] + areas[1:])
    return peel_tube_glottis(areas, flow_out, walls, viscous, hop, kinetic, end_extra)[0]


def peel_tube_glottis(areas, flow_out, walls=None, viscous=0.0, hop=0.0, kinetic=True, end_extra=0.0):
    """peel_tube's flow entering and the first node's pressure."""
    areas = numpy.asarray(areas, float)
    if walls is None:
        walls = numpy.zeros((len(areas), 9))
    link = 0.5 * (areas[:-1] + areas[1:])
    return _peel_tube(areas, numpy.asarray(walls, float), numpy.asarray(flow_out, float),
                      viscous * viscous_per_cm(link) * NODE_LENGTH, float(hop), bool(kinetic), float(end_extra))


# A moving tube: per step and node, sustain's GEODUMP row: area, the air
# the frame draws from the node, its link's viscous resistance, and the
# lining (WALLDUMP's nine). Tube::step steps the linings on the last step's
# constants, then recompiles the geometry, then the air on the new one; the
# links scatter on the geometry of the step that made their waves.
GEOMETRY_COLUMNS = 12


def read_geometry(path, count):
    """sustain's GEODUMP as [step, node, column]."""
    return numpy.fromfile(path, dtype=numpy.float64).reshape(-1, count, GEOMETRY_COLUMNS)


@njit(cache=True)
def _radiation_gain(area, end_extra):
    return _radiation_gain_scaled(area, end_extra, 1.0, 1.0)


@njit(cache=True)
def _radiation_gain_scaled(area, end_extra, mass_scale, resistance_scale):
    """The open end's load (Tube's, with its radiation scales)."""
    rho_c = 1.14e-3 * 35000.0
    resistance = resistance_scale * 128.0 * rho_c / (9.0 * numpy.pi ** 2 * area)
    inertance = mass_scale * 8.0 * 1.14e-3 / (3.0 * numpy.pi * numpy.sqrt(numpy.pi * area)) + 1.14e-3 * end_extra / area
    return resistance, 0.5 / 96000.0 / inertance


@njit(cache=True)
def _kinetic_term(flow, first_area, second_area):
    """The Borda-Carnot half-term K U|U| / 2 of a link (as _series)."""
    if flow >= 0.0:
        narrow, wide = first_area, second_area
    else:
        narrow, wide = second_area, first_area
    if wide <= narrow:
        return 0.0
    difference = 1.0 / max(narrow, 1e-4) - 1.0 / wide
    return 0.5 * (0.5 * 1.14e-3 * difference * difference) * flow * abs(flow)


# The reversible jet (VoiceBody::reversible_jet), its constants written in
# as the layout's defaults (first oral node 10, Strouhal 0.2, frication band
# 3, corner ceiling 12 kHz, air viscosity 1.86e-4 P); the layout's seed.
JET_SEED = 73021
# The links' viscous glide in the engine's phase (_moving_geometry); env
# ENGINE_GLIDE_PHASE=0 restores the replica as it was before 2026-10-20.
ENGINE_GLIDE_PHASE = int(os.environ.get("ENGINE_GLIDE_PHASE", "1"))


@njit(cache=True)
def _split_mix(key):
    key = key + numpy.uint64(0x9E3779B97F4A7C15)
    key = (key ^ (key >> numpy.uint64(30))) * numpy.uint64(0xBF58476D1CE4E5B9)
    key = (key ^ (key >> numpy.uint64(27))) * numpy.uint64(0x94D049BB133111EB)
    return key ^ (key >> numpy.uint64(31))


@njit(cache=True)
def _hashed_gauss(seed, link, step):
    total = 0.0
    for draw in range(4):
        key = numpy.uint64(seed) ^ ((numpy.uint64(step) * numpy.uint64(4096) + numpy.uint64(link)) * numpy.uint64(4)
                                    + numpy.uint64(draw))
        total += float(_split_mix(key) >> numpy.uint64(11)) * (1.0 / 9007199254740992.0)
    return (total - 2.0) * numpy.sqrt(3.0)


@njit(cache=True)
def _narrowest(areas):
    """The reversible jet's narrowest node (-1: none) and its area."""
    narrowest = -1
    least = 0.6
    for node in range(10, areas.shape[0]):
        if areas[node] < least:
            least = areas[node]
            narrowest = node
    return narrowest, max(least, 1e-4)


@njit(cache=True)
def _jet_source(link, previous_flow, narrowest, area, state, gain, seed, step, closing=0.0):
    """One link's reversible-jet source for `step` (its shaping filter,
    state[4], advanced): VoiceBody::reversible_jet for that link (`closing`
    its reversible_jet_closing, sustain REVJETCLOSE)."""
    velocity = abs(previous_flow) / area
    slit = max(area / 1.2, 0.01)
    corner = min(max(0.2 * velocity / slit, 1500.0), 12000.0)
    corner_gain = 1.0 - numpy.exp(-2.0 * numpy.pi * corner / 96000.0)
    floor_gain = 1.0 - numpy.exp(-2.0 * numpy.pi * corner / (3.0 * 96000.0))
    pole = 1.0 - corner_gain
    gain_squared = corner_gain * corner_gain
    pole_term = 1.0 - pole * pole
    passed = gain_squared * gain_squared * (1.0 + pole * pole) / (pole_term * pole_term * pole_term)
    white = _hashed_gauss(seed, link, step)
    state[0] += floor_gain * (white - state[0])
    once = white - state[0]
    state[1] += floor_gain * (once - state[1])
    state[2] += corner_gain * ((once - state[1]) - state[2])
    state[3] += corner_gain * (state[2] - state[3])
    if narrowest < 0 or link < narrowest:
        return 0.0
    diameter = numpy.sqrt(4.0 * area / numpy.pi)
    reynolds = velocity * diameter * 1.14e-3 / 1.86e-4
    onset = 1800.0 / max(reynolds, 1.0)
    turbulent = max(1.0 - onset * onset, 0.0)
    fade = min(max((0.6 - area) / 0.3, 0.0), 1.0)
    if closing > 0.0:
        closed = 0.25 * closing
        fade *= min(max((area - closed) / (closing - closed), 0.0), 1.0)
    dynamic = 0.5 * 1.14e-3 * velocity * velocity * turbulent * fade * fade * (velocity / 2000.0)
    distance = (link - narrowest) * (35000.0 / 96000.0)
    strength = gain * dynamic * numpy.sqrt(area) * numpy.exp(-distance / 1.0)
    return strength * state[3] / numpy.sqrt(passed)


@njit(cache=True)
def _forward_moving(geometry, flow_in, hop, kinetic, end_extra, jet_gain, jet_seed, jet_gate=None):
    steps = flow_in.shape[0]
    n = geometry.shape[1]
    rho_c = 1.14e-3 * 35000.0
    f = numpy.zeros(n - 1)
    b = numpy.zeros(n - 1)
    f_smooth = numpy.zeros(n - 1)
    b_smooth = numpy.zeros(n - 1)
    q = numpy.zeros(n)
    v = numpy.zeros(n)
    bq = numpy.zeros(n)
    bv = numpy.zeros(n)
    p = numpy.zeros(n)
    y = numpy.zeros(n - 1)
    inertance_flow = 0.0
    previous = 0.0
    out = numpy.zeros(steps)
    link_flow = numpy.zeros(n - 1)
    sources = numpy.zeros(n - 1)
    jet_states = numpy.zeros((n - 1, 4))
    for step in range(steps):
        if jet_gain > 0.0:
            narrowest, least = _narrowest(geometry[max(step - 1, 0), :, 0])
            gain = jet_gain
            if jet_gate is not None:
                gain = jet_gain * jet_gate[min(step, jet_gate.shape[0] - 1)]
            for i in range(n - 1):
                sources[i] = _jet_source(i, link_flow[i], narrowest, least, jet_states[i], gain, jet_seed, step)
        before = geometry[max(step - 1, 0), :, 3:]
        walls = geometry[step, :, 3:]
        _lining_springs(before, q, v, bq, bv)
        for i in range(n - 1):
            y[i] = 0.5 * (geometry[step, i, 0] + geometry[step, i + 1, 0]) / rho_c
        p[0] = _junction(0, walls, flow_in[step] + 2.0 * y[0] * b[0] - geometry[step, 0, 1], y[0], v, bv)
        for i in range(1, n - 1):
            p[i] = _junction(i, walls, 2.0 * y[i - 1] * f[i - 1] + 2.0 * y[i] * b[i] - geometry[step, i, 1],
                             y[i - 1] + y[i], v, bv)
        resistance, gain = _radiation_gain(geometry[step, n - 1, 0], end_extra)
        known = 2.0 * y[n - 2] * f[n - 2] - inertance_flow - gain * previous - geometry[step, n - 1, 1]
        p[n - 1] = _junction(n - 1, walls, known, y[n - 2] + 1.0 / resistance + gain, v, bv)
        inertance_flow += gain * (p[n - 1] + previous)
        previous = p[n - 1]
        out[step] = p[n - 1] / resistance + inertance_flow
        for i in range(n - 1):
            f_smooth[i] = (1.0 - hop) * (p[i] - b[i]) + hop * f_smooth[i]
            b_smooth[i] = (1.0 - hop) * (p[i + 1] - f[i]) + hop * b_smooth[i]
        for i in range(n - 1):
            link_flow[i], f[i], b[i] = _series(f_smooth[i], b_smooth[i], 1.0 / y[i], geometry[step, i, 2], kinetic,
                                               geometry[step, i, 0], geometry[step, i + 1, 0], sources[i])
    return out


@njit(cache=True)
def _peel_moving_side(geometry, flow_out, hop, kinetic, end_extra, jet_gain, jet_seed, port_node, side_in, side_y,
                      stop_node, jet_gate=None):
    """_peel_moving with a side port at port_node (the nose's: side_in the wave
    it sends into the node, side_y its admittance, per step); with stop_node
    >= 0, returns once that node's pressure is known (pressure then holds it)."""
    count = flow_out.shape[0]
    n = geometry.shape[1]
    rho_c = 1.14e-3 * 35000.0
    q = numpy.zeros(n)
    v = numpy.zeros(n)
    bq = numpy.zeros(n)
    bv = numpy.zeros(n)
    pressure = numpy.zeros(count)
    arriving = numpy.zeros(count)
    flow_l = 0.0
    previous = 0.0
    for t in range(count):
        resistance, gain = _radiation_gain(geometry[t, n - 1, 0], end_extra)
        y = 0.5 * (geometry[t, n - 2, 0] + geometry[t, n - 1, 0]) / rho_c
        pressure[t] = (flow_out[t] - flow_l - gain * previous) / (1.0 / resistance + gain)
        _lining_one(n - 1, geometry[max(t - 1, 0), :, 3:], q, v, bq, bv)
        middle = _wall_middle(n - 1, geometry[t, :, 3:], pressure[t], v, bv)
        known = pressure[t] * (y + 1.0 / resistance + gain) + geometry[t, n - 1, 3] * middle
        arriving[t] = (known + flow_l + gain * previous + geometry[t, n - 1, 1]) / (2.0 * y)
        flow_l += gain * (pressure[t] + previous)
        previous = pressure[t]
    rightward_in = numpy.zeros(count)
    leftward_out = numpy.zeros(count)
    # The reversible jet's narrowest node and area before each step.
    narrowest = numpy.full(count, -1)
    least = numpy.ones(count)
    if jet_gain > 0.0:
        for t in range(count):
            narrowest[t], least[t] = _narrowest(geometry[max(t - 1, 0), :, 0])
    jet_state = numpy.zeros(4)
    for j in range(n - 2, -1, -1):
        # Link j: the flow it carried as it scattered at step t - 1, on that
        # step's geometry; from it the wave that had entered from node j
        # (less half its jet source, known from its own flow the step before).
        smooth = 0.0
        previous_flow = 0.0
        jet_state[:] = 0.0
        for t in range(count):
            raw = (pressure[t - 1] - arriving[t - 1]) if t > 0 else 0.0
            smooth = (1.0 - hop) * raw + hop * smooth
            made = max(t - 1, 0)
            first = geometry[made, j, 0]
            second = geometry[made, j + 1, 0]
            z = rho_c / (0.5 * (first + second))
            flow = (arriving[t] - smooth) / z
            entered = smooth + (z + 0.5 * geometry[made, j, 2]) * flow
            if kinetic:
                entered += _kinetic_term(flow, first, second)
            if jet_gain > 0.0 and t > 0:
                gain = jet_gain
                if jet_gate is not None:
                    gain = jet_gain * jet_gate[min(made, jet_gate.shape[0] - 1)]
                entered -= 0.5 * _jet_source(j, previous_flow, narrowest[made], least[made], jet_state, gain,
                                             jet_seed, made)
            previous_flow = flow
            rightward_in[t] = entered
            leftward_out[t] = entered - z * flow
        for t in range(count):
            following = rightward_in[t + 1] if t + 1 < count else 0.0
            pressure[t] = (following - hop * rightward_in[t]) / (1.0 - hop) + leftward_out[t]
        if j == stop_node:
            return arriving, pressure
        for t in range(count):
            _lining_one(j, geometry[max(t - 1, 0), :, 3:], q, v, bq, bv)
            middle = _wall_middle(j, geometry[t, :, 3:], pressure[t], v, bv)
            right = 0.5 * (geometry[t, j, 0] + geometry[t, j + 1, 0]) / rho_c
            if j > 0:
                left = 0.5 * (geometry[t, j - 1, 0] + geometry[t, j, 0]) / rho_c
                side = side_y[t] if j == port_node else 0.0
                known = pressure[t] * (left + right + side) + geometry[t, j, 3] * middle
                if side > 0.0:
                    known -= 2.0 * side * side_in[t]
                arriving[t] = (known + geometry[t, j, 1] - 2.0 * right * leftward_out[t]) / (2.0 * left)
            else:
                arriving[t] = pressure[t] * right + geometry[t, 0, 3] * middle + geometry[t, 0, 1] - 2.0 * right * leftward_out[t]
    return arriving, pressure


@njit(cache=True)
def _peel_moving(geometry, flow_out, hop, kinetic, end_extra, jet_gain, jet_seed, jet_gate=None):
    none = numpy.zeros(1)
    return _peel_moving_side(geometry, flow_out, hop, kinetic, end_extra, jet_gain, jet_seed, -1, none, none, -1,
                             jet_gate)


@njit(cache=True)
def _peel_nose(areas, walls, targets, glide_start, hop, kinetic, flow_out, port_y, hop_zero=0.0):
    """The nose peeled from its nostrils' flow: the wave the port delivered
    into its first node each step (0 while the port is shut) and that node's
    pressure."""
    count = flow_out.shape[0]
    m = areas.shape[0]
    rho_c = 1.14e-3 * 35000.0
    q = numpy.zeros(m)
    v = numpy.zeros(m)
    bq = numpy.zeros(m)
    bv = numpy.zeros(m)
    pressure = numpy.zeros(count)
    arriving = numpy.zeros(count)
    resistance, gain = _radiation_gain(areas[m - 1], 0.0)
    y_end = 0.5 * (areas[m - 2] + areas[m - 1]) / rho_c
    flow_l = 0.0
    previous = 0.0
    for t in range(count):
        pressure[t] = (flow_out[t] - flow_l - gain * previous) / (1.0 / resistance + gain)
        _lining_one(m - 1, walls, q, v, bq, bv)
        middle = _wall_middle(m - 1, walls, pressure[t], v, bv)
        known = pressure[t] * (y_end + 1.0 / resistance + gain) + walls[m - 1, 0] * middle
        arriving[t] = (known + flow_l + gain * previous) / (2.0 * y_end)
        flow_l += gain * (pressure[t] + previous)
        previous = pressure[t]
    rightward_in = numpy.zeros(count)
    leftward_out = numpy.zeros(count)
    for j in range(m - 2, -1, -1):
        smooth = 0.0
        raw_before = 0.0
        shelf = (1.0 - hop) / (1.0 - hop_zero)
        z = rho_c / (0.5 * (areas[j] + areas[j + 1]))
        for t in range(count):
            raw = (pressure[t - 1] - arriving[t - 1]) if t > 0 else 0.0
            smooth = shelf * (raw - hop_zero * raw_before) + hop * smooth
            raw_before = raw
            made = max(t - 1, 0)
            share = min((glide_start + made + 1) / 16.0, 1.0)
            flow = (arriving[t] - smooth) / z
            entered = smooth + (z + 0.5 * share * targets[j]) * flow
            if kinetic:
                entered += _kinetic_term(flow, areas[j], areas[j + 1])
            rightward_in[t] = entered
            leftward_out[t] = entered - z * flow
        unsmoothed = 0.0
        for t in range(count):
            following = rightward_in[t + 1] if t + 1 < count else 0.0
            unsmoothed = (following - hop * rightward_in[t]) / shelf + hop_zero * unsmoothed
            pressure[t] = unsmoothed + leftward_out[t]
        for t in range(count):
            _lining_one(j, walls, q, v, bq, bv)
            middle = _wall_middle(j, walls, pressure[t], v, bv)
            right = 0.5 * (areas[j] + areas[j + 1]) / rho_c
            if j > 0:
                left = 0.5 * (areas[j - 1] + areas[j]) / rho_c
                known = pressure[t] * (left + right) + walls[j, 0] * middle
                arriving[t] = (known - 2.0 * right * leftward_out[t]) / (2.0 * left)
            else:
                known = pressure[t] * (right + port_y[t]) + walls[0, 0] * middle
                arriving[t] = (known - 2.0 * right * leftward_out[t]) / (2.0 * port_y[t]) if port_y[t] > 0.0 else 0.0
    return arriving, pressure


CONTROL_STEPS = 12      # the body's control tick (units.hpp kControlSteps)
FRAME_CORNER = 250.0    # Hz, the pushed air's two poles (Tube frame_gain_)
FRAME_OMEGA = 2.0 * numpy.pi * 160.0   # the frame's critically damped follower
FAST_FRAME = 0.5        # cm/s: a node faster than this has exact viscous links


def material_table(count):
    """Each node's material row (MATERIALS order) as lining_constants assigns them."""
    total = count * NODE_LENGTH
    out = numpy.zeros((count, 9))
    for node in range(count):
        x = (node + 0.5) * NODE_LENGTH
        out[node] = MATERIALS["firm" if x < 2.0 else ("lip" if x > total - 1.0 else "soft")]
    return out


@njit(cache=True)
def _compile_lining(radius, material, row):
    """Tube::compile_node into row[3:12] (support factor 1)."""
    surface = 2.0 * numpy.pi * radius * (35000.0 / 96000.0)
    cover_mass = material[0] * material[1] * surface
    body_mass = material[5] * material[6] * surface
    hoop_k = material[2] * surface / (radius * radius)
    coupling_k = material[4] * surface
    support_k = material[7] * surface
    reduced = cover_mass * body_mass / (cover_mass + body_mass)
    row[3] = surface
    row[4] = cover_mass
    row[5] = body_mass
    row[6] = hoop_k
    row[7] = coupling_k
    row[8] = support_k
    row[9] = 2.0 * material[3] * numpy.sqrt(reduced * (hoop_k + coupling_k))
    row[10] = 2.0 * material[8] * numpy.sqrt(body_mass * support_k)
    row[11] = radius


@njit(cache=True)
def _viscous_link(area, scale):
    safe = max(area, 1e-4)
    return scale * max(0.03 * safe ** -2.6, 0.0065 / safe * numpy.sqrt(2250.0)) * (35000.0 / 96000.0)


@njit(cache=True)
def _moving_geometry(rest, materials, viscous, engine_phase=False):
    """Tube::refresh_geometry for a schedule of rest areas (per step, as set
    before the step; the body's tick every CONTROL_STEPS ramps toward them):
    the GEODUMP rows. Static admittance, no conserved span, no quiet nodes,
    no area limit, the hop loss by constant."""
    steps = rest.shape[0]
    n = rest.shape[1]
    h = 1.0 / 96000.0
    omega = 2.0 * numpy.pi * 160.0
    gain = 1.0 - numpy.exp(-2.0 * numpy.pi * 250.0 * h)
    out = numpy.zeros((steps, n, 12))
    target = numpy.sqrt(numpy.maximum(rest[0], 1e-4) / numpy.pi)
    radius = target.copy()
    goal = target.copy()
    ramp = numpy.zeros(n)
    rate = numpy.zeros(n)
    low = numpy.zeros(n)
    lower = numpy.zeros(n)
    fast = numpy.zeros(n, dtype=numpy.bool_)
    resistance = numpy.zeros(n - 1)
    resistance_step = numpy.zeros(n - 1)
    ramp_left = 0
    countdown = 0
    if engine_phase:
        # The links' viscous glide as Tube's constructor leaves it: one glide
        # step already taken from nothing (a sixteenth of the law's), the
        # same step again to come, countdown 15. Started at 0 with countdown
        # 0 the resistances ran a step behind GEODUMP's (up to 0.019) and the
        # truth's equation error floored at -160.2 dB; in phase 8.9e-16 and
        # -244.3 dB (/aa iy ae/).
        for i in range(n - 1):
            link = 0.5 * (max(rest[0, i], 1e-4) + max(rest[0, i + 1], 1e-4))
            resistance_step[i] = _viscous_link(link, viscous) / 16.0
            resistance[i] = resistance_step[i]
        countdown = 15
    for step in range(steps):
        if step % 12 == 0:
            for i in range(n):
                goal[i] = numpy.sqrt(max(rest[step, i], 1e-4) / numpy.pi)
                ramp[i] = (goal[i] - target[i]) / 12.0
            ramp_left = 12
        if ramp_left > 0:
            for i in range(n):
                target[i] += ramp[i]
            ramp_left -= 1
            if ramp_left == 0:
                for i in range(n):
                    target[i] = goal[i]
        for i in range(n):
            offset = target[i] - radius[i]
            if abs(offset) >= 1e-9 or abs(rate[i]) >= 1e-6:
                rate[i] += h * (omega * omega * offset - 2.0 * omega * rate[i])
                radius[i] = max(radius[i] + h * rate[i], 1e-4)
                fast[i] = abs(rate[i]) > 0.5
            else:
                rate[i] = 0.0
                fast[i] = False
            _compile_lining(radius[i], materials[i], out[step, i])
            pushed = out[step, i, 3] * rate[i]
            low[i] += gain * (pushed - low[i])
            lower[i] += gain * (low[i] - lower[i])
            out[step, i, 1] = lower[i]
            out[step, i, 0] = max(numpy.pi * radius[i] * radius[i], 1e-4)
        for i in range(n - 1):
            link = 0.5 * (out[step, i, 0] + out[step, i + 1, 0])
            if fast[i] or fast[i + 1]:
                resistance[i] = _viscous_link(link, viscous)
                resistance_step[i] = 0.0
            else:
                if countdown == 0:
                    resistance_step[i] = (_viscous_link(link, viscous) - resistance[i]) / 16.0
                resistance[i] += resistance_step[i]
            out[step, i, 2] = resistance[i]
        countdown = 15 if countdown == 0 else countdown - 1
    return out


def moving_geometry(rest, viscous=1.0):
    """The moving tube's geometry (GEODUMP rows) for rest areas per step."""
    rest = numpy.ascontiguousarray(rest, float)
    return _moving_geometry(rest, material_table(rest.shape[1]), float(viscous), ENGINE_GLIDE_PHASE != 0)


def velum_port(times, areas, steps):
    """The velum port's area per step for sustain's VELUMSEQ (set at every
    96th step, linear between) as VoiceBody::exchange_branches glides it:
    two stages of 1 ms, under 1e-6 cm^2 shut."""
    times = numpy.asarray(times, float)
    areas = numpy.asarray(areas, float)
    gain = 1.0 - numpy.exp(-H / 0.001)
    guide = 0.0
    now = 0.0
    target = 0.0
    out = numpy.zeros(steps)
    for step in range(steps):
        if step % 96 == 0:
            target = max(float(numpy.interp(step / STEP_RATE, times, areas)), 0.0)
        guide += gain * (target - guide)
        now += gain * (guide - now)
        if now < 1e-6:
            now = 0.0
        out[step] = now
    return out


def keyframe_schedule(times, tubes, steps):
    """Rest areas per step for sustain's AREASEQ: set at every 96th step,
    linear between the keyframes (times in s, ascending; tubes one row
    each), held before the first and after the last."""
    times = numpy.asarray(times, float)
    tubes = numpy.asarray(tubes, float)
    ticks = (numpy.arange(steps) // 96 * 96) / STEP_RATE
    out = numpy.empty((steps, tubes.shape[1]))
    for node in range(tubes.shape[1]):
        out[:, node] = numpy.interp(ticks, times, tubes[:, node])
    return out


def glide_schedule(first, second, steps, start, duration):
    """Rest areas per step for sustain's GLIDE: set at every 96th step,
    linear from `first` to `second` over [start, start + duration] (s)."""
    first = numpy.asarray(first, float)
    second = numpy.asarray(second, float)
    times = (numpy.arange(steps) // 96 * 96) / STEP_RATE
    share = numpy.clip((times - start) / duration, 0.0, 1.0)
    return first[None, :] + share[:, None] * (second - first)[None, :]


# The nose (sustain REVNOSE, NOSEDUMP, PORTDUMP): a fixed tube opening from
# a vocal node through the velum's port, radiating at the nostrils. Each
# step the port, a one-hop series element, joins the vocal node's side wave
# and the nose's first wave of the step before (VoiceBody::exchange_branches);
# then the vocal tube steps with the port as a side branch at that node; then
# the nose steps (skipped once the port has been shut 50 ms); the sound is
# both openings' flows.
QUIET_STEPS = 4800


def read_nose(path):
    """sustain's NOSEDUMP: (port node, areas, link resistance targets per
    the viscous law, linings [node, 9], hop low-pass, viscous scale, the
    links' resistance after the first step)."""
    with open(path) as handle:
        head = handle.readline().split()
        rows = numpy.array([[float(v) for v in line.split()] for line in handle if line.strip()])
    port_node = int(head[0])
    hop = float(head[2])
    scale = float(head[3])
    areas = rows[:, 0]
    links = 0.5 * (areas[:-1] + areas[1:])
    targets = scale * viscous_per_cm(links) * NODE_LENGTH
    return port_node, areas, targets, rows[:, 3:12], hop, scale, rows[:-1, 2]


@njit(cache=True)
def _nose_step(areas, walls, targets, glide_steps, hop, kinetic, incoming, port_area, f, b, f_smooth, b_smooth,
               q, v, bq, bv, p, state, rad_mass, rad_resistance, hop_zero, f_raw, b_raw):
    """One step of the nose (Tube::step for a fixed tube): state[0] the
    radiation's inertance flow, state[1] its last pressure; returns (flow
    radiated, wave leaving toward the port)."""
    m = areas.shape[0]
    rho_c = 1.14e-3 * 35000.0
    _lining_springs(walls, q, v, bq, bv)
    share = min(glide_steps / 16.0, 1.0)
    y_port = port_area / rho_c
    y0 = 0.5 * (areas[0] + areas[1]) / rho_c
    admittance = y0 + y_port
    if admittance <= 0.0:
        admittance = 1e-4 / rho_c
    p[0] = _junction(0, walls, 2.0 * y_port * incoming + 2.0 * y0 * b[0], admittance, v, bv)
    for i in range(1, m - 1):
        yl = 0.5 * (areas[i - 1] + areas[i]) / rho_c
        yr = 0.5 * (areas[i] + areas[i + 1]) / rho_c
        p[i] = _junction(i, walls, 2.0 * yl * f[i - 1] + 2.0 * yr * b[i], yl + yr, v, bv)
    resistance, gain = _radiation_gain_scaled(areas[m - 1], 0.0, rad_mass, rad_resistance)
    yl = 0.5 * (areas[m - 2] + areas[m - 1]) / rho_c
    known = 2.0 * yl * f[m - 2] - state[0] - gain * state[1]
    p[m - 1] = _junction(m - 1, walls, known, yl + 1.0 / resistance + gain, v, bv)
    state[0] += gain * (p[m - 1] + state[1])
    state[1] = p[m - 1]
    radiated = p[m - 1] / resistance + state[0]
    leaving = p[0] - incoming
    shelf = (1.0 - hop) / (1.0 - hop_zero)
    for i in range(m - 1):
        wave = p[i] - b[i]
        f_smooth[i] = shelf * (wave - hop_zero * f_raw[i]) + hop * f_smooth[i]
        f_raw[i] = wave
        wave = p[i + 1] - f[i]
        b_smooth[i] = shelf * (wave - hop_zero * b_raw[i]) + hop * b_smooth[i]
        b_raw[i] = wave
    for i in range(m - 1):
        link = 0.5 * (areas[i] + areas[i + 1])
        flow, f[i], b[i] = _series(f_smooth[i], b_smooth[i], rho_c / link, share * targets[i], kinetic,
                                   areas[i], areas[i + 1], 0.0)
    return radiated, leaving


@njit(cache=True)
def _forward_branched(geometry, flow_in, hop, kinetic, end_extra, port, port_node, nose_areas, nose_walls,
                      nose_targets, nose_hop, glide_start, seal, rad_mass, rad_resistance, port_scale, hop_zero):
    """_forward_moving with the nose: the flow radiated by lips and nostrils."""
    steps = flow_in.shape[0]
    n = geometry.shape[1]
    m = nose_areas.shape[0]
    rho_c = 1.14e-3 * 35000.0
    f = numpy.zeros(n - 1)
    b = numpy.zeros(n - 1)
    f_smooth = numpy.zeros(n - 1)
    b_smooth = numpy.zeros(n - 1)
    q = numpy.zeros(n)
    v = numpy.zeros(n)
    bq = numpy.zeros(n)
    bv = numpy.zeros(n)
    p = numpy.zeros(n)
    y = numpy.zeros(n - 1)
    nf = numpy.zeros(m - 1)
    nb = numpy.zeros(m - 1)
    nf_smooth = numpy.zeros(m - 1)
    nb_smooth = numpy.zeros(m - 1)
    nq = numpy.zeros(m)
    nv = numpy.zeros(m)
    nbq = numpy.zeros(m)
    nbv = numpy.zeros(m)
    np_ = numpy.zeros(m)
    nose_state = numpy.zeros(2)
    nose_f_raw = numpy.zeros(m - 1)
    nose_b_raw = numpy.zeros(m - 1)
    inertance_flow = 0.0
    previous = 0.0
    side_out = 0.0
    nose_left = 0.0
    nose_radiated = 0.0
    quiet = 0
    glide_steps = glide_start
    out = numpy.zeros(steps)
    nostrils = numpy.zeros(steps)
    delivered_wave = numpy.zeros(steps)
    returned_wave = numpy.zeros(steps)
    port_pressure = numpy.zeros(steps)
    for step in range(steps):
        # The port.
        area = port[step]
        into_owner = 0.0
        into_child = 0.0
        open_area = 0.0
        if area >= 1e-6:
            open_area = area
            z_port = rho_c / area
            r_port = port_scale * max(0.03 * max(area, 1e-4) ** -2.6, 0.0065 / max(area, 1e-4) * numpy.sqrt(2250.0)) * (35000.0 / 96000.0)
            flow, into_child, into_owner = _series(side_out, nose_left, z_port, r_port, kinetic,
                                                   geometry[max(step - 1, 0), port_node, 0], nose_areas[0], 0.0)
        # The vocal tube, the port a side branch at port_node.
        before = geometry[max(step - 1, 0), :, 3:]
        walls = geometry[step, :, 3:]
        _lining_springs(before, q, v, bq, bv)
        for i in range(n - 1):
            y[i] = 0.5 * (geometry[step, i, 0] + geometry[step, i + 1, 0]) / rho_c
        y_side = open_area / rho_c
        for i in range(n):
            if i == 0:
                known = flow_in[step] + 2.0 * y[0] * b[0] - geometry[step, 0, 1]
                admittance = y[0]
            elif i < n - 1:
                known = 2.0 * y[i - 1] * f[i - 1] + 2.0 * y[i] * b[i] - geometry[step, i, 1]
                admittance = y[i - 1] + y[i]
            else:
                resistance, gain = _radiation_gain(geometry[step, n - 1, 0], end_extra)
                known = 2.0 * y[n - 2] * f[n - 2] - inertance_flow - gain * previous - geometry[step, n - 1, 1]
                admittance = y[n - 2] + 1.0 / resistance + gain
            if i == port_node and open_area > 0.0:
                known += 2.0 * y_side * into_owner
                admittance += y_side
            p[i] = _junction(i, walls, known, admittance, v, bv)
        side_out = p[port_node] - into_owner if open_area > 0.0 else 0.0
        resistance, gain = _radiation_gain(geometry[step, n - 1, 0], end_extra)
        inertance_flow += gain * (p[n - 1] + previous)
        previous = p[n - 1]
        lips = p[n - 1] / resistance + inertance_flow
        for i in range(n - 1):
            f_smooth[i] = (1.0 - hop) * (p[i] - b[i]) + hop * f_smooth[i]
            b_smooth[i] = (1.0 - hop) * (p[i + 1] - f[i]) + hop * b_smooth[i]
        for i in range(n - 1):
            if seal > 0.0 and (geometry[step, i, 0] <= seal or geometry[step, i + 1, 0] <= seal):
                f[i] = b_smooth[i]
                b[i] = f_smooth[i]
                continue
            flow, f[i], b[i] = _series(f_smooth[i], b_smooth[i], 1.0 / y[i], geometry[step, i, 2], kinetic,
                                       geometry[step, i, 0], geometry[step, i + 1, 0], 0.0)
        # The nose (its first step with the port as it now stands).
        quiet = 0 if area > 0.0 else quiet + 1
        if quiet <= QUIET_STEPS:
            glide_steps += 1
            nose_radiated, nose_left = _nose_step(nose_areas, nose_walls, nose_targets, glide_steps, nose_hop, kinetic,
                                                  into_child, open_area, nf, nb, nf_smooth, nb_smooth, nq, nv, nbq,
                                                  nbv, np_, nose_state, rad_mass, rad_resistance, hop_zero, nose_f_raw,
                                                  nose_b_raw)
        out[step] = lips + nose_radiated
        nostrils[step] = nose_radiated
        delivered_wave[step] = into_child
        returned_wave[step] = into_owner
        port_pressure[step] = p[port_node]
    return out, nostrils, delivered_wave, returned_wave, port_pressure


def forward_branched(geometry, flow_in, port, nose, hop=0.0, kinetic=True, end_extra=0.0, seal=0.0, radiation=(1.0, 1.0),
                     port_scale=1.0, nose_hop_zero=0.0, nose_hop_value=None):
    """The flow radiated by lips and nostrils (nose: read_nose's tuple; port:
    the port's area per step)."""
    port_node, areas, targets, walls, nose_hop, scale, first = nose
    # The nose's links glide from nothing to their viscous resistance over 16
    # steps; how far after the first step, from the dump.
    glide_start = int(round(16.0 * float(numpy.mean(first / targets)))) - 1
    return _forward_branched(numpy.ascontiguousarray(geometry, float), numpy.asarray(flow_in, float), float(hop),
                             bool(kinetic), float(end_extra), numpy.asarray(port, float), int(port_node),
                             numpy.asarray(areas, float), numpy.ascontiguousarray(walls, float),
                             numpy.asarray(targets, float), float(nose[4] if nose_hop_value is None else nose_hop_value), glide_start,
                             float(seal), float(radiation[0]), float(radiation[1]), float(port_scale), float(nose_hop_zero))[0]


def forward_branched_parts(geometry, flow_in, port, nose, hop=0.0, kinetic=True, end_extra=0.0, seal=0.0,
                           radiation=(1.0, 1.0), port_scale=1.0, nose_hop_zero=0.0, nose_hop_value=None):
    """forward_branched's sound, the nostrils' part of it, the port's waves
    into the nose and back into the vocal node, and the port node's
    pressure, per step."""
    port_node, areas, targets, walls, nose_hop, scale, first = nose
    glide_start = int(round(16.0 * float(numpy.mean(first / targets)))) - 1
    return _forward_branched(numpy.ascontiguousarray(geometry, float), numpy.asarray(flow_in, float), float(hop),
                             bool(kinetic), float(end_extra), numpy.asarray(port, float), int(port_node),
                             numpy.asarray(areas, float), numpy.ascontiguousarray(walls, float),
                             numpy.asarray(targets, float), float(nose[4] if nose_hop_value is None else nose_hop_value), glide_start,
                             float(seal), float(radiation[0]), float(radiation[1]), float(port_scale), float(nose_hop_zero))


def forward_moving(geometry, flow_in, hop=0.0, kinetic=True, end_extra=0.0, jet=0.0, seed=JET_SEED, gate=None):
    """The flow radiated by a moving tube (GEODUMP geometry) for the flow
    entering; `jet` the reversible jet's gain (sustain REVJET), `gate` its
    gate per step (jet_gate_steps) or None."""
    return _forward_moving(numpy.ascontiguousarray(geometry, float), numpy.asarray(flow_in, float), float(hop),
                           bool(kinetic), float(end_extra), float(jet), int(seed),
                           None if gate is None else numpy.asarray(gate, float))


def peel_moving(geometry, flow_out, hop=0.0, kinetic=True, end_extra=0.0, jet=0.0, seed=JET_SEED, gate=None):
    """The flow entering a moving tube, peeled back from the flow radiated."""
    return peel_moving_glottis(geometry, flow_out, hop, kinetic, end_extra, jet, seed, gate)[0]


def peel_moving_glottis(geometry, flow_out, hop=0.0, kinetic=True, end_extra=0.0, jet=0.0, seed=JET_SEED, gate=None):
    """peel_moving's flow entering and the first node's pressure."""
    return _peel_moving(numpy.ascontiguousarray(geometry, float), numpy.asarray(flow_out, float), float(hop),
                        bool(kinetic), float(end_extra), float(jet), int(seed),
                        None if gate is None else numpy.asarray(gate, float))


def jet_gate_steps(times, gates, steps):
    """The reversible jet's gate per step for sustain's JETSEQ (set at every
    96th step, linear between, held outside), as VoiceBody's jet reads it."""
    times = numpy.asarray(times, float)
    gates = numpy.asarray(gates, float)
    out = numpy.zeros(steps)
    for tick in range(0, steps, 96):
        out[tick:tick + 96] = float(numpy.interp(tick / STEP_RATE, times, gates))
    return out


def snr(reference, estimate):
    reference = numpy.asarray(reference)
    error = reference - numpy.asarray(estimate)
    return 10.0 * numpy.log10(numpy.sum(reference ** 2) / max(numpy.sum(error ** 2), 1e-30))


def main():
    areas = numpy.loadtxt(sys.argv[1])
    pairs = numpy.fromfile(sys.argv[2], dtype=numpy.float32).reshape(-1, 2).astype(float)
    flow_in, flow_out = pairs[:, 0], pairs[:, 1]
    replica = forward(areas, flow_in)
    print("nodes %d, length %.2f cm, steps %d" % (len(areas), len(areas) * SOUND_SPEED / STEP_RATE, len(flow_in)))
    print("forward: replica's lip flow against the engine's: %.1f dB" % snr(flow_out, replica))
    recovered = peel(areas, flow_out)
    # (The last nodes-many steps need sound not yet radiated. The flow's
    # mean passes the tube unchanged; the voice itself is its swing.)
    usable = len(flow_in) - len(areas)
    print("backward (peeled, %d steps): glottal flow against the engine's: %.1f dB (swing %.1f dB)" % (
        usable, snr(flow_in[:usable], recovered[:usable]),
        snr(flow_in[:usable] - flow_in[:usable].mean(), recovered[:usable] - recovered[:usable].mean())))


if __name__ == "__main__":
    main()
