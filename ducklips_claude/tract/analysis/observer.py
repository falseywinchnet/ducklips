"""The reversible engine where the sound cannot reveal the inside: peel where
it can, run forward with the folds where it cannot.

A stop's closure with the velum shut hides everything behind the seal: the
sound carries nothing of it. The glottal flow there is whatever our folds
make (sustain FORCEWAVE with NaN: free-running), and the tube runs forward
on it (this module's forward step, the reversible core's replica). Before
the closure the tube is peeled node by node as ever, up to the first moment
the closure would reach the peel's lookahead (node 0's: the seal node's
index steps before the seal); from the release on it is peeled again, each
node starting from the state the forward run left.

State rows per node (a forward step's whole state after a step): 0 the
forward wave arriving at the next node, 1 the backward wave arriving at this
node (both on this node's outer link, after its scatter), 2 the outer link's
forward smoothed wave, 3 its backward smoothed wave, 4 the outer link's last
unsmoothed forward wave, 5 its last unsmoothed backward wave (the shelf's),
6-9 the lining (q, v, bq, bv), 10 the node's pressure; row 0's 11, 12 the
open end's radiation (inertance flow, last pressure); with the reversible
jet, per link: 13 its flow at the step, 14-17 its source's shaping filter
after it, 18 the source it carried then."""
import os
import sys

import numpy
from numba import njit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import marching
import reverse_core

RHO_C = 1.14e-3 * 35000.0
STATE_COLUMNS = 19


@njit(cache=True)
def _lining_all(walls, q, v, bq, bv):
    for i in range(walls.shape[0]):
        marching._lining(walls, i, q, v, bq, bv)


@njit(cache=True)
def _junction_with(walls, i, known, admittance, v, bv):
    """A node's junction with its lining (Tube::step_air): the pressure, the
    cover's velocity updated."""
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
def _forward_steps(geometry, flow_in, t0, t1, state, hop, hop_zero, kinetic, end_extra, seal, lips_out, port_pressure,
                   jet_gain=0.0, jet_seed=0, jet_closing=0.0, jet_gate=None):
    """Steps t0..t1 - 1 of the tube forward (reverse_core._forward_moving's
    step, with the shelf, the seal and the reversible jet) from `state`; the
    lips' flow and the first node's pressure per step into lips_out,
    port_pressure."""
    n = geometry.shape[1]
    q = state[:, 6]
    v = state[:, 7]
    bq = state[:, 8]
    bv = state[:, 9]
    y = numpy.zeros(n - 1)
    shelf = (1.0 - hop) / (1.0 - hop_zero)
    for step in range(t0, t1):
        if jet_gain > 0.0:
            # Each link's source from its own flow the step before (its
            # filter advanced even while sealed, as VoiceBody's).
            narrowest, least = reverse_core._narrowest(geometry[max(step - 1, 0), :, 0])
            gain = jet_gain
            if jet_gate is not None:
                gain = jet_gain * jet_gate[min(step, jet_gate.shape[0] - 1)]
            for i in range(n - 1):
                state[i, 18] = reverse_core._jet_source(i, state[i, 13], narrowest, least, state[i, 14:18], gain,
                                                        jet_seed, step, jet_closing)
        before = geometry[max(step - 1, 0), :, 3:]
        walls = geometry[step, :, 3:]
        _lining_all(before, q, v, bq, bv)
        for i in range(n - 1):
            y[i] = 0.5 * (geometry[step, i, 0] + geometry[step, i + 1, 0]) / RHO_C
        for i in range(n):
            if i == 0:
                known = flow_in[step] + 2.0 * y[0] * state[0, 1] - geometry[step, 0, 1]
                admittance = y[0]
            elif i < n - 1:
                known = 2.0 * y[i - 1] * state[i - 1, 0] + 2.0 * y[i] * state[i, 1] - geometry[step, i, 1]
                admittance = y[i - 1] + y[i]
            else:
                resistance, gain = marching._radiation(geometry[step, n - 1, 0], end_extra)
                known = 2.0 * y[n - 2] * state[n - 2, 0] - state[0, 11] - gain * state[0, 12] - geometry[step, n - 1, 1]
                admittance = y[n - 2] + 1.0 / resistance + gain
            state[i, 10] = _junction_with(walls, i, known, admittance, v, bv)
        resistance, gain = marching._radiation(geometry[step, n - 1, 0], end_extra)
        state[0, 11] += gain * (state[n - 1, 10] + state[0, 12])
        state[0, 12] = state[n - 1, 10]
        lips_out[step] = state[n - 1, 10] / resistance + state[0, 11]
        port_pressure[step] = state[0, 10]
        for i in range(n - 1):
            wave = state[i, 10] - state[i, 1]
            state[i, 2] = shelf * (wave - hop_zero * state[i, 4]) + hop * state[i, 2]
            state[i, 4] = wave
            wave = state[i + 1, 10] - state[i, 0]
            state[i, 3] = shelf * (wave - hop_zero * state[i, 5]) + hop * state[i, 3]
            state[i, 5] = wave
        for i in range(n - 1):
            first = geometry[step, i, 0]
            second = geometry[step, i + 1, 0]
            if seal > 0.0 and (first <= seal or second <= seal):
                state[i, 0] = state[i, 3]
                state[i, 1] = state[i, 2]
                state[i, 13] = 0.0
                continue
            z = 1.0 / y[i]
            drive = 2.0 * (state[i, 2] - state[i, 3])
            if jet_gain > 0.0:
                drive += state[i, 18]
            if drive >= 0.0:
                narrow, wide = first, second
            else:
                narrow, wide = second, first
            k = 0.0
            if kinetic and wide > narrow:
                difference = 1.0 / max(narrow, 1e-4) - 1.0 / wide
                k = 0.5 * 1.14e-3 * difference * difference
            linear = 2.0 * z + geometry[step, i, 2]
            magnitude = 2.0 * abs(drive) / (linear + numpy.sqrt(linear * linear + 4.0 * k * abs(drive)))
            flow = magnitude if drive >= 0.0 else -magnitude
            arriving_next = state[i, 3] + z * flow
            state[i, 1] = state[i, 2] - z * flow
            state[i, 0] = arriving_next
            state[i, 13] = flow


@njit(cache=True)
def _peel_states(geometry, sound, t_start, init, record, hop, hop_zero, kinetic, end_extra, flow_in_out, pressure_out,
                 recorded, jet_gain=0.0, jet_seed=0, jet_closing=0.0, jet_gate=None):
    """The tube peeled node by node over [t_start, len(sound)), starting (if
    t_start > 0) from `init`, the forward state after step t_start - 1. The
    flow entering and the first node's pressure into flow_in_out,
    pressure_out (from t_start on); if record >= t_start, the forward state
    after step `record` into `recorded`."""
    count = sound.shape[0]
    n = geometry.shape[1]
    shelf = (1.0 - hop) / (1.0 - hop_zero)
    q = numpy.zeros(n)
    v = numpy.zeros(n)
    bq = numpy.zeros(n)
    bv = numpy.zeros(n)
    if t_start > 0:
        for i in range(n):
            q[i] = init[i, 6]
            v[i] = init[i, 7]
            bq[i] = init[i, 8]
            bv[i] = init[i, 9]
    pressure = numpy.zeros(count + 1)
    arriving = numpy.zeros(count + 1)
    flow_l = init[0, 11] if t_start > 0 else 0.0
    previous = init[0, 12] if t_start > 0 else 0.0
    # The open end.
    for t in range(t_start, count):
        resistance, gain = marching._radiation(geometry[t, n - 1, 0], end_extra)
        y = 0.5 * (geometry[t, n - 2, 0] + geometry[t, n - 1, 0]) / RHO_C
        pressure[t] = (sound[t] - flow_l - gain * previous) / (1.0 / resistance + gain)
        marching._lining(geometry[max(t - 1, 0), :, 3:], n - 1, q, v, bq, bv)
        middle = marching._middle(geometry[t, :, 3:], n - 1, pressure[t], v, bv)
        known = pressure[t] * (y + 1.0 / resistance + gain) + geometry[t, n - 1, 3] * middle
        arriving[t] = (known + flow_l + gain * previous + geometry[t, n - 1, 1]) / (2.0 * y)
        flow_l += gain * (pressure[t] + previous)
        previous = pressure[t]
        if t == record:
            recorded[n - 1, 6] = q[n - 1]
            recorded[n - 1, 7] = v[n - 1]
            recorded[n - 1, 8] = bq[n - 1]
            recorded[n - 1, 9] = bv[n - 1]
            recorded[n - 1, 10] = pressure[t]
            recorded[0, 11] = flow_l
            recorded[0, 12] = previous
    rightward_in = numpy.zeros(count + 1)
    leftward_out = numpy.zeros(count + 1)
    # The reversible jet's narrowest node and area before each step.
    narrowest = numpy.full(count, -1)
    least = numpy.ones(count)
    if jet_gain > 0.0:
        for t in range(t_start, count):
            narrowest[t], least[t] = reverse_core._narrowest(geometry[max(t - 1, 0), :, 0])
    jet_state = numpy.zeros(4)
    for j in range(n - 2, -1, -1):
        # Link j, from its state before t_start (or from rest).
        smooth = init[j, 3] if t_start > 0 else 0.0
        raw_before = init[j, 5] if t_start > 0 else 0.0
        previous_flow = init[j, 13] if t_start > 0 else 0.0
        for c in range(4):
            jet_state[c] = init[j, 14 + c] if t_start > 0 else 0.0
        source = 0.0
        for t in range(t_start, count):
            made = max(t - 1, 0)
            if t == t_start and t_start > 0:
                pass   # its smoothed backward wave of step t - 1 is the state's
            else:
                raw = (pressure[t - 1] - arriving[t - 1]) if t > 0 else 0.0
                smooth = shelf * (raw - hop_zero * raw_before) + hop * smooth
                raw_before = raw
            first = geometry[made, j, 0]
            second = geometry[made, j + 1, 0]
            z = RHO_C / (0.5 * (first + second))
            flow = (arriving[t] - smooth) / z
            entered = smooth + (z + 0.5 * geometry[made, j, 2]) * flow
            if kinetic:
                entered += marching._kinetic(flow, first, second)
            if jet_gain > 0.0 and t > 0:
                # The source this link carried at step `made`: the state's at
                # the resumed step (that scatter is already in it), else from
                # its own flow the step before (less half: 2 (a - b) + S).
                if t == t_start:
                    source = init[j, 18]
                else:
                    gain = jet_gain
                    if jet_gate is not None:
                        gain = jet_gain * jet_gate[min(made, jet_gate.shape[0] - 1)]
                    source = reverse_core._jet_source(j, previous_flow, narrowest[made], least[made], jet_state,
                                                      gain, jet_seed, made, jet_closing)
                entered -= 0.5 * source
            previous_flow = flow
            rightward_in[t] = entered
            leftward_out[t] = entered - z * flow
            if t == record + 1:
                # After step `record`: its smoothed and unsmoothed backward
                # waves, set at this iteration (from step `record`'s values).
                recorded[j, 3] = smooth
                recorded[j, 5] = raw_before
                recorded[j, 13] = flow
                for c in range(4):
                    recorded[j, 14 + c] = jet_state[c]
                recorded[j, 18] = source
        if record >= t_start and record + 1 < count:
            # The forward state's scattered waves after step `record`: the
            # link's scatter at `record` made them (iteration record + 1).
            recorded[j, 0] = arriving[record + 1]
            recorded[j, 1] = leftward_out[record + 1]
            recorded[j, 2] = rightward_in[record + 1]
        unsmoothed = init[j, 4] if t_start > 0 else 0.0
        for t in range(t_start, count):
            following = rightward_in[t + 1] if t + 1 < count else 0.0
            unsmoothed = (following - hop * rightward_in[t]) / shelf + hop_zero * unsmoothed
            pressure[t] = unsmoothed + leftward_out[t]
            if t == record:
                recorded[j, 4] = unsmoothed
        for t in range(t_start, count):
            marching._lining(geometry[max(t - 1, 0), :, 3:], j, q, v, bq, bv)
            middle = marching._middle(geometry[t, :, 3:], j, pressure[t], v, bv)
            right = 0.5 * (geometry[t, j, 0] + geometry[t, j + 1, 0]) / RHO_C
            if j > 0:
                left = 0.5 * (geometry[t, j - 1, 0] + geometry[t, j, 0]) / RHO_C
                known = pressure[t] * (left + right) + geometry[t, j, 3] * middle
                arriving[t] = (known + geometry[t, j, 1] - 2.0 * right * leftward_out[t]) / (2.0 * left)
            else:
                flow_in_out[t] = pressure[t] * right + geometry[t, 0, 3] * middle + geometry[t, 0, 1] - 2.0 * right * leftward_out[t]
                pressure_out[t] = pressure[t]
            if t == record:
                recorded[j, 6] = q[j]
                recorded[j, 7] = v[j]
                recorded[j, 8] = bq[j]
                recorded[j, 9] = bv[j]
                recorded[j, 10] = pressure[t]


def closures(geometry, seal, hide=0.0):
    """The sealed stretches: (first sealed step, last sealed step, seal node).
    With `hide` above the seal, the stretches narrower than it (the peel
    through a narrow open node amplifies a wrong hypothesis's mismatch: at
    the seal area, /aa d aa/'s closure 0.05 ms late cost -23.6 dB against
    -209.5 at the truth)."""
    if seal <= 0.0:
        return []
    least = geometry[:, :, 0].min(axis=1)
    sealed = least <= max(seal, hide)
    out = []
    t = 0
    count = len(sealed)
    while t < count:
        if sealed[t]:
            start = t
            while t < count and sealed[t]:
                t += 1
            # Only a narrowing that seals: one that stays open (a fricative's,
            # 0.08 cm^2) peels well (/aa s aa/ exact, its times smooth).
            if numpy.any(least[start:t] <= seal):
                out.append((start, t - 1, int(numpy.argmin(geometry[start, :, 0]))))
        else:
            t += 1
    return out


def implied(geometry, sound, seal, hop, end_extra, free_flow, hop_zero=0.0, hide=0.0, jet=0.0,
            seed=reverse_core.JET_SEED, jet_closing=0.0, jet_gate=None):
    """The flow entering the vocal tube and its first node's pressure, per
    step: peeled where the sound reveals them, else the folds' own flow (from
    free_flow(forced waves), NaN where the folds run free) with the tube run
    forward on it. Returns (flow, pressure, revealed mask)."""
    geometry = numpy.ascontiguousarray(geometry, float)
    count = len(sound)
    n = geometry.shape[1]
    flow = numpy.zeros(count)
    pressure = numpy.zeros(count)
    revealed = numpy.ones(count, bool)
    modelled = numpy.zeros(count)
    state = numpy.zeros((n, STATE_COLUMNS))
    recorded = numpy.zeros((n, STATE_COLUMNS))
    start = 0
    area0 = geometry[:, 0, 0]
    for first, last, node in closures(geometry, seal, hide):
        hidden_from = max(first - node, start)
        hidden_to = min(last + 2, count)
        _peel_states(geometry, sound, start, state, hidden_from - 1, hop, hop_zero, True, end_extra, flow, pressure,
                     recorded, jet, seed, jet_closing, jet_gate)
        forced = 0.5 * (pressure - flow * RHO_C / area0)
        forced[hidden_from:] = numpy.nan
        observed = flow.copy()
        observed[~revealed] = numpy.nan
        observed[hidden_from:] = numpy.nan
        folds = free_flow(forced, observed)
        state = recorded.copy()
        lips = numpy.zeros(count)
        first_pressure = numpy.zeros(count)
        _forward_steps(geometry, folds, hidden_from, hidden_to, state, hop, hop_zero, True, end_extra, seal, lips,
                       first_pressure, jet, seed, jet_closing, jet_gate)
        flow[hidden_from:hidden_to] = folds[hidden_from:hidden_to]
        pressure[hidden_from:hidden_to] = first_pressure[hidden_from:hidden_to]
        modelled[hidden_from:hidden_to] = lips[hidden_from:hidden_to]
        revealed[hidden_from:hidden_to] = False
        start = hidden_to
    _peel_states(geometry, sound, start, state, -1, hop, hop_zero, True, end_extra, flow, pressure, recorded, jet, seed,
                 jet_closing, jet_gate)
    return flow, pressure, revealed, modelled


# The nose for the forward windows: rows per nose node as the vocal state's
# (STATE_COLUMNS); the port's carried values in port_state: 0 the vocal
# node's side wave of the step before, 1 the nose's first-node wave of the
# step before, 2 steps since the port was last open (QUIET), 3 the nose's
# glide steps (its links' viscous law reached over 16 steps).
QUIET_STEPS = 4800


@njit(cache=True)
def _nose_steps_isolated(areas, walls, targets, nstate, port_state, t0, t1, nose_hop, nose_hop_zero, out):
    """The nose, its port shut, run alone over t0..t1 - 1 (its ringing after
    the velum closes); its nostrils' flow into out."""
    for step in range(t0, t1):
        port_state[2] += 1.0
        if port_state[2] <= QUIET_STEPS:
            port_state[3] += 1.0
            out[step], port_state[1] = _nose_step_state(areas, walls, targets, port_state[3], nose_hop, nose_hop_zero,
                                                        0.0, 0.0, nstate)
        else:
            out[step] = out[step - 1] if step > 0 else 0.0


@njit(cache=True)
def _nose_step_state(areas, walls, targets, glide_steps, hop, hop_zero, incoming, port_area, nstate):
    """One step of the nose (reverse_core._nose_step) on a state laid out as
    the vocal one's. Returns (nostrils' flow, first node's wave toward the
    port)."""
    m = areas.shape[0]
    q = nstate[:, 6]
    v = nstate[:, 7]
    bq = nstate[:, 8]
    bv = nstate[:, 9]
    _lining_all(walls, q, v, bq, bv)
    share = min(glide_steps / 16.0, 1.0)
    y_port = port_area / RHO_C
    for i in range(m):
        if i == 0:
            y0 = 0.5 * (areas[0] + areas[1]) / RHO_C
            admittance = y0 + y_port
            known = 2.0 * y_port * incoming + 2.0 * y0 * nstate[0, 1]
        elif i < m - 1:
            yl = 0.5 * (areas[i - 1] + areas[i]) / RHO_C
            yr = 0.5 * (areas[i] + areas[i + 1]) / RHO_C
            admittance = yl + yr
            known = 2.0 * yl * nstate[i - 1, 0] + 2.0 * yr * nstate[i, 1]
        else:
            resistance, gain = marching._radiation(areas[m - 1], 0.0)
            yl = 0.5 * (areas[m - 2] + areas[m - 1]) / RHO_C
            admittance = yl + 1.0 / resistance + gain
            known = 2.0 * yl * nstate[m - 2, 0] - nstate[0, 11] - gain * nstate[0, 12]
        if admittance <= 0.0:
            admittance = 1e-4 / RHO_C
        nstate[i, 10] = _junction_with(walls, i, known, admittance, v, bv)
    resistance, gain = marching._radiation(areas[m - 1], 0.0)
    nstate[0, 11] += gain * (nstate[m - 1, 10] + nstate[0, 12])
    nstate[0, 12] = nstate[m - 1, 10]
    radiated = nstate[m - 1, 10] / resistance + nstate[0, 11]
    leaving = nstate[0, 10] - incoming
    shelf = (1.0 - hop) / (1.0 - hop_zero)
    for i in range(m - 1):
        wave = nstate[i, 10] - nstate[i, 1]
        nstate[i, 2] = shelf * (wave - hop_zero * nstate[i, 4]) + hop * nstate[i, 2]
        nstate[i, 4] = wave
        wave = nstate[i + 1, 10] - nstate[i, 0]
        nstate[i, 3] = shelf * (wave - hop_zero * nstate[i, 5]) + hop * nstate[i, 3]
        nstate[i, 5] = wave
    for i in range(m - 1):
        link = 0.5 * (areas[i] + areas[i + 1])
        z = RHO_C / link
        drive = 2.0 * (nstate[i, 2] - nstate[i, 3])
        if drive >= 0.0:
            narrow, wide = areas[i], areas[i + 1]
        else:
            narrow, wide = areas[i + 1], areas[i]
        k = 0.0
        if wide > narrow:
            difference = 1.0 / max(narrow, 1e-4) - 1.0 / wide
            k = 0.5 * 1.14e-3 * difference * difference
        linear = 2.0 * z + share * targets[i]
        magnitude = 2.0 * abs(drive) / (linear + numpy.sqrt(linear * linear + 4.0 * k * abs(drive)))
        flow = magnitude if drive >= 0.0 else -magnitude
        arriving_next = nstate[i, 3] + z * flow
        nstate[i, 1] = nstate[i, 2] - z * flow
        nstate[i, 0] = arriving_next
    return radiated, leaving


@njit(cache=True)
def _forward_branched_steps(geometry, flow_in, t0, t1, state, nose_areas, nose_walls, nose_targets, nstate, port_state,
                            port, port_node, hop, hop_zero, nose_hop, nose_hop_zero, kinetic, end_extra, seal, lips_out,
                            nose_out, first_pressure):
    """_forward_steps with the nose through the velum's port (VoiceBody's
    order: the port's exchange, the vocal tube with it as a side branch, the
    nose); the lips' and nostrils' flows and the first node's pressure."""
    n = geometry.shape[1]
    q = state[:, 6]
    v = state[:, 7]
    bq = state[:, 8]
    bv = state[:, 9]
    y = numpy.zeros(n - 1)
    shelf = (1.0 - hop) / (1.0 - hop_zero)
    for step in range(t0, t1):
        area = port[step]
        into_owner = 0.0
        into_child = 0.0
        open_area = 0.0
        if area >= 1e-6:
            open_area = area
            z_port = RHO_C / area
            r_port = max(0.03 * max(area, 1e-4) ** -2.6, 0.0065 / max(area, 1e-4) * numpy.sqrt(2250.0)) * (35000.0 / 96000.0)
            owner = geometry[max(step - 1, 0), port_node, 0]
            drive = 2.0 * (port_state[0] - port_state[1])
            if drive >= 0.0:
                narrow, wide = owner, nose_areas[0]
            else:
                narrow, wide = nose_areas[0], owner
            k = 0.0
            if kinetic and wide > narrow:
                difference = 1.0 / max(narrow, 1e-4) - 1.0 / wide
                k = 0.5 * 1.14e-3 * difference * difference
            linear = 2.0 * z_port + r_port
            magnitude = 2.0 * abs(drive) / (linear + numpy.sqrt(linear * linear + 4.0 * k * abs(drive)))
            flow = magnitude if drive >= 0.0 else -magnitude
            into_child = port_state[1] + z_port * flow
            into_owner = port_state[0] - z_port * flow
        before = geometry[max(step - 1, 0), :, 3:]
        walls = geometry[step, :, 3:]
        _lining_all(before, q, v, bq, bv)
        for i in range(n - 1):
            y[i] = 0.5 * (geometry[step, i, 0] + geometry[step, i + 1, 0]) / RHO_C
        y_side = open_area / RHO_C
        for i in range(n):
            if i == 0:
                known = flow_in[step] + 2.0 * y[0] * state[0, 1] - geometry[step, 0, 1]
                admittance = y[0]
            elif i < n - 1:
                known = 2.0 * y[i - 1] * state[i - 1, 0] + 2.0 * y[i] * state[i, 1] - geometry[step, i, 1]
                admittance = y[i - 1] + y[i]
            else:
                resistance, gain = marching._radiation(geometry[step, n - 1, 0], end_extra)
                known = 2.0 * y[n - 2] * state[n - 2, 0] - state[0, 11] - gain * state[0, 12] - geometry[step, n - 1, 1]
                admittance = y[n - 2] + 1.0 / resistance + gain
            if i == port_node and open_area > 0.0:
                known += 2.0 * y_side * into_owner
                admittance += y_side
            state[i, 10] = _junction_with(walls, i, known, admittance, v, bv)
        port_state[0] = state[port_node, 10] - into_owner if open_area > 0.0 else 0.0
        resistance, gain = marching._radiation(geometry[step, n - 1, 0], end_extra)
        state[0, 11] += gain * (state[n - 1, 10] + state[0, 12])
        state[0, 12] = state[n - 1, 10]
        lips_out[step] = state[n - 1, 10] / resistance + state[0, 11]
        first_pressure[step] = state[0, 10]
        for i in range(n - 1):
            wave = state[i, 10] - state[i, 1]
            state[i, 2] = shelf * (wave - hop_zero * state[i, 4]) + hop * state[i, 2]
            state[i, 4] = wave
            wave = state[i + 1, 10] - state[i, 0]
            state[i, 3] = shelf * (wave - hop_zero * state[i, 5]) + hop * state[i, 3]
            state[i, 5] = wave
        for i in range(n - 1):
            first = geometry[step, i, 0]
            second = geometry[step, i + 1, 0]
            if seal > 0.0 and (first <= seal or second <= seal):
                state[i, 0] = state[i, 3]
                state[i, 1] = state[i, 2]
                continue
            z = 1.0 / y[i]
            drive = 2.0 * (state[i, 2] - state[i, 3])
            if drive >= 0.0:
                narrow, wide = first, second
            else:
                narrow, wide = second, first
            k = 0.0
            if kinetic and wide > narrow:
                difference = 1.0 / max(narrow, 1e-4) - 1.0 / wide
                k = 0.5 * 1.14e-3 * difference * difference
            linear = 2.0 * z + geometry[step, i, 2]
            magnitude = 2.0 * abs(drive) / (linear + numpy.sqrt(linear * linear + 4.0 * k * abs(drive)))
            flow = magnitude if drive >= 0.0 else -magnitude
            arriving_next = state[i, 3] + z * flow
            state[i, 1] = state[i, 2] - z * flow
            state[i, 0] = arriving_next
        port_state[2] = 0.0 if area > 0.0 else port_state[2] + 1.0
        if port_state[2] <= QUIET_STEPS:
            port_state[3] += 1.0
            nose_out[step], port_state[1] = _nose_step_state(nose_areas, nose_walls, nose_targets, port_state[3],
                                                             nose_hop, nose_hop_zero, into_child, open_area, nstate)
        else:
            nose_out[step] = nose_out[step - 1] if step > 0 else 0.0


def hidden_windows(geometry, port, port_node, seal):
    """Stretches the sound cannot reveal for the single-path peel: each seal,
    and each time the velum's port is open, from when it first reaches node
    0's lookahead (the node's index steps early) to one step past its end."""
    count = geometry.shape[0]
    spans = []
    for first, last, node in closures(geometry, seal):
        spans.append((max(first - node, 0), min(last + 2, count)))
    opened = port >= 1e-6
    t = 0
    while t < count:
        if opened[t]:
            begin = t
            while t < count and opened[t]:
                t += 1
            spans.append((max(begin - port_node, 0), min(t + 1, count)))
        else:
            t += 1
    spans.sort()
    merged = []
    for begin, end in spans:
        if merged and begin <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((begin, end))
    return merged


def implied_branched(geometry, sound, port, nose, seal, hop, end_extra, free_flow, nose_hop, nose_hop_zero=0.0,
                     hop_zero=0.0):
    """implied() with the nose: the windows (hidden_windows) run forward with
    our folds' flow through the vocal tube, the port and the nose; outside
    them the nose, its port shut, rings on alone and its flow is taken off
    the sound before the single-path peel."""
    geometry = numpy.ascontiguousarray(geometry, float)
    port_node, nose_areas, targets, nose_walls, dumped_hop, scale, first_resistance = nose
    nose_areas = numpy.ascontiguousarray(nose_areas)
    nose_walls = numpy.ascontiguousarray(nose_walls)
    glide_start = int(round(16.0 * float(numpy.mean(first_resistance / targets)))) - 1
    count = len(sound)
    n = geometry.shape[1]
    m = len(nose_areas)
    flow = numpy.zeros(count)
    pressure = numpy.zeros(count)
    revealed = numpy.ones(count, bool)
    nostrils = numpy.zeros(count)
    state = numpy.zeros((n, STATE_COLUMNS))
    nstate = numpy.zeros((m, STATE_COLUMNS))
    port_state = numpy.array([0.0, 0.0, 0.0, float(glide_start)])
    recorded = numpy.zeros((n, STATE_COLUMNS))
    modelled = numpy.zeros(count)
    area0 = geometry[:, 0, 0]
    start = 0
    for begin, end in hidden_windows(geometry, port, port_node, seal):
        _nose_steps_isolated(nose_areas, nose_walls, targets, nstate, port_state, start, begin, nose_hop, nose_hop_zero,
                             nostrils)
        _peel_states(geometry, sound - nostrils, start, state, begin - 1, hop, hop_zero, True, end_extra, flow, pressure,
                     recorded)
        forced = 0.5 * (pressure - flow * RHO_C / area0)
        forced[begin:] = numpy.nan
        folds = free_flow(forced)
        state = recorded.copy()
        lips = numpy.zeros(count)
        first_pressure = numpy.zeros(count)
        _forward_branched_steps(geometry, folds, begin, end, state, nose_areas, nose_walls, targets, nstate, port_state,
                                port, port_node, hop, hop_zero, nose_hop, nose_hop_zero, True, end_extra, seal, lips,
                                nostrils, first_pressure)
        flow[begin:end] = folds[begin:end]
        pressure[begin:end] = first_pressure[begin:end]
        modelled[begin:end] = lips[begin:end] + nostrils[begin:end]
        revealed[begin:end] = False
        start = end
    _nose_steps_isolated(nose_areas, nose_walls, targets, nstate, port_state, start, count, nose_hop, nose_hop_zero,
                         nostrils)
    _peel_states(geometry, sound - nostrils, start, state, -1, hop, hop_zero, True, end_extra, flow, pressure, recorded)
    return flow, pressure, revealed, modelled

