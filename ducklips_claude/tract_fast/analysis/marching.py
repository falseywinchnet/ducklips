"""Time-major peels and the lips/nostrils split solved step by step.

reverse_core's peels go node by node over the whole sound. Run instead one
diagonal per new sample (node j at time tau_j = k - (last - j): each node
needs its outer neighbour one step later), a peel takes in the sound sample
by sample. The nose's share of the sound at step k then fixes, alone, the
junction's condition at step k - m (m the nose's nodes): two trial diagonals
of the nose (share 0 and 1; the condition is nearly linear in it, a Newton
step finishes) solve it. One pass, O(steps x nodes). (Matrix-free GMRES on
the whole condition stalled at 0.92 of its start after 204 products: the
condition is a causal shift of the share.)"""
import os
import sys
import weakref

import numpy
from numba import njit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kernels_jet
import reverse_core

RHO_C = 1.14e-3 * 35000.0
NODE_LENGTH = 35000.0 / 96000.0


@njit(cache=True)
def _lining(walls, i, q, v, bq, bv):
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


@njit(cache=True)
def _middle(walls, i, pressure, v, bv):
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
def _lining_at(walls_at, row, offset, i, q, v, bq, bv):
    """_lining on walls_at[row, :, offset:offset + 9] read in place (a
    view made per node and step costs numba two atomic reference counts;
    offset 3 reads a GEODUMP geometry's linings without copying them)."""
    h = 1.0 / 96000.0
    quarter = 0.25 * h * h
    mc = walls_at[row, i, offset + 1]
    if mc <= 0.0:
        return
    mb = walls_at[row, i, offset + 2]
    hoop = walls_at[row, i, offset + 3]
    coupling = walls_at[row, i, offset + 4]
    support = walls_at[row, i, offset + 5]
    cc = walls_at[row, i, offset + 6]
    bc = walls_at[row, i, offset + 7]
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
def _middle_at(walls_at, row, offset, i, pressure, v, bv):
    """_middle on walls_at[row, :, offset:offset + 9] read in place."""
    h = 1.0 / 96000.0
    mc = walls_at[row, i, offset + 1]
    if mc <= 0.0:
        return 0.0
    surface = walls_at[row, i, offset]
    cc = walls_at[row, i, offset + 6]
    following = (mc * v[i] + h * surface * pressure - h * cc * (0.5 * v[i] - bv[i])) / (mc + 0.5 * h * cc)
    middle = 0.5 * (v[i] + following)
    v[i] = following
    return middle


@njit(cache=True)
def _kinetic(flow, first, second):
    if flow >= 0.0:
        narrow, wide = first, second
    else:
        narrow, wide = second, first
    if wide <= narrow:
        return 0.0
    difference = 1.0 / max(narrow, 1e-4) - 1.0 / wide
    return 0.5 * (0.5 * 1.14e-3 * difference * difference) * flow * abs(flow)


@njit(cache=True)
def _radiation(area, end_extra):
    return _radiation_scaled(area, end_extra, 1.0, 1.0)


@njit(cache=True)
def _radiation_scaled(area, end_extra, mass_scale, resistance_scale):
    resistance = resistance_scale * 128.0 * RHO_C / (9.0 * numpy.pi ** 2 * area)
    inertance = mass_scale * 8.0 * 1.14e-3 / (3.0 * numpy.pi * numpy.sqrt(numpy.pi * area)) + 1.14e-3 * end_extra / area
    return resistance, 0.5 / 96000.0 / inertance


# A time-major peel's state, one row per node (columns): 0 q, 1 v, 2 bq,
# 3 bv (lining), 4 pressure at its current time, 5 the wave arriving at its
# current time, 6 pressure one step earlier, 7 the arriving wave one step
# earlier, 8 the hop's smoothed wave on its outer link, 9 the wave that
# entered that link at the node's current time, 10 the wave it sent back
# then; and in row 0's spare columns the open end's radiation: 11 the
# inertance flow, 12 the last end pressure; 13 the outer link's last
# unsmoothed leftward wave and 14 its last unsmoothed rightward wave (the
# hop filter as a shelf); with the reversible jet, per outer link: 15 its
# flow at its last step, 16-19 its source's shaping filter, 20 the source
# it carried then. The open end is the last node.
COLUMNS = 21


@njit(cache=True)
def _row(array, t):
    """Geometry rows by time; a fixed tube's (fewer rows) held at its last."""
    return min(max(t, 0), array.shape[0] - 1)


@njit(cache=True)
def _diagonal(state, k, sound, areas_at, walls_at, frame_at, resistance_at, hop, end_extra, side_y, stop, hop_zero=0.0,
              jet=None, walls_offset=0):
    """Take in sound[k]: node j moves to time k - (last - j), from the open
    end down to node `stop` (0: all the way). Returns, at node 0, the first
    node's junction left over for its inner port (the flow entering, for the
    glottis; twice the port's admittance times the wave it delivered, for
    the nose's port, side_y its admittance), else 0. Geometry rows by time
    (_row): areas_at [t, node], walls_at [t, node, 9], frame_at [t, node],
    resistance_at [t, link] (walls_at's linings from column walls_offset:
    a GEODUMP geometry with 3). `jet`: the reversible jet (jet_table's), or
    None."""
    n = state.shape[0]
    last = n - 1
    left_over = 0.0
    q = state[:, 0]
    v = state[:, 1]
    bq = state[:, 2]
    bv = state[:, 3]
    t = k
    if t >= 0:
        g = _row(areas_at, t)
        area = areas_at[g, last]
        resistance, gain = _radiation(area, end_extra)
        y = 0.5 * (areas_at[g, last - 1] + area) / RHO_C
        state[last, 6] = state[last, 4]
        state[last, 7] = state[last, 5]
        pressure = (sound[t] - state[0, 11] - gain * state[0, 12]) / (1.0 / resistance + gain)
        _lining_at(walls_at, _row(walls_at, t - 1), walls_offset, last, q, v, bq, bv)
        middle = _middle_at(walls_at, _row(walls_at, t), walls_offset, last, pressure, v, bv)
        known = pressure * (y + 1.0 / resistance + gain) + walls_at[_row(walls_at, t), last, walls_offset] * middle
        state[last, 4] = pressure
        state[last, 5] = (known + state[0, 11] + gain * state[0, 12] + frame_at[_row(frame_at, t), last]) / (2.0 * y)
        state[0, 11] += gain * (pressure + state[0, 12])
        state[0, 12] = pressure
    for j in range(last - 1, stop - 1, -1):
        tau = k - (last - j)
        link_time = tau + 1
        if link_time < 0:
            break
        made = max(link_time - 1, 0)
        ga = _row(areas_at, made)
        first = areas_at[ga, j]
        second = areas_at[ga, j + 1]
        z = RHO_C / (0.5 * (first + second))
        raw = (state[j + 1, 6] - state[j + 1, 7]) if link_time > 0 else 0.0
        shelf = (1.0 - hop) / (1.0 - hop_zero)
        smooth = shelf * (raw - hop_zero * state[j, 13]) + hop * state[j, 8]
        state[j, 13] = raw
        state[j, 8] = smooth
        flow = (state[j + 1, 5] - smooth) / z
        entered = smooth + (z + 0.5 * resistance_at[_row(resistance_at, made), j]) * flow + _kinetic(flow, first, second)
        if jet is not None:
            # The source the link carried at step `made`, from its own flow
            # the step before (reverse_core's peel: 2 (a - b) + S).
            if link_time > 0:
                row = _row(jet[4], made)
                source = kernels_jet.jet_source(j, state[j, 15], jet[3][row], jet[4][row], state, j, 16, jet[7],
                                                jet[0] * jet[5][row], jet[1], made, jet[2], jet[6], jet[8])
                entered -= 0.5 * source
                state[j, 20] = source
            state[j, 15] = flow
        sent = entered - z * flow
        earlier_entered = state[j, 9]
        earlier_sent = state[j, 10]
        state[j, 9] = entered
        state[j, 10] = sent
        if tau < 0:
            continue
        unsmoothed = (entered - hop * earlier_entered) / shelf + hop_zero * state[j, 14]
        state[j, 14] = unsmoothed
        pressure = unsmoothed + earlier_sent
        state[j, 6] = state[j, 4]
        state[j, 7] = state[j, 5]
        state[j, 4] = pressure
        if j == stop and stop > 0:
            state[j, 11] = earlier_sent   # the wave its outer link sent into it at tau
            break
        gt = _row(areas_at, tau)
        now = _row(walls_at, tau)
        _lining_at(walls_at, _row(walls_at, tau - 1), walls_offset, j, q, v, bq, bv)
        middle = _middle_at(walls_at, now, walls_offset, j, pressure, v, bv)
        right = 0.5 * (areas_at[gt, j] + areas_at[gt, j + 1]) / RHO_C
        frame = frame_at[_row(frame_at, tau), j]
        if j > 0:
            left = 0.5 * (areas_at[gt, j - 1] + areas_at[gt, j]) / RHO_C
            known = pressure * (left + right) + walls_at[now, j, walls_offset] * middle
            state[j, 5] = (known + frame - 2.0 * right * earlier_sent) / (2.0 * left)
        else:
            known = pressure * (right + side_y[_row(side_y, tau)]) + walls_at[now, 0, walls_offset] * middle
            left_over = known + frame - 2.0 * right * earlier_sent
            state[0, 5] = left_over
    return left_over


# The nose run forward, skewed: with the port's wave known to step k, node
# j moves to time k + j (each node needs its inner neighbour one step
# earlier), so the nostrils' flow is known m - 1 steps ahead of the port.
# State rows per node: 0 q, 1 v, 2 bq, 3 bv, 4 pressure, 5 the forward wave
# arriving from the inner link, 6 the backward wave arriving from the outer
# link, 7 the outer link's forward smoothed wave, 8 its backward smoothed
# wave; row 0 col 9, 10: the nostrils' inertance flow and last pressure.
NOSE_COLUMNS = 13   # 11, 12: the outer link's last unsmoothed forward and backward waves (the shelf's)


@njit(cache=True)
def _nose_forward(state, k, delivered, port_area, areas, walls, targets, glide_start, hop, rad_mass=1.0,
                  rad_resistance=1.0, hop_zero=0.0):
    """Advance the skewed nose by the port's wave delivered at step k.
    Returns (the wave the nose's first node sent back at step k, the
    nostrils' flow at step k + m - 1)."""
    m = areas.shape[0]
    leaving = 0.0
    radiated = 0.0
    q = state[:, 0]
    v = state[:, 1]
    bq = state[:, 2]
    bv = state[:, 3]
    for j in range(m):
        tau = k + j
        if j >= 1:
            # Link j - 1 scatters at tau - 1 (its inner smoothed wave just
            # set by node j - 1 at tau - 1, its outer one by node j then).
            link = 0.5 * (areas[j - 1] + areas[j])
            share = min((glide_start + (tau - 1) + 1) / 16.0, 1.0)
            z = RHO_C / link
            drive = 2.0 * (state[j - 1, 7] - state[j - 1, 8])
            if drive >= 0.0:
                narrow, wide = areas[j - 1], areas[j]
            else:
                narrow, wide = areas[j], areas[j - 1]
            kinetic = 0.0
            if wide > narrow:
                difference = 1.0 / max(narrow, 1e-4) - 1.0 / wide
                kinetic = 0.5 * 1.14e-3 * difference * difference
            linear = 2.0 * z + share * targets[j - 1]
            magnitude = 2.0 * abs(drive) / (linear + numpy.sqrt(linear * linear + 4.0 * kinetic * abs(drive)))
            flow = magnitude if drive >= 0.0 else -magnitude
            state[j, 5] = state[j - 1, 8] + z * flow     # forward wave reaching node j at tau
            state[j - 1, 6] = state[j - 1, 7] - z * flow # backward wave reaching node j - 1 at tau
        _lining(walls, j, q, v, bq, bv)
        if j == 0:
            y_port = port_area / RHO_C
            y0 = 0.5 * (areas[0] + areas[1]) / RHO_C
            admittance = y0 + y_port
            known = 2.0 * y_port * delivered + 2.0 * y0 * state[0, 6]
        elif j < m - 1:
            yl = 0.5 * (areas[j - 1] + areas[j]) / RHO_C
            yr = 0.5 * (areas[j] + areas[j + 1]) / RHO_C
            admittance = yl + yr
            known = 2.0 * yl * state[j, 5] + 2.0 * yr * state[j, 6]
        else:
            resistance, gain = _radiation_scaled(areas[m - 1], 0.0, rad_mass, rad_resistance)
            yl = 0.5 * (areas[m - 2] + areas[m - 1]) / RHO_C
            admittance = yl + 1.0 / resistance + gain
            known = 2.0 * yl * state[j, 5] - state[0, 9] - gain * state[0, 10]
        if admittance <= 0.0:
            admittance = 1e-4 / RHO_C
        h = 1.0 / 96000.0
        mc = walls[j, 1]
        if mc > 0.0:
            surface = walls[j, 0]
            cc = walls[j, 6]
            resist = surface * surface / admittance + cc
            push = surface * known / admittance + cc * bv[j]
            following = ((mc - 0.5 * h * resist) * v[j] + h * push) / (mc + 0.5 * h * resist)
            middle = 0.5 * (v[j] + following)
            v[j] = following
            pressure = (known - surface * middle) / admittance
        else:
            pressure = known / admittance
        state[j, 4] = pressure
        if j == 0:
            leaving = pressure - delivered
        if j == m - 1:
            resistance, gain = _radiation_scaled(areas[m - 1], 0.0, rad_mass, rad_resistance)
            state[0, 9] += gain * (pressure + state[0, 10])
            state[0, 10] = pressure
            radiated = pressure / resistance + state[0, 9]
        # Each hop's filter: g (x - hop_zero x') + hop y' (one-pole when
        # hop_zero is 0; a shelf, its loss capped at high frequency, else).
        shelf = (1.0 - hop) / (1.0 - hop_zero)
        if j < m - 1:
            wave = pressure - state[j, 6]
            state[j, 7] = shelf * (wave - hop_zero * state[j, 11]) + hop * state[j, 7]
            state[j, 11] = wave
        if j >= 1:
            wave = pressure - state[j, 5]
            state[j - 1, 8] = shelf * (wave - hop_zero * state[j - 1, 12]) + hop * state[j - 1, 8]
            state[j - 1, 12] = wave
    return leaving, radiated


@njit(cache=True)
def _march(sound, areas, walls, frame, resistance, hop, end_extra, port, port_node, nose_areas, nose_walls,
           nose_targets, glide_start, nose_hop, rad_mass=1.0, rad_resistance=1.0, port_scale=1.0, hop_zero=0.0):
    """The nose's share of the sound and the wave the port returned into the
    vocal node, step by step: the oral peel (time-major) gives the port
    node's pressure, the port steps forward, the skewed nose runs forward
    from it and gives the nostrils' flow ahead of when the oral peel needs
    it."""
    count = sound.shape[0]
    n = areas.shape[1]
    m = nose_areas.shape[0]
    oral = numpy.zeros((n, COLUMNS))
    nose = numpy.zeros((m, NOSE_COLUMNS))
    share = numpy.zeros(count + m)
    lips = numpy.zeros(count)
    returned = numpy.zeros(count)
    no_side = numpy.zeros(1)
    oral_lag = n - 1 - port_node
    side_out = 0.0
    nose_leaving = 0.0
    for t in range(count):
        # The oral peel brings the port node to time t - 1.
        ko = t - 1 + oral_lag
        if 0 <= ko < count:
            lips[ko] = sound[ko] - share[ko]
            _diagonal(oral, ko, lips, areas, walls, frame, resistance, hop, end_extra, no_side, port_node)
        if t >= 1:
            open_before = port[t - 1] >= 1e-6
            side_out = oral[port_node, 4] - returned[t - 1] if open_before and ko >= 0 else 0.0
        # The port, forward.
        delivered = 0.0
        open_area = 0.0
        if port[t] >= 1e-6:
            open_area = port[t]
            z = RHO_C / port[t]
            r_port = port_scale * max(0.03 * max(port[t], 1e-4) ** -2.6, 0.0065 / max(port[t], 1e-4) * numpy.sqrt(2250.0)) * NODE_LENGTH
            owner = areas[min(max(t - 1, 0), areas.shape[0] - 1), port_node]
            drive = 2.0 * (side_out - nose_leaving)
            if drive >= 0.0:
                narrow, wide = owner, nose_areas[0]
            else:
                narrow, wide = nose_areas[0], owner
            kinetic = 0.0
            if wide > narrow:
                difference = 1.0 / max(narrow, 1e-4) - 1.0 / wide
                kinetic = 0.5 * 1.14e-3 * difference * difference
            linear = 2.0 * z + r_port
            magnitude = 2.0 * abs(drive) / (linear + numpy.sqrt(linear * linear + 4.0 * kinetic * abs(drive)))
            flow = magnitude if drive >= 0.0 else -magnitude
            delivered = nose_leaving + z * flow
            returned[t] = side_out - z * flow
        # The nose, forward, skewed.
        nose_leaving, radiated = _nose_forward(nose, t, delivered, open_area, nose_areas, nose_walls, nose_targets,
                                               glide_start, nose_hop, rad_mass, rad_resistance, hop_zero)
        share[t + m - 1] = radiated
    return share[:count], returned


@njit(cache=True)
def _port_step(area, incoming, leaving_before, owner_area, nose_area):
    """The port's series element read backward: from the wave it delivered
    into the nose's first node at a step and the wave that node had sent the
    step before, the port's flow, the wave the vocal node sent the step
    before and the wave the port returned into it."""
    z = RHO_C / area
    resistance = max(0.03 * max(area, 1e-4) ** -2.6, 0.0065 / max(area, 1e-4) * numpy.sqrt(2250.0)) * NODE_LENGTH
    flow = (incoming - leaving_before) / z
    sent_before = leaving_before + 0.5 * ((2.0 * z + resistance) * flow) + _kinetic(flow, owner_area, nose_area)
    return sent_before, sent_before - z * flow


@njit(cache=True)
def _split_march(sound, areas, walls, frame, resistance, hop, end_extra, port, port_node, nose_areas, nose_walls,
                 nose_frame, nose_resistance, nose_hop):
    """The nose's share of the sound, step by step, and the wave the port
    returned into the vocal node each step."""
    count = sound.shape[0]
    n = areas.shape[1]
    m = nose_areas.shape[1]
    oral = numpy.zeros((n, COLUMNS))
    nose = numpy.zeros((m, COLUMNS))
    trial = numpy.zeros((m, COLUMNS))
    share = numpy.zeros(count)
    lips = sound.copy()
    returned = numpy.zeros(count)
    leaving = numpy.zeros(count)
    incoming = numpy.zeros(count)
    port_y = port / RHO_C
    no_side = numpy.zeros(1)
    oral_lag = n - 1 - port_node
    for k in range(count):
        tau = k - m
        # The oral peel to the port node at time tau (its sound known by then).
        ko = tau + oral_lag
        if ko >= 0:
            lips[ko] = sound[ko] - share[ko]
            _diagonal(oral, ko, lips, areas, walls, frame, resistance, hop, end_extra, no_side, port_node)
        if tau < 0:
            _diagonal(nose, k, share, nose_areas, nose_walls, nose_frame, nose_resistance, nose_hop, 0.0, port_y, 0)
            continue
        oral_pressure = oral[port_node, 4]
        opened = port[tau + 1] >= 1e-6
        values = numpy.zeros(3)
        conditions = numpy.zeros(3)
        values[1] = 1.0
        for attempt in range(3):
            if attempt == 2:
                slope = (conditions[1] - conditions[0]) / (values[1] - values[0])
                values[2] = values[0] - conditions[0] / slope if slope != 0.0 else 0.0
            trial[:, :] = nose[:, :]
            share[k] = values[attempt]
            left_over = _diagonal(trial, k, share, nose_areas, nose_walls, nose_frame, nose_resistance, nose_hop, 0.0,
                                  port_y, 0)
            if opened:
                delivered = left_over / (2.0 * port_y[tau + 1])
                sent_before, back = _port_step(port[tau + 1], delivered, leaving[tau],
                                               areas[min(tau, areas.shape[0] - 1), port_node], nose_areas[0, 0])
                conditions[attempt] = oral_pressure - returned[tau] - sent_before
            else:
                conditions[attempt] = left_over
        nose[:, :] = trial[:, :]
        delivered = left_over / (2.0 * port_y[tau + 1]) if opened else 0.0
        incoming[tau + 1] = delivered
        leaving[tau + 1] = nose[0, 4] - delivered
        if opened:
            sent_before, back = _port_step(port[tau + 1], delivered, leaving[tau],
                                           areas[min(tau, areas.shape[0] - 1), port_node], nose_areas[0, 0])
            returned[tau + 1] = back
    return share, returned


@njit(cache=True)
def _hybrid_march(sound, areas, walls, frame, resistance, hop, end_extra, port, port_node, nose_areas, nose_walls,
                  nose_targets, glide_start, nose_hop, shut_area, open_port):
    """The split by whichever path dominates, step by step: the inverted path
    must be the strong one. While the mouth in front of the port is shut
    (its narrowest node under shut_area) and the port is open (over
    open_port), the nose is peeled and each nostril sample solved from the
    junction's condition; otherwise the oral path is peeled and the nose run
    forward from the port. Both nose representations run on together, so
    either can take over at any step. Returns the nose's share, the wave the
    port returned into the vocal node, and the mode per step (1: nose
    peeled)."""
    count = sound.shape[0]
    n = areas.shape[1]
    m = nose_areas.shape[0]
    oral = numpy.zeros((n, COLUMNS))
    forward = numpy.zeros((m, NOSE_COLUMNS))
    peel = numpy.zeros((m, COLUMNS))
    trial = numpy.zeros((m, COLUMNS))
    share = numpy.zeros(count + m)
    lips = numpy.zeros(count)
    returned = numpy.zeros(count)
    mode = numpy.zeros(count)
    leaving = numpy.zeros(count + 1)
    port_y = port / RHO_C
    nose_areas_rows = numpy.zeros((1, m))
    nose_areas_rows[0, :] = nose_areas
    nose_walls_rows = numpy.zeros((1, m, 9))
    nose_walls_rows[0, :, :] = nose_walls
    nose_frame = numpy.zeros((1, m))
    rows = numpy.zeros((17, m - 1))
    for r in range(17):
        for j in range(m - 1):
            rows[r, j] = min((glide_start + r + 1) / 16.0, 1.0) * nose_targets[j]
    no_side = numpy.zeros(1)
    oral_lag = n - 1 - port_node
    nose_leaving = 0.0
    values = numpy.zeros(3)
    conditions = numpy.zeros(3)
    for t in range(count):
        k = t + m - 1
        ko = t - 1 + oral_lag
        if 0 <= ko < count:
            lips[ko] = sound[ko] - share[ko]
            _diagonal(oral, ko, lips, areas, walls, frame, resistance, hop, end_extra, no_side, port_node)
        side_out = 0.0
        if t >= 1 and port[t - 1] >= 1e-6 and ko >= 0:
            side_out = oral[port_node, 4] - returned[t - 1]
        opened = port[t] >= 1e-6
        g = min(t, areas.shape[0] - 1)
        narrowest = 1e9
        for node in range(port_node + 1, n):
            narrowest = min(narrowest, areas[g, node])
        owner = areas[min(max(t - 1, 0), areas.shape[0] - 1), port_node]
        delivered = 0.0
        if opened and narrowest < shut_area and port[t] > open_port and t >= 1:
            # The nose peeled; this step's nostril sample from the condition.
            mode[t] = 1.0
            oral_pressure = oral[port_node, 4] if ko >= 0 else 0.0
            values[0] = 0.0
            values[1] = 1.0
            for attempt in range(3):
                if attempt == 2:
                    slope = (conditions[1] - conditions[0]) / (values[1] - values[0])
                    values[2] = values[0] - conditions[0] / slope if slope != 0.0 else 0.0
                trial[:, :] = peel[:, :]
                share[k] = values[attempt]
                left_over = _diagonal(trial, k, share, nose_areas_rows, nose_walls_rows, nose_frame, rows, nose_hop,
                                      0.0, port_y, 0)
                arrived = left_over / (2.0 * port_y[t])
                sent_before, back = _port_step(port[t], arrived, leaving[t - 1], owner, nose_areas[0])
                conditions[attempt] = oral_pressure - returned[t - 1] - sent_before
            peel[:, :] = trial[:, :]
            delivered = left_over / (2.0 * port_y[t])
            sent_before, back = _port_step(port[t], delivered, leaving[t - 1], owner, nose_areas[0])
            returned[t] = back
            leaving[t] = peel[0, 4] - delivered
            nose_leaving, radiated = _nose_forward(forward, t, delivered, port[t], nose_areas, nose_walls,
                                                   nose_targets, glide_start, nose_hop)
        else:
            # The oral path peeled; the port and the nose forward.
            if opened:
                z = RHO_C / port[t]
                r_port = max(0.03 * max(port[t], 1e-4) ** -2.6, 0.0065 / max(port[t], 1e-4) * numpy.sqrt(2250.0)) * NODE_LENGTH
                drive = 2.0 * (side_out - nose_leaving)
                if drive >= 0.0:
                    narrow, wide = owner, nose_areas[0]
                else:
                    narrow, wide = nose_areas[0], owner
                kinetic = 0.0
                if wide > narrow:
                    difference = 1.0 / max(narrow, 1e-4) - 1.0 / wide
                    kinetic = 0.5 * 1.14e-3 * difference * difference
                linear = 2.0 * z + r_port
                magnitude = 2.0 * abs(drive) / (linear + numpy.sqrt(linear * linear + 4.0 * kinetic * abs(drive)))
                flow = magnitude if drive >= 0.0 else -magnitude
                delivered = nose_leaving + z * flow
                returned[t] = side_out - z * flow
            nose_leaving, radiated = _nose_forward(forward, t, delivered, port[t] if opened else 0.0, nose_areas,
                                                   nose_walls, nose_targets, glide_start, nose_hop)
            share[k] = radiated
            _diagonal(peel, k, share, nose_areas_rows, nose_walls_rows, nose_frame, rows, nose_hop, 0.0, port_y, 0)
            leaving[t] = peel[0, 4] - delivered if t >= 0 else 0.0
    return share[:count], returned, mode



# The oral side in front of the port run forward, skewed outward from the
# port node (node port_node + i at time tau + i), the port node's pressure
# given. Rows as the nose's (NOSE_COLUMNS) plus: 11 the row's previous
# pressure, 12 its previous forward wave, 13 its outer link's backward
# smoothed wave one step earlier; row 0's 9, 10 the open end's radiation;
# with the reversible jet, per outer link: 14 its flow at its last step,
# 15-18 its source's shaping filter, 19 the source it carried then. A sealed
# node's links turn every wave back.
ORAL_FORWARD_COLUMNS = 20


@njit(cache=True)
def _oral_forward(state, tau, port_pressure, port_node, areas, walls, frame, resistance, hop, end_extra, seal,
                  first_node=-1, jet=None, walls_offset=0):
    """Advance the front part by the port node's pressure at tau. Returns
    (the wave arriving at the port node at tau from its outer link, the lips'
    flow at tau + (last - port_node))."""
    n = areas.shape[1]
    last = n - 1
    v0 = port_node
    q = state[:, 0]
    v = state[:, 1]
    bq = state[:, 2]
    bv = state[:, 3]
    tg = min(max(tau, 0), areas.shape[0] - 1)
    arriving_port = state[v0, 6]
    start = v0 + 1 if first_node < 0 else first_node
    if start == v0 + 1:
        state[v0, 11] = state[v0, 4]
        state[v0, 4] = port_pressure
        state[v0, 7] = (1.0 - hop) * (port_pressure - state[v0, 6]) + hop * state[v0, 7]
    lips = 0.0
    for j in range(start, n):
        t = tau + (j - v0)
        g = min(max(t, 0), areas.shape[0] - 1)
        gp = min(max(t - 1, 0), areas.shape[0] - 1)
        # Link j - 1 scatters at t - 1.
        first = areas[gp, j - 1]
        second = areas[gp, j]
        source = 0.0
        if jet is not None and t - 1 >= 0:
            # Its source at step t - 1 from its own flow the step before (the
            # filter advanced even while sealed, as VoiceBody's).
            # (Its geometry's row clamped as the areas' are: the front runs
            # ahead of the sound, and past the table numba read garbage.)
            row = _row(jet[4], t - 1)
            source = kernels_jet.jet_source(j - 1, state[j - 1, 14], jet[3][row], jet[4][row], state, j - 1, 15,
                                            jet[7], jet[0] * jet[5][row], jet[1], t - 1, jet[2], jet[6], jet[8])
            state[j - 1, 19] = source
        if seal > 0.0 and (first <= seal or second <= seal):
            state[j, 12] = state[j, 5]
            state[j, 5] = state[j - 1, 8]
            state[j - 1, 6] = state[j - 1, 7]
            if jet is not None:
                state[j - 1, 14] = 0.0
        else:
            z = RHO_C / (0.5 * (first + second))
            drive = 2.0 * (state[j - 1, 7] - state[j - 1, 8]) + source
            if drive >= 0.0:
                narrow, wide = first, second
            else:
                narrow, wide = second, first
            kinetic = 0.0
            if wide > narrow:
                difference = 1.0 / max(narrow, 1e-4) - 1.0 / wide
                kinetic = 0.5 * 1.14e-3 * difference * difference
            linear = 2.0 * z + resistance[gp, j - 1]
            magnitude = 2.0 * abs(drive) / (linear + numpy.sqrt(linear * linear + 4.0 * kinetic * abs(drive)))
            flow = magnitude if drive >= 0.0 else -magnitude
            state[j, 12] = state[j, 5]
            state[j, 5] = state[j - 1, 8] + z * flow
            state[j - 1, 6] = state[j - 1, 7] - z * flow
            if jet is not None:
                state[j - 1, 14] = flow
        _lining_at(walls, gp, walls_offset, j, q, v, bq, bv)
        yl = 0.5 * (areas[g, j - 1] + areas[g, j]) / RHO_C
        if j < last:
            yr = 0.5 * (areas[g, j] + areas[g, j + 1]) / RHO_C
            admittance = yl + yr
            known = 2.0 * yl * state[j, 5] + 2.0 * yr * state[j, 6] - frame[g, j]
        else:
            resistance_end, gain = _radiation(areas[g, last], end_extra)
            admittance = yl + 1.0 / resistance_end + gain
            known = 2.0 * yl * state[j, 5] - state[0, 9] - gain * state[0, 10] - frame[g, j]
        h = 1.0 / 96000.0
        mc = walls[g, j, walls_offset + 1]
        if mc > 0.0:
            surface = walls[g, j, walls_offset]
            cc = walls[g, j, walls_offset + 6]
            resist = surface * surface / admittance + cc
            push = surface * known / admittance + cc * bv[j]
            following = ((mc - 0.5 * h * resist) * v[j] + h * push) / (mc + 0.5 * h * resist)
            middle = 0.5 * (v[j] + following)
            v[j] = following
            pressure = (known - surface * middle) / admittance
        else:
            pressure = known / admittance
        state[j, 11] = state[j, 4]
        state[j, 4] = pressure
        if j == last:
            resistance_end, gain = _radiation(areas[g, last], end_extra)
            state[0, 9] += gain * (pressure + state[0, 10])
            state[0, 10] = pressure
            lips = pressure / resistance_end + state[0, 9]
        else:
            state[j, 7] = (1.0 - hop) * (pressure - state[j, 6]) + hop * state[j, 7]
        state[j - 1, 13] = state[j - 1, 8]
        state[j - 1, 8] = (1.0 - hop) * (pressure - state[j, 5]) + hop * state[j - 1, 8]
    return arriving_port, lips


@njit(cache=True)
def _sealed_front(areas, t, port_node, seal):
    if seal <= 0.0:
        return False
    g = min(max(t, 0), areas.shape[0] - 1)
    for node in range(port_node + 1, areas.shape[1]):
        if areas[g, node] <= seal:
            return True
    return False


@njit(cache=True)
def _sealed_counts(areas, port_node, seal):
    """How many of the steps before each have a seal in front of the port
    (_sealed_front per step, summed: entry g + 1 counts steps 0..g)."""
    rows = areas.shape[0]
    out = numpy.zeros(rows + 1, numpy.int64)
    for g in range(rows):
        out[g + 1] = out[g] + (1 if _sealed_front(areas, g, port_node, seal) else 0)
    return out


@njit(cache=True)
def _peel_from_front(front, oral, port_node):
    """The time-major oral peel's state from the front part's forward state
    (their nodes at the same times: before a diagonal the peel's node j sits
    where the front's does, port_node + i at the front's time + i)."""
    n = oral.shape[0]
    for j in range(port_node, n):
        for c in range(4):
            oral[j, c] = front[j, c]
        oral[j, 4] = front[j, 4]
        oral[j, 6] = front[j, 11]
        oral[j, 5] = front[j, 5]
        oral[j, 7] = front[j, 12]
        if j < n - 1:
            oral[j, 8] = front[j, 13]
            oral[j, 9] = front[j, 7]
            oral[j, 10] = front[j, 6]
    oral[0, 11] = front[0, 9]
    oral[0, 12] = front[0, 10]
    for j in range(port_node, n - 1):
        oral[j, 15] = front[j, 14]
        for c in range(4):
            oral[j, 16 + c] = front[j, 15 + c]
        oral[j, 20] = front[j, 19]


@njit(cache=True)
def _march_run(t_begin, t_end, sound, areas, walls, frame, resistance, hop, end_extra, port, port_node, nose_areas,
               nose_walls, nose_targets, glide_start, nose_hop, nose_hop_zero, seal, oral, front, ahead_state, forward,
               peel, share,
               lips, returned, port_pressure, from_front, mode, leaving, delivered_out, front_lips, scalars, override,
               pressure_guess, window_waves, in_window, sealed_count, walls_offset, jet=None):
    """Steps t_begin..t_end - 1 of the split (see march_split), every state
    passed in and updated. scalars: 0 the nose's last leaving wave (forward),
    1 the ahead copy's base time. override[k], where finite, is taken as the
    nostrils' sample at k in the sealed mode; pressure_guess[t], where finite,
    as the port node's pressure at t for the part in front of the seal run
    ahead (a release window's unknowns: past the release the lips hear it
    before the march has reached it). Where in_window[t], the step is taken
    from window_waves (sustain PORTWAVES of a render whose folds run free
    there: the sound cannot decide those steps, our folds do): the port's
    waves, the port node's pressure, its outer link's wave and the
    nostrils' flow drive every part's state through. `jet`: the vocal tube's
    reversible jet (_diagonal's), or None. sealed_count: _sealed_counts'
    (whether a seal lies in front of the port, per step, summed); walls:
    the linings from column walls_offset (the geometry itself with 3)."""
    count = sound.shape[0]
    n = areas.shape[1]
    m = nose_areas.shape[0]
    port_y = port / RHO_C
    nose_areas_rows = numpy.zeros((1, m))
    nose_areas_rows[0, :] = nose_areas
    nose_walls_rows = numpy.zeros((1, m, 9))
    nose_walls_rows[0, :, :] = nose_walls
    nose_frame = numpy.zeros((1, m))
    rows = numpy.zeros((17, m - 1))
    for r in range(17):
        for j in range(m - 1):
            rows[r, j] = min((glide_start + r + 1) / 16.0, 1.0) * nose_targets[j]
    no_side = numpy.zeros(1)
    oral_lag = n - 1 - port_node
    nose_leaving = scalars[0]
    ahead_base = int(scalars[1])
    for t in range(t_begin, t_end):
        k = t + m - 1
        ko = t - 1 + oral_lag
        # The oral peel is good while no seal falls in its reach (its port
        # node at t - 1, its open end at t - 1 + oral_lag); past that, the
        # nostrils carry the sound and the lips' own flow comes from the
        # part in front of the seal, run on alone (sealed, it has no input).
        # (Whether any step in that reach has a seal in front of the port:
        # from the per-step counts, not a scan of reach x nodes per step.)
        lo = max(t - 1, 0)
        hi = min(t - 1 + oral_lag, count - 1)
        crossing = False
        if hi >= lo:
            bottom = min(lo, areas.shape[0] - 1)
            top = min(hi, areas.shape[0] - 1)
            crossing = sealed_count[top + 1] - sealed_count[bottom] > 0
        sealed = crossing and port[t] >= 1e-6 and t >= 1
        windowed = in_window[t]
        if sealed and not windowed and mode[t - 1] != 1.0:
            # (Entered from the open mode, or from a window the folds decided:
            # the front part is current in both.)
            ahead_state[:, :] = front[:, :]
            ahead_base = t - 2
        if not sealed and not windowed and t >= 1 and mode[t - 1] != 0.0:
            _peel_from_front(front, oral, port_node)
        if windowed:
            mode[t] = 2.0
            delivered = window_waves[t, 0]
            returned[t] = window_waves[t, 1]
            if t >= 1:
                port_pressure[t - 1] = window_waves[t - 1, 2]
                from_front[t - 1] = window_waves[t - 2, 3] if t >= 2 else 0.0
                arrived, lips_front = _oral_forward(front, t - 1, port_pressure[t - 1], port_node, areas, walls,
                                                    frame, resistance, hop, end_extra, seal, -1, jet, walls_offset)
                if t - 1 + oral_lag < front_lips.shape[0]:
                    front_lips[t - 1 + oral_lag] = lips_front
            if k < count:
                share[k] = window_waves[k, 4]
            nose_leaving, radiated = _nose_forward(forward, t, delivered, port[t] if port[t] >= 1e-6 else 0.0,
                                                   nose_areas, nose_walls, nose_targets, glide_start, nose_hop, 1.0,
                                                   1.0, nose_hop_zero)
            _diagonal(peel, k, share, nose_areas_rows, nose_walls_rows, nose_frame, rows, nose_hop, 0.0, port_y, 0,
                      nose_hop_zero)
            leaving[t] = peel[0, 4] - delivered
            delivered_out[t] = delivered
            continue
        owner = areas[min(max(t - 1, 0), areas.shape[0] - 1), port_node]
        delivered = 0.0
        if sealed:
            mode[t] = 1.0
            lips_now = 0.0
            while ahead_base + oral_lag < k:
                ahead_base += 1
                if ahead_base < pressure_guess.shape[0] and numpy.isfinite(pressure_guess[ahead_base]):
                    held = pressure_guess[ahead_base]
                else:
                    held = port_pressure[min(ahead_base, t - 2)] if t >= 2 else 0.0
                arrived, lips_now = _oral_forward(ahead_state, ahead_base, held, port_node, areas, walls, frame,
                                                  resistance, hop, end_extra, seal, port_node + 1, jet, walls_offset)
            if k < override.shape[0] and numpy.isfinite(override[k]):
                share[k] = override[k]
            else:
                share[k] = sound[min(k, count - 1)] - lips_now
            left_over = _diagonal(peel, k, share, nose_areas_rows, nose_walls_rows, nose_frame, rows, nose_hop, 0.0,
                                  port_y, 0, nose_hop_zero)
            delivered = left_over / (2.0 * port_y[t])
            sent_before, back = _port_step(port[t], delivered, leaving[t - 1], owner, nose_areas[0])
            returned[t] = back
            leaving[t] = peel[0, 4] - delivered
            port_pressure[t - 1] = sent_before + returned[t - 1]
            from_front[t - 1], lips_front = _oral_forward(front, t - 1, port_pressure[t - 1], port_node, areas, walls,
                                                          frame, resistance, hop, end_extra, seal, -1, jet,
                                                          walls_offset)
            if t - 1 + oral_lag < front_lips.shape[0]:
                front_lips[t - 1 + oral_lag] = lips_front
            nose_leaving, radiated = _nose_forward(forward, t, delivered, port[t], nose_areas, nose_walls,
                                                   nose_targets, glide_start, nose_hop, 1.0, 1.0, nose_hop_zero)
        else:
            mode[t] = 0.0
            if 0 <= ko < count:
                lips[ko] = sound[ko] - share[ko]
                _diagonal(oral, ko, lips, areas, walls, frame, resistance, hop, end_extra, no_side, port_node, 0.0, jet,
                          walls_offset)
            side_out = 0.0
            if t >= 1 and ko >= 0:
                port_pressure[t - 1] = oral[port_node, 4]
                from_front[t - 1] = oral[port_node, 11]
                if port[t - 1] >= 1e-6:
                    side_out = port_pressure[t - 1] - returned[t - 1]
                arrived, lips_front = _oral_forward(front, t - 1, port_pressure[t - 1], port_node, areas, walls, frame,
                                                    resistance, hop, end_extra, seal, -1, jet, walls_offset)
                if t - 1 + oral_lag < front_lips.shape[0]:
                    front_lips[t - 1 + oral_lag] = lips_front
            if port[t] >= 1e-6:
                z = RHO_C / port[t]
                r_port = max(0.03 * max(port[t], 1e-4) ** -2.6, 0.0065 / max(port[t], 1e-4) * numpy.sqrt(2250.0)) * NODE_LENGTH
                drive = 2.0 * (side_out - nose_leaving)
                if drive >= 0.0:
                    narrow, wide = owner, nose_areas[0]
                else:
                    narrow, wide = nose_areas[0], owner
                kinetic = 0.0
                if wide > narrow:
                    difference = 1.0 / max(narrow, 1e-4) - 1.0 / wide
                    kinetic = 0.5 * 1.14e-3 * difference * difference
                linear = 2.0 * z + r_port
                magnitude = 2.0 * abs(drive) / (linear + numpy.sqrt(linear * linear + 4.0 * kinetic * abs(drive)))
                flow = magnitude if drive >= 0.0 else -magnitude
                delivered = nose_leaving + z * flow
                returned[t] = side_out - z * flow
            nose_leaving, radiated = _nose_forward(forward, t, delivered, port[t] if port[t] >= 1e-6 else 0.0,
                                                   nose_areas, nose_walls, nose_targets, glide_start, nose_hop, 1.0,
                                                   1.0, nose_hop_zero)
            share[k] = radiated
            _diagonal(peel, k, share, nose_areas_rows, nose_walls_rows, nose_frame, rows, nose_hop, 0.0, port_y, 0,
                      nose_hop_zero)
            leaving[t] = peel[0, 4] - delivered
        delivered_out[t] = delivered
    scalars[0] = nose_leaving
    scalars[1] = ahead_base


# The vocal tube below the port, peeled in skewed chunks (as
# reverse_core.SidePeel: node j's values at t need node j + 1's to t + 1
# only, so a chunk's boundary carries a few in-flight samples and each
# link's and node's state; the values are the whole peel's to the bit).
# Per link (BELOW_LINK): 0 smoothed wave, 1 last flow, 2-5 jet filter.
BELOW_LINK = 6
BELOW_CHUNK = 512


@njit(cache=True)
def _below_top(geometry, port_node, pressure, from_front, returned, port_y, arriving, lining, t0, t1):
    """The port node's junction for steps t0..t1 - 1 (its pressure given)."""
    q = lining[0]
    v = lining[1]
    bq = lining[2]
    bv = lining[3]
    j = port_node
    for t in range(t0, t1):
        _lining_at(geometry, max(t - 1, 0), 3, j, q, v, bq, bv)
        middle = _middle_at(geometry, t, 3, j, pressure[t], v, bv)
        right = 0.5 * (geometry[t, j, 0] + geometry[t, j + 1, 0]) / RHO_C
        left = 0.5 * (geometry[t, j - 1, 0] + geometry[t, j, 0]) / RHO_C
        known = pressure[t] * (left + right + port_y[t]) + geometry[t, j, 3] * middle - 2.0 * port_y[t] * returned[t]
        arriving[t] = (known + geometry[t, j, 1] - 2.0 * right * from_front[t]) / (2.0 * left)


@njit(cache=True)
def _below_link(geometry, hop, kinetic, jet, pressure, arriving, rightward_in, leftward_out, lining, links, j, a0, a1,
                p0, p1, count):
    """Link j over a0..a1 - 1, node j over p0..p1 - 1 (_peel_below's loops)."""
    q = lining[0]
    v = lining[1]
    bq = lining[2]
    bv = lining[3]
    smooth = links[j, 0]
    previous_flow = links[j, 1]
    for t in range(a0, a1):
        raw = (pressure[t - 1] - arriving[t - 1]) if t > 0 else 0.0
        smooth = (1.0 - hop) * raw + hop * smooth
        made = max(t - 1, 0)
        first = geometry[made, j, 0]
        second = geometry[made, j + 1, 0]
        z = RHO_C / (0.5 * (first + second))
        flow = (arriving[t] - smooth) / z
        entered = smooth + (z + 0.5 * geometry[made, j, 2]) * flow
        if kinetic:
            entered += _kinetic(flow, first, second)
        if jet is not None:
            if t > 0:
                row = _row(jet[4], made)
                entered -= 0.5 * kernels_jet.jet_source(j, previous_flow, jet[3][row], jet[4][row], links, j, 2,
                                                        jet[7], jet[0] * jet[5][row], jet[1], made, jet[2], jet[6], jet[8])
            previous_flow = flow
        rightward_in[t] = entered
        leftward_out[t] = entered - z * flow
    links[j, 0] = smooth
    links[j, 1] = previous_flow
    for t in range(p0, p1):
        following = rightward_in[t + 1] if t + 1 < count else 0.0
        pressure[t] = (following - hop * rightward_in[t]) / (1.0 - hop) + leftward_out[t]
    for t in range(p0, p1):
        _lining_at(geometry, max(t - 1, 0), 3, j, q, v, bq, bv)
        middle = _middle_at(geometry, t, 3, j, pressure[t], v, bv)
        right = 0.5 * (geometry[t, j, 0] + geometry[t, j + 1, 0]) / RHO_C
        if j > 0:
            left = 0.5 * (geometry[t, j - 1, 0] + geometry[t, j, 0]) / RHO_C
            known = pressure[t] * (left + right) + geometry[t, j, 3] * middle
            arriving[t] = (known + geometry[t, j, 1] - 2.0 * right * leftward_out[t]) / (2.0 * left)
        else:
            arriving[t] = pressure[t] * right + geometry[t, 0, 3] * middle + geometry[t, 0, 1] - 2.0 * right * leftward_out[t]


@njit(cache=True)
def _below_snapshot(index, boundary, pressure, arriving, rightward_in, leftward_out, lining, links, snap_links,
                    snap_lining, snap_waves):
    snap_links[index, :, :] = links
    snap_lining[index, :, :] = lining
    span = snap_waves.shape[2]
    for c in range(span):
        t = boundary - span + c
        if t >= 0:
            snap_waves[index, 0, c] = pressure[t]
            snap_waves[index, 1, c] = arriving[t]
            snap_waves[index, 2, c] = rightward_in[t]
            snap_waves[index, 3, c] = leftward_out[t]


@njit(cache=True)
def _below_advance(geometry, port_node, from_front, returned, port_y, hop, kinetic, jet, pressure, arriving,
                   rightward_in, leftward_out, lining, links, begin, target, chunk, snap_every, snap_links, snap_lining,
                   snap_waves):
    """Chunks from boundary `begin` to `target` (count: complete): the port
    node, then each link and node down to the glottis."""
    count = pressure.shape[0]
    last_link = port_node - 1
    boundary = begin
    while boundary < target:
        if snap_every > 0 and boundary % snap_every == 0 and boundary > begin:
            _below_snapshot(boundary // snap_every, boundary, pressure, arriving, rightward_in, leftward_out, lining,
                            links, snap_links, snap_lining, snap_waves)
        following = min(boundary + chunk, target)
        final = following >= count
        _below_top(geometry, port_node, pressure, from_front, returned, port_y, arriving, lining, boundary, following)
        for j in range(last_link, -1, -1):
            lag = last_link - j
            a0 = max(boundary - lag, 0)
            a1 = count if final else following - lag
            if a1 <= a0:
                break
            p0 = max(a0 - 1, 0)
            p1 = count if final else a1 - 1
            _below_link(geometry, hop, kinetic, jet, pressure, arriving, rightward_in, leftward_out, lining, links, j,
                        a0, a1, p0, p1, count)
        boundary = following
    if snap_every > 0 and boundary % snap_every == 0 and boundary > begin and boundary < count:
        _below_snapshot(boundary // snap_every, boundary, pressure, arriving, rightward_in, leftward_out, lining,
                        links, snap_links, snap_lining, snap_waves)


class BelowPeel:
    """_peel_below's problem run in skewed chunks; with `record`, its state
    kept at every chunk boundary so that a peel whose inputs (port node's
    pressure, its outer link's wave, the port's returned wave) first differ
    at step k resumes from the last boundary at or before k (resumed). The
    geometry, port admittance and jet are the same for the resumed one."""

    def __init__(self, geometry, port_node, port_pressure, from_front, returned, port_y, hop, kinetic, jet=None,
                 record=False, chunk=BELOW_CHUNK):
        self.geometry = numpy.ascontiguousarray(geometry, float)
        self.port_node = int(port_node)
        self.port_y = numpy.ascontiguousarray(port_y, float)
        self.hop = float(hop)
        self.kinetic = bool(kinetic)
        self.jet = jet
        self.chunk = int(chunk)
        self.from_front = numpy.ascontiguousarray(from_front, float)
        self.returned = numpy.ascontiguousarray(returned, float)
        count = len(port_pressure)
        self.count = count
        n = self.geometry.shape[1]
        self.pressure = numpy.array(port_pressure, float)
        self.arriving = numpy.zeros(count)
        self.rightward_in = numpy.zeros(count)
        self.leftward_out = numpy.zeros(count)
        self.lining = numpy.zeros((4, n))
        self.links = numpy.zeros((max(self.port_node, 1), BELOW_LINK))
        self.boundary = 0
        slots = (count + self.chunk - 1) // self.chunk + 1 if record else 0
        self.snap_every = self.chunk if record else 0
        self.snap_links = numpy.zeros((slots, self.links.shape[0], BELOW_LINK))
        self.snap_lining = numpy.zeros((slots, 4, n))
        self.snap_waves = numpy.zeros((slots, 4, self.port_node + 2 if record else 1))

    def advance(self, target):
        """Carry the peel to `target` (count: complete); node 0's flow is
        then final below target - port_node (all, at count)."""
        target = min(int(target), self.count)
        if target <= self.boundary:
            return
        if target < self.count:
            target = min(max(self.boundary, (target + self.chunk - 1) // self.chunk * self.chunk), self.count)
        _below_advance(self.geometry, self.port_node, self.from_front, self.returned, self.port_y, self.hop,
                       self.kinetic, self.jet, self.pressure, self.arriving, self.rightward_in, self.leftward_out,
                       self.lining, self.links, self.boundary, target, self.chunk, self.snap_every, self.snap_links,
                       self.snap_lining, self.snap_waves)
        self.boundary = target

    def resumed(self, port_pressure, from_front, returned, first_changed, geometry=None, port_y=None, jet=None):
        """A peel for inputs equal to this one's before step first_changed,
        from the last recorded boundary at or before it (this one recorded).
        geometry, port_y, jet: another hypothesis's, equal to these before
        first_changed (None: these)."""
        other = BelowPeel.__new__(BelowPeel)
        other.__dict__.update(self.__dict__)
        other.from_front = numpy.ascontiguousarray(from_front, float)
        other.returned = numpy.ascontiguousarray(returned, float)
        if geometry is not None:
            other.geometry = numpy.ascontiguousarray(geometry, float)
        if port_y is not None:
            other.port_y = numpy.ascontiguousarray(port_y, float)
        if jet is not None:
            other.jet = jet
        index = min(max(int(first_changed), 0), self.count, self.boundary) // self.chunk
        index = min(index, (self.count - 1) // self.chunk)
        boundary = index * self.chunk
        other.boundary = boundary
        other.pressure = self.pressure.copy()
        other.arriving = self.arriving.copy()
        other.rightward_in = self.rightward_in.copy()
        other.leftward_out = self.leftward_out.copy()
        span = self.snap_waves.shape[2]
        for column in range(span):
            t = boundary - span + column
            if t >= 0:
                other.pressure[t] = self.snap_waves[index, 0, column]
                other.arriving[t] = self.snap_waves[index, 1, column]
                other.rightward_in[t] = self.snap_waves[index, 2, column]
                other.leftward_out[t] = self.snap_waves[index, 3, column]
        # From the boundary on, the slots hold the port node's own pressure
        # until the link below replaces them.
        other.pressure[boundary:] = numpy.asarray(port_pressure, float)[boundary:]
        other.links = self.snap_links[index].copy()
        other.lining = self.snap_lining[index].copy()
        other.snap_links = self.snap_links.copy()
        other.snap_lining = self.snap_lining.copy()
        other.snap_waves = self.snap_waves.copy()
        return other


class BelowMemo:
    """The last _peel_below's inputs and recorded run, to resume the next
    call from where its inputs first differ (march_residual's windows each
    peel the same geometry with the march a little further on). Used only
    for a read-only geometry (moving_geometry's), the very same object:
    its values cannot have changed in place."""

    def __init__(self, geometry, port_node, port_pressure, from_front, returned, port_y, hop, kinetic, jet, run):
        self.geometry = weakref.ref(geometry)
        self.port_node = port_node
        self.port_pressure = numpy.array(port_pressure, float)
        self.from_front = numpy.array(from_front, float)
        self.returned = numpy.array(returned, float)
        self.port_y = numpy.array(port_y, float)
        self.hop = hop
        self.kinetic = kinetic
        self.jet = jet
        self.jet_arrays = None if jet is None else tuple(numpy.array(jet[index]) for index in (3, 4, 5))
        self.run = run

    def matches(self, geometry, port_node, port_pressure, port_y, hop, kinetic, jet) -> bool:
        if self.geometry() is not geometry or geometry.flags.writeable:
            return False
        if port_node != self.port_node or hop != self.hop or kinetic != self.kinetic:
            return False
        if len(port_pressure) != len(self.port_pressure) or not _same_bits(port_y, self.port_y):
            return False
        if (jet is None) != (self.jet is None):
            return False
        if jet is not None:
            if jet is not self.jet:
                return False
            for index, kept in zip((3, 4, 5), self.jet_arrays):
                if not _same_bits(jet[index], kept):
                    return False
        return True


def _same_bits(first, second) -> bool:
    a = numpy.ascontiguousarray(first)
    b = numpy.ascontiguousarray(second)
    return a.shape == b.shape and a.dtype == b.dtype and bool(numpy.array_equal(a.view(numpy.uint8), b.view(numpy.uint8)))


def _same_or_contiguous(array):
    """The array itself when already C-contiguous float64, else a contiguous
    float64 copy. (numpy.ascontiguousarray(a, float) of a numba-made array
    returns a new view object every call, so the memo's identity test never
    held: measured, no resume across the window renders.)"""
    if isinstance(array, numpy.ndarray) and array.dtype == numpy.float64 and array.flags.c_contiguous:
        return array
    return numpy.ascontiguousarray(array, float)


_BELOW_MEMO = []   # at most one BelowMemo (the last call's)


def _peel_below(geometry, port_node, port_pressure, from_front, returned, port_y, hop, kinetic, jet=None, until=None):
    """The vocal tube peeled from the port node to the glottis, the port
    node's pressure, the wave its outer link sent into it and the port's
    returned wave given per step. Returns the flow entering. `jet`: the
    reversible jet (jet_table's), or None. With `until`, only the flow and
    pressure before step `until` are wanted: the peel stops a chunk past
    what they need, and later entries are NOT the whole peel's (a window's
    render reads only those before it). (Skewed chunks; on the same
    read-only geometry as the call before, resumed from where the inputs
    first differ: the whole peel's values to the bit.)"""
    geometry = _same_or_contiguous(geometry)
    hop = float(hop)
    kinetic = bool(kinetic)
    port_node = int(port_node)
    memo = _BELOW_MEMO[0] if _BELOW_MEMO else None
    run = None
    if memo is not None and memo.matches(geometry, port_node, port_pressure, port_y, hop, kinetic, jet):
        first_changed = min(reverse_core.first_difference(port_pressure, memo.port_pressure),
                            reverse_core.first_difference(from_front, memo.from_front),
                            reverse_core.first_difference(returned, memo.returned))
        run = memo.run.resumed(port_pressure, from_front, returned, first_changed, geometry)
    if run is None:
        run = BelowPeel(geometry, port_node, port_pressure, from_front, returned, port_y, hop, kinetic, jet,
                        record=not geometry.flags.writeable)
    if until is None:
        run.advance(run.count)
    else:
        run.advance(max(int(until), 0) + port_node + 1)
    _BELOW_MEMO[:] = []
    arriving = run.arriving.copy()
    pressure = run.pressure.copy()
    if not geometry.flags.writeable:
        # (The memo holds the geometry weakly only: kept alive by it, the
        # last word's 170-460 MB geometry would outlive its residual.)
        run.geometry = None
        _BELOW_MEMO.append(BelowMemo(geometry, port_node, port_pressure, from_front, returned, port_y, hop, kinetic,
                                     jet, run))
    return arriving, pressure



class March:
    """The split's whole state (march_split), so a stretch can be run again
    from a snapshot."""

    def __init__(self, count, n, m):
        self.arrays = {
            "oral": numpy.zeros((n, COLUMNS)), "front": numpy.zeros((n, ORAL_FORWARD_COLUMNS)),
            "ahead_state": numpy.zeros((n, ORAL_FORWARD_COLUMNS)), "forward": numpy.zeros((m, NOSE_COLUMNS)),
            "peel": numpy.zeros((m, COLUMNS)), "share": numpy.zeros(count + m), "lips": numpy.zeros(count),
            "returned": numpy.zeros(count), "port_pressure": numpy.zeros(count), "from_front": numpy.zeros(count),
            "mode": numpy.zeros(count), "leaving": numpy.zeros(count + 1), "delivered_out": numpy.zeros(count),
            "front_lips": numpy.zeros(count + n), "scalars": numpy.zeros(2)}

    def snapshot(self):
        copy = March.__new__(March)
        copy.arrays = {key: value.copy() for key, value in self.arrays.items()}
        return copy

    def restore(self, other):
        for key, value in other.arrays.items():
            self.arrays[key][...] = value

    def small(self):
        """Copies of the small per-node states (SNAPSHOT_KEYS)."""
        return {key: self.arrays[key].copy() for key in SNAPSHOT_KEYS}

    def restore_from(self, arrays, small, step, oral_lag, m):
        """The state before `step` of a run whose final arrays are `arrays`
        and whose small states there are `small`: each per-step array's
        entries written before `step` (final by then) kept, the rest zero as
        a run reaching `step` would have them."""
        for key, value in arrays.items():
            self.arrays[key][...] = value
        for key, value in small.items():
            self.arrays[key][...] = value
        frontier = {"port_pressure": step - 1, "from_front": step - 1, "returned": step, "mode": step,
                    "leaving": step, "delivered_out": step, "share": step + m - 1, "lips": step - 1 + oral_lag,
                    "front_lips": step - 1 + oral_lag}
        for key, written in frontier.items():
            self.arrays[key][max(written, 0):] = 0.0


ORDER = ("oral", "front", "ahead_state", "forward", "peel", "share", "lips", "returned", "port_pressure",
         "from_front", "mode", "leaving", "delivered_out", "front_lips", "scalars")


def march_split(sound, geometry, port, nose, hop, end_extra, seal, newton=4, nose_hop=None, nose_hop_zero=0.0):
    """The nostrils' share of the sound, the port's returned wave, the port
    node's pressure and the wave its outer link sent into it, per step.
    Open, the oral path is peeled and the port and nose run forward; sealed,
    the nose is peeled. At each release the lips hear the port node's
    pressure (through the part in front of the seal) up to ~9 steps before
    the march reaches it: those pressures across the release window are
    solved together (Newton), each guess against what the march then
    computes."""
    port_node, nose_areas, targets, nose_walls, dumped_hop, scale, first = nose
    pole = dumped_hop if nose_hop is None else nose_hop
    glide_start = int(round(16.0 * float(numpy.mean(first / targets)))) - 1
    geometry = numpy.ascontiguousarray(geometry, float)
    areas = numpy.ascontiguousarray(geometry[:, :, 0])
    walls = geometry   # (its linings read in place from column 3: no copy of three quarters of it)
    frame = numpy.ascontiguousarray(geometry[:, :, 1])
    resistance = numpy.ascontiguousarray(geometry[:, :-1, 2])
    nose_areas = numpy.ascontiguousarray(nose_areas)
    nose_walls = numpy.ascontiguousarray(nose_walls)
    count = len(sound)
    n = areas.shape[1]
    m = len(nose_areas)
    oral_lag = n - 1 - port_node
    sealed_count = _sealed_counts(areas, int(port_node), float(seal))
    state = March(count, n, m)
    override = numpy.full(count + m, numpy.nan)
    guess_array = numpy.full(count, numpy.nan)
    fixed = (sound, areas, walls, frame, resistance, hop, end_extra, port, port_node, nose_areas, nose_walls, targets,
             glide_start, pole, nose_hop_zero, seal)

    no_waves = numpy.zeros((count, 6))
    no_window = numpy.zeros(count, bool)

    def run(begin, end):
        _march_run(begin, end, *fixed, *[state.arrays[key] for key in ORDER], override, guess_array, no_waves,
                   no_window, sealed_count, 3)

    sealed = (areas[:, port_node + 1:].min(axis=1) <= seal) if seal > 0.0 else numpy.zeros(count, bool)
    releases = [t for t in range(1, count) if sealed[t - 1] and not sealed[t]]
    t = 0
    for release in releases:
        last = release - 1                        # the last sealed step
        begin = max(last - (m - 1) + 1, t)        # the first step whose nostril sample passes it
        # The lips hear the port node oral_lag steps late; the last nostril
        # sample the sealed mode takes is last + m, so the pressures wanted
        # run to last + m - oral_lag; the march computes them by its step
        # last + m - oral_lag + 1 (the oral peel's, after the release).
        reach = last + m - oral_lag
        end = min(reach + 2, count)
        unknown = numpy.arange(max(begin - 2, 0), min(reach + 1, count))
        run(t, begin)
        # The copy run ahead was driven by held pressures through the seal
        # (harmless to the sealed lips, wrong for the cavity behind them):
        # from the front part's exact state, guessed pressures from here.
        state.arrays["ahead_state"][...] = state.arrays["front"]
        state.arrays["scalars"][1] = begin - 2
        start = state.snapshot()
        known = state.arrays["port_pressure"][max(begin - 3, 0)]
        guess = numpy.full(len(unknown), known)

        def residual(values):
            state.restore(start)
            guess_array[:] = numpy.nan
            guess_array[unknown] = values
            run(begin, end)
            return state.arrays["port_pressure"][unknown] - values

        for round_index in range(newton):
            base = residual(guess)
            jacobian = numpy.zeros((len(unknown), len(unknown)))
            step = 1e-6 * max(float(numpy.max(numpy.abs(guess))), 1.0)
            for column in range(len(unknown)):
                moved = guess.copy()
                moved[column] += step
                jacobian[:, column] = (residual(moved) - base) / step
            guess = guess - numpy.linalg.solve(jacobian, base)
        state.restore(start)
        guess_array[:] = numpy.nan
        guess_array[unknown] = guess
        run(begin, end)
        guess_array[:] = numpy.nan
        t = end
    run(t, count)
    return (state.arrays["share"][:count].copy(), state.arrays["returned"].copy(),
            state.arrays["port_pressure"].copy(), state.arrays["from_front"].copy(), state.arrays["mode"].copy())



@njit(cache=True)
def _narrowest_table(areas, narrowest, least):
    for t in range(areas.shape[0]):
        narrowest[t], least[t] = reverse_core._narrowest(areas[max(t - 1, 0)])


def jet_table(areas, jet, gate=None):
    """The reversible jet for the marches: None without one, else (gain,
    seed, closing, narrowest node per step, its area per step, its gate per
    step, the white noise [link, step] (kernels_jet.noise_table), the
    filter coefficients' cache), each step's narrowing from the step
    before's areas (VoiceBody::reversible_jet's); `gate` per step (1
    without)."""
    if jet is None or jet[0] <= 0.0:
        return None
    areas = numpy.ascontiguousarray(areas, float)
    count = areas.shape[0]
    narrowest = numpy.zeros(count, numpy.int64)
    least = numpy.zeros(count)
    _narrowest_table(areas, narrowest, least)
    gates = numpy.ones(count)
    if gate is not None:
        size = min(count, len(gate))
        gates[:size] = gate[:size]
        gates[size:] = gate[size - 1] if size > 0 else 1.0
    # (The front runs ahead of the sound by up to a tube and a nose; past
    # the table the noise is hashed as before.)
    noise = kernels_jet.noise_table(int(jet[1]), areas.shape[1] - 1, count + 2 * areas.shape[1] + 64)
    # The ninth: the engine's JETBOUND (reverse_core.JET_BOUND; 0 off).
    return (float(jet[0]), int(jet[1]), float(jet[2]), narrowest, least, gates, noise, kernels_jet.new_cache(1)[0],
            float(reverse_core.JET_BOUND))


def march_windows(sound, geometry, port, nose, hop, end_extra, seal, windows, render_window, nose_hop, nose_hop_zero,
                  jet=None):
    """The split with our folds deciding the windows: the march up to each
    window, then render_window(state arrays, begin) -> PORTWAVES of a render
    forced with the waves implied so far and free from the window on; the
    window's steps taken from it. Returns (share, returned, port pressure,
    from_front, mode, window waves, window mask)."""
    return _march_windows_run(sound, geometry, port, nose, hop, end_extra, seal, windows, render_window, nose_hop,
                              nose_hop_zero, jet, False, None, 0)[:7]


MARCH_EVERY = 512   # steps between a recorded march's state snapshots
# The small per-node states a snapshot keeps (the per-step arrays are the
# run's outputs: their entries before a step are final by then).
SNAPSHOT_KEYS = ("oral", "front", "ahead_state", "forward", "peel", "scalars")


class MarchRecord:
    """A march_windows run kept for march_windows_resume: its windows, each
    window's render, the small states every MARCH_EVERY steps, its arrays."""

    def __init__(self, windows, renders, snapshots, arrays, lookahead, every):
        self.windows = windows
        self.renders = renders
        self.snapshots = snapshots
        self.arrays = arrays
        self.lookahead = lookahead
        self.every = every


def march_windows_record(sound, geometry, port, nose, hop, end_extra, seal, windows, render_window, nose_hop,
                         nose_hop_zero, jet=None):
    """march_windows, and its MarchRecord."""
    return _march_windows_run(sound, geometry, port, nose, hop, end_extra, seal, windows, render_window, nose_hop,
                              nose_hop_zero, jet, True, None, 0)


def march_windows_resume(record, sound, geometry, port, nose, hop, end_extra, seal, windows, render_window, nose_hop,
                         nose_hop_zero, jet=None, first_changed=0):
    """march_windows for inputs equal to the record's before step
    first_changed (sound, geometry rows, port, jet gate; and render_window
    causal: its waves before a step depend on its inputs and the hypothesis
    before it). The march resumes from the last snapshot a lookahead (the
    tube's and the nose's node counts) before first_changed; a window whose
    render could not have changed is taken from the record (render_window
    is not called for it). Returns march_windows' seven arrays and a new
    MarchRecord: the same to the bit as march_windows_record here."""
    return _march_windows_run(sound, geometry, port, nose, hop, end_extra, seal, windows, render_window, nose_hop,
                              nose_hop_zero, jet, True, record, int(first_changed))


def _march_windows_run(sound, geometry, port, nose, hop, end_extra, seal, windows, render_window, nose_hop,
                       nose_hop_zero, jet, keep, base, first_changed):
    port_node, nose_areas, targets, nose_walls, dumped_hop, scale, first = nose
    glide_start = int(round(16.0 * float(numpy.mean(first / targets)))) - 1
    geometry = numpy.ascontiguousarray(geometry, float)
    areas = numpy.ascontiguousarray(geometry[:, :, 0])
    walls = geometry   # (its linings read in place from column 3: no copy of three quarters of it)
    frame = numpy.ascontiguousarray(geometry[:, :, 1])
    resistance = numpy.ascontiguousarray(geometry[:, :-1, 2])
    nose_areas = numpy.ascontiguousarray(nose_areas)
    nose_walls = numpy.ascontiguousarray(nose_walls)
    count = len(sound)
    n = areas.shape[1]
    m = len(nose_areas)
    oral_lag = n - 1 - port_node
    windows = [(int(begin), int(end)) for begin, end in windows]
    sealed_count = _sealed_counts(areas, int(port_node), float(seal))
    state = March(count, n, m)
    override = numpy.full(count + m, numpy.nan)
    guess_array = numpy.full(count, numpy.nan)
    waves = numpy.zeros((count, 6))
    in_window = numpy.zeros(count, bool)
    fixed = (sound, areas, walls, frame, resistance, hop, end_extra, port, port_node, nose_areas, nose_walls, targets,
             glide_start, nose_hop, nose_hop_zero, seal)
    jet_steps = jet_table(areas, jet, None if jet is None or len(jet) < 4 else jet[3])
    # A step reads the sound, geometry, port and jet at most this far ahead
    # (the oral peel's open end, the nose peel's port admittance).
    lookahead = n + m + 2
    every = MARCH_EVERY
    snapshots = {}
    renders = {}
    marching_run = MarchRunner(fixed, state, override, guess_array, waves, in_window, sealed_count, jet_steps,
                               snapshots if keep else None, every)
    t = 0
    reuse = set()
    if base is not None:
        begin_at = max(first_changed - lookahead, 0) // every * every
        for index, window in enumerate(windows):
            if (index < len(base.windows) and base.windows[index] == window and
                    first_changed >= max(window[1] + m + 2, window[0] + lookahead)):
                reuse.add(index)
        moved = True
        while moved:
            moved = False
            for index, (begin, end) in enumerate(windows):
                if end > begin_at and index not in reuse and begin < begin_at:
                    begin_at = begin // every * every
                    moved = True
        begin_at = min(begin_at, max(base.snapshots) if base.snapshots else 0)
        while begin_at not in base.snapshots and begin_at > 0:
            begin_at -= every
        state.restore_from(base.arrays, base.snapshots[begin_at], begin_at, oral_lag, m)
        for boundary, kept in base.snapshots.items():
            if boundary <= begin_at:
                snapshots[boundary] = kept
        for index, (begin, end) in enumerate(windows):
            if begin < begin_at:
                dumped = base.renders[index]
                renders[index] = dumped
                size = min(len(dumped), count)
                waves[max(begin - 2, 0):size, :] = dumped[max(begin - 2, 0):size, :]
                in_window[begin:end] = True
        t = begin_at
    elif keep:
        snapshots[0] = state.small()
    for index, (begin, end) in enumerate(windows):
        if end <= t:
            continue
        if begin >= t:
            marching_run.run(t, begin)
            if index in reuse:
                dumped = base.renders[index]
            else:
                dumped = render_window(state.arrays, begin)
            renders[index] = dumped
            size = min(len(dumped), count)
            # (From two steps before: the window's first step reads the port
            # node's pressure of the step before and its outer wave of the step
            # before that. Copied from the window on, they read zero, and the
            # march after a window was lost: 44.4 dB where it had been 289.6.)
            waves[max(begin - 2, 0):size, :] = dumped[max(begin - 2, 0):size, :]
            in_window[begin:end] = True
            marching_run.run(begin, end)
        else:
            marching_run.run(t, end)
        t = end
    marching_run.run(t, count)
    outputs = (state.arrays["share"][:count].copy(), state.arrays["returned"].copy(),
               state.arrays["port_pressure"].copy(), state.arrays["from_front"].copy(), state.arrays["mode"].copy(),
               waves, in_window)
    record = None
    if keep:
        record = MarchRecord(windows, renders, snapshots, {key: value.copy() for key, value in state.arrays.items()},
                             lookahead, every)
    return outputs + (record,)


class MarchRunner:
    """_march_run over a stretch, split at every `every`-th step to keep the
    small states there (snapshots: None, none kept)."""

    def __init__(self, fixed, state, override, guess_array, waves, in_window, sealed_count, jet_steps, snapshots,
                 every):
        self.fixed = fixed
        self.state = state
        self.override = override
        self.guess_array = guess_array
        self.waves = waves
        self.in_window = in_window
        self.sealed_count = sealed_count
        self.jet_steps = jet_steps
        self.snapshots = snapshots
        self.every = every

    def _march(self, begin, end):
        _march_run(begin, end, *self.fixed, *[self.state.arrays[key] for key in ORDER], self.override,
                   self.guess_array, self.waves, self.in_window, self.sealed_count, 3, self.jet_steps)

    def run(self, begin, end):
        if end <= begin:
            return
        if self.snapshots is None:
            self._march(begin, end)
            return
        t = begin
        while t < end:
            boundary = min((t // self.every + 1) * self.every, end)
            self._march(t, boundary)
            t = boundary
            if t % self.every == 0 and t < len(self.state.arrays["mode"]):
                self.snapshots[t] = self.state.small()
