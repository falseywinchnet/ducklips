"""The reversible jet's source (VoiceBody::reversible_jet) as reverse_core's
_jet_source computes it, with its invariant work hoisted; bit-identical.

Two equivalences, both exact (the same expressions on the same inputs):
- the white noise is a pure function of (seed, link, step): a table filled
  once per process by the same _hashed_gauss serves every residual (the
  hypothesis never changes it);
- the shaping filter's coefficients are a pure function of the corner
  frequency, which sits on its 1500 Hz floor for nearly every link and step:
  they are recomputed only when the corner changes (per-link cache).
With the gain exactly zero (the gate shut, most of a word) and a finite
velocity, the source is (gain * 0.0) * state[3] / sqrt(passed), the
reference's value to the bit (strength = gain * dynamic * ... is a zero of
gain's sign when dynamic is finite); the Reynolds number, the two square
roots and the distance's exponential are then skipped. A non-finite or
huge velocity (a runaway hypothesis) takes the full expression."""
import numpy
from numba import njit

# A cache row per link: the corner it holds, its two filter gains, and the
# square root of the filter's passed power.
CACHE_COLUMNS = 4
VELOCITY_FINITE = 1e100   # below this, the dynamic pressure (velocity cubed) is finite


@njit(cache=True)
def _split_mix(key):
    key = key + numpy.uint64(0x9E3779B97F4A7C15)
    key = (key ^ (key >> numpy.uint64(30))) * numpy.uint64(0xBF58476D1CE4E5B9)
    key = (key ^ (key >> numpy.uint64(27))) * numpy.uint64(0x94D049BB133111EB)
    return key ^ (key >> numpy.uint64(31))


@njit(cache=True)
def hashed_gauss(seed, link, step):
    """reverse_core._hashed_gauss (the same expression)."""
    total = 0.0
    for draw in range(4):
        key = numpy.uint64(seed) ^ ((numpy.uint64(step) * numpy.uint64(4096) + numpy.uint64(link)) * numpy.uint64(4)
                                    + numpy.uint64(draw))
        total += float(_split_mix(key) >> numpy.uint64(11)) * (1.0 / 9007199254740992.0)
    return (total - 2.0) * numpy.sqrt(3.0)


@njit(cache=True)
def _fill_noise(seed, links, steps):
    out = numpy.empty((links, steps))
    for link in range(links):
        for step in range(steps):
            out[link, step] = hashed_gauss(seed, link, step)
    return out


_NOISE = {}


def noise_table(seed: int, links: int, steps: int) -> numpy.ndarray:
    """The jet's white noise [link, step] for at least `steps` steps (kept
    per process and seed; grown when a longer sound asks)."""
    key = (int(seed), int(links))
    table = _NOISE.get(key)
    if table is None or table.shape[1] < steps:
        size = max(int(steps), 0 if table is None else table.shape[1])
        table = _fill_noise(int(seed), int(links), size)
        _NOISE[key] = table
    return table


@njit(cache=True)
def new_cache(links):
    """A fresh per-link coefficient cache (NaN: nothing held)."""
    out = numpy.empty((links, CACHE_COLUMNS))
    out[:, :] = numpy.nan
    return out


@njit(cache=True)
def white_at(noise, seed, link, step):
    if link < noise.shape[0] and step < noise.shape[1]:
        return noise[link, step]
    return hashed_gauss(seed, link, step)


@njit(cache=True)
def jet_source(link, previous_flow, narrowest, area, filters, row, base, cache, gain, seed, step, closing, noise,
               bound=0.0, offset=0):
    """reverse_core._jet_source to the bit: its 4 filter values at
    filters[row, base:base + 4] (read in place: a view per call costs numba
    two atomic reference counts), `cache` a CACHE_COLUMNS row (any link's:
    it holds a pure function of the corner), `noise` noise_table's."""
    velocity = abs(previous_flow) / area
    if bound > 0.0:
        # The engine's JETBOUND: the fastest the lungs can drive air.
        velocity = min(velocity, bound)
    slit = max(area / 1.2, 0.01)
    corner = min(max(0.2 * velocity / slit, 1500.0), 12000.0)
    if corner == cache[0]:
        corner_gain = cache[1]
        floor_gain = cache[2]
        root_passed = cache[3]
    else:
        corner_gain = 1.0 - numpy.exp(-2.0 * numpy.pi * corner / 96000.0)
        floor_gain = 1.0 - numpy.exp(-2.0 * numpy.pi * corner / (3.0 * 96000.0))
        pole = 1.0 - corner_gain
        gain_squared = corner_gain * corner_gain
        pole_term = 1.0 - pole * pole
        passed = gain_squared * gain_squared * (1.0 + pole * pole) / (pole_term * pole_term * pole_term)
        root_passed = numpy.sqrt(passed)
        cache[0] = corner
        cache[1] = corner_gain
        cache[2] = floor_gain
        cache[3] = root_passed
    # `offset`: the link's global index less its local one (a peel of the
    # tube's front only); the noise, the narrowest node and the distance are
    # the engine's, in global numbering.
    whole = link + offset
    white = white_at(noise, seed, whole, step)
    filters[row, base] += floor_gain * (white - filters[row, base])
    once = white - filters[row, base]
    filters[row, base + 1] += floor_gain * (once - filters[row, base + 1])
    filters[row, base + 2] += corner_gain * ((once - filters[row, base + 1]) - filters[row, base + 2])
    filters[row, base + 3] += corner_gain * (filters[row, base + 2] - filters[row, base + 3])
    if narrowest < 0 or whole < narrowest:
        return 0.0
    if gain == 0.0 and velocity <= VELOCITY_FINITE:
        return (gain * 0.0) * filters[row, base + 3] / root_passed
    diameter = numpy.sqrt(4.0 * area / numpy.pi)
    reynolds = velocity * diameter * 1.14e-3 / 1.86e-4
    onset = 1800.0 / max(reynolds, 1.0)
    turbulent = max(1.0 - onset * onset, 0.0)
    fade = min(max((0.6 - area) / 0.3, 0.0), 1.0)
    if closing > 0.0:
        closed = 0.25 * closing
        fade *= min(max((area - closed) / (closing - closed), 0.0), 1.0)
    dynamic = 0.5 * 1.14e-3 * velocity * velocity * turbulent * fade * fade * (velocity / 2000.0)
    distance = (whole - narrowest) * (35000.0 / 96000.0)
    strength = gain * dynamic * numpy.sqrt(area) * numpy.exp(-distance / 1.0)
    return strength * filters[row, base + 3] / root_passed
