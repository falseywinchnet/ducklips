"""A control field as a Fourier feature pyramid in time (after Zhao et al.,
"Fourier Feature Pyramids for Physics-Informed Neural Networks",
arXiv 2605.24278), without its decoder: the field is the sum of its levels,
so it is linear in its coefficients and every solve is an exact least
squares.

Each level is a grid of samples over the word's ticks, coarse to fine, each
twice the last; a level is read at a tick by Fourier interpolation of its
grid. Time is not periodic, so the interpolation is the cosine one (the
grid's even extension: no seam at the word's ends). The finest level is the
1 ms tick grid itself (the engine reads its controls only at ticks), so the
pyramid spans every field the engine can be given: the band a word needs is
measured by how many levels it takes, never chosen.

Space is not a pyramid: a seal is one node wide (sequence_recover.places),
not bandlimited along the tube; each node is a channel."""
from typing import List

import numpy


def level_sizes(ticks: int, levels: int) -> List[int]:
    """The grid sizes coarse to fine: the finest the tick count, each coarser
    half the next (rounded up), never below one sample."""
    sizes: List[int] = []
    size = ticks
    for _ in range(levels):
        sizes.append(max(size, 1))
        size = (size + 1) // 2
    sizes.reverse()
    return sizes


def interpolation(ticks: int, size: int) -> numpy.ndarray:
    """The matrix (ticks x size) reading a level's grid of `size` samples at
    the word's ticks: the grid's cosine transform (DCT-II, orthonormal),
    then its modes evaluated at the ticks. Samples sit at the centres of
    `size` equal cells over the word, the ticks at the centres of `ticks`
    cells, so a level as fine as the ticks reads back its own samples."""
    cells = numpy.arange(size, dtype=numpy.float64)
    modes = numpy.arange(size, dtype=numpy.float64)
    # Orthonormal DCT-II of the grid: coefficient k from sample j.
    forward = numpy.cos(numpy.pi * numpy.outer(modes, cells + 0.5) / size)
    weights = numpy.full(size, numpy.sqrt(2.0 / size))
    weights[0] = numpy.sqrt(1.0 / size)
    forward *= weights[:, None]
    # The modes at the ticks (the same normalisation, inverse transform).
    positions = (numpy.arange(ticks, dtype=numpy.float64) + 0.5) / ticks
    read = numpy.cos(numpy.pi * numpy.outer(positions, modes)) * weights[None, :]
    return read @ forward


class Pyramid:
    """The levels' reading matrices for a word of `ticks` ticks."""

    def __init__(self, ticks: int, levels: int) -> None:
        self.ticks: int = ticks
        self.sizes: List[int] = level_sizes(ticks, levels)
        self.reads: List[numpy.ndarray] = [interpolation(ticks, size) for size in self.sizes]

    def field(self, grids: List[numpy.ndarray]) -> numpy.ndarray:
        """The field at the ticks (ticks x channels): the sum of the levels
        read (grids coarse to fine, each size x channels)."""
        total = numpy.zeros((self.ticks, grids[0].shape[1]))
        for read, grid in zip(self.reads, grids):
            total += read @ grid
        return total

    def coarse_to_fine(self, target: numpy.ndarray, used: int) -> List[numpy.ndarray]:
        """The first `used` levels' grids by exact least squares, coarse to
        fine: each level takes what the coarser ones left (the paper's
        coarse-to-fine order); the levels past `used` are zero."""
        grids: List[numpy.ndarray] = []
        left = numpy.asarray(target, dtype=numpy.float64).copy()
        for index, read in enumerate(self.reads):
            if index >= used:
                grids.append(numpy.zeros((self.sizes[index], left.shape[1])))
                continue
            grid = numpy.linalg.lstsq(read, left, rcond=None)[0]
            grids.append(grid)
            left -= read @ grid
        return grids
