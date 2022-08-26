"""Utility classes/functions."""
import numpy as np
from tqdm.auto import tqdm

import meep


class TqdmWrappedIterable:
    def __init__(self, iterable, **kwargs):
        self._iterable = iterable
        self._kwargs = kwargs

    def __iter__(self):
        yield from tqdm(self._iterable, **self._kwargs)

    def __len__(self):
        return len(self._iterable)


def _array_to_vec3s(arr):
    return [meep.Vector3(vec[0], vec[1], vec[2]) for vec in arr]


def interpolate_path(band_sequence, k_interpolation_points):
    band_points = []
    for start, end in zip(band_sequence[:-1], band_sequence[1:]):
        for t in np.linspace(0, 1, k_interpolation_points, endpoint=False):
            band_points.append(start * (1 - t) + end * t)
    band_points.append(end)
    return np.asarray(band_points)

def monkhorst_pack(size):
    """Uniform sampling of "cubic" Brillouin zone with `size[i]` samples along dimension `i`."""
    if np.less_equal(size, 0).any():
        raise ValueError('size entries must be non-negative: %s' % list(size))
    if len(size) == 1:
        idxs = np.array([np.array([i]) for i in range(0,size[0])])
    elif len(size) == 2:
        idxs = np.array([np.array([i,j]) for i in range(0,size[0]) for j in range(0,size[1])])
    elif len(size) == 3:
        idxs = np.array([np.array([i,j,k]) for i in range(0,size[0]) for j in range(0,size[1]) for k in range(0,size[2])])
    else:
        raise ValueError('unsupported dimension of size: %s' % len(size))
    return (idxs + 0.5) / size - 0.5