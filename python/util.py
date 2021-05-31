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
