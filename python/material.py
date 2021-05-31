"""This file defines an example system whose density of states will be calculated with MPB.

This is analogous to the file ``3Dexample.ctl``.
"""
import numpy as np

import meep
from meep import mpb
from util import TqdmWrappedIterable, _array_to_vec3s, interpolate_path

basis_size = np.ones(3) * np.sqrt(3) / 2
basis = np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]])
basis1, basis2, basis3 = _array_to_vec3s(basis)
lattice = meep.Lattice(
    basis_size=basis_size, basis1=basis1, basis2=basis2, basis3=basis3
)

eps = 16
eps_background = 1
dielectric = meep.Medium(epsilon=eps)
air = meep.Medium(epsilon=eps_background)
isovalue = -1.1


def gyroid_material(h, k, l, x, y, z):
    # This is copied directly from a Scheme implementation
    prefactor = np.cos(2 * np.pi * (h + k + l) * 0.25)
    value_xyz = (
        np.sin(2 * np.pi * (h * x + l / 4))
        * np.sin(2 * np.pi * (k * y + h / 4))
        * np.sin(2 * np.pi * (l * z + k / 4))
    )
    value_yzx = (
        np.sin(2 * np.pi * (h * y + l / 4))
        * np.sin(2 * np.pi * (k * z + h / 4))
        * np.sin(2 * np.pi * (l * x + k / 4))
    )
    value_zxy = (
        np.sin(2 * np.pi * (h * z + l / 4))
        * np.sin(2 * np.pi * (k * x + h / 4))
        * np.sin(2 * np.pi * (l * y + k / 4))
    )
    return prefactor * (value_xyz + value_yzx + value_zxy)


def eps_func(p_lattice):
    p_cartesian = meep.lattice_to_cartesian(p_lattice, lattice)
    x, y, z = p_cartesian
    return (
        gyroid_material(1, 1, 0, x, y, z) < isovalue
        or gyroid_material(1, 1, 0, -x, -y, -z) < isovalue
    )


def default_material(p_lattice):
    in_dielectric = eps_func(p_lattice)
    return dielectric if in_dielectric else air


# We create a set of k points that will be evaluated
G = np.array([0, 0, 0])
H = np.array([0.5, -0.5, 0.5])
P = np.array([0.75, -0.25, -0.25])
N = np.array([0.5, 0, -0.5])

k_point_sequence = [H, G, N, P, G]
k_point_labels = ["H", "$\Gamma$", "N", "P", "$\Gamma$"]
k_point_interpolation = 10  # 20
k_point_band_path = _array_to_vec3s(
    interpolate_path(k_point_sequence, k_point_interpolation)
)
