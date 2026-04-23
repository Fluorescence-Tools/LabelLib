import numpy as np
import pytest


ll = pytest.importorskip("LabelLib")


def _inputs():
    atoms = np.array(
        [
            [0.0, -4.0, 22.0, 1.5],
            [9.0, 0.0, 0.0, 3.0],
            [9.0, 8.0, 0.0, 3.0],
            [0.0, -10.5, 0.0, 1.5],
            [5.0, 0.0, 0.0, 2.0],
            [0.0, -4.0, -10.5, 1.3],
        ],
        dtype=np.float32,
    )
    source = np.array([0.0, -4.0, 0.0], dtype=np.float32)
    return atoms, source


def _checksum(grid):
    return float(np.sum(np.asarray(grid.grid, dtype=np.float64)))


def test_dye_density_arr_accepts_4xN_and_Nx4():
    atoms_nx4, source = _inputs()
    atoms_4xn = atoms_nx4.T

    g1 = ll.dyeDensityAV1_arr(atoms_4xn, source, 20.0, 2.0, 3.5, 0.9)
    g2 = ll.dyeDensityAV1_arr(atoms_nx4, source, 20.0, 2.0, 3.5, 0.9)

    assert np.isclose(_checksum(g1), _checksum(g2), rtol=0.0, atol=1e-4)


def test_add_weights_arr_accepts_5xN_and_Nx5():
    atoms_nx4, source = _inputs()
    av = ll.dyeDensityAV1_arr(atoms_nx4, source, 20.0, 2.0, 3.5, 0.9)

    xyzrq_nx5 = np.hstack([atoms_nx4[:2], np.full((2, 1), 2.0, dtype=np.float32)])
    xyzrq_5xn = xyzrq_nx5.T

    g1 = ll.addWeights_arr(av, xyzrq_5xn)
    g2 = ll.addWeights_arr(av, xyzrq_nx5)

    assert np.isclose(_checksum(g1), _checksum(g2), rtol=0.0, atol=1e-4)


def test_arr_distance_sampling_returns_numpy_vector_shape():
    atoms_nx4, source = _inputs()
    av = ll.dyeDensityAV1_arr(atoms_nx4, source, 20.0, 2.0, 3.5, 0.9)
    distances = ll.sampleDistanceDistInv_arr(av, av, 2000)

    assert isinstance(distances, np.ndarray)
    assert distances.ndim == 1
    assert distances.shape == (2000,)


def test_grid_array_helpers_shapes():
    atoms_nx4, source = _inputs()
    av = ll.dyeDensityAV1_arr(atoms_nx4, source, 20.0, 2.0, 3.5, 0.9)

    assert av.origin_arr().shape == (3,)
    assert av.shape_arr().shape == (3,)
    assert av.grid_arr().ndim == 1


def test_arr_functions_reject_wrong_shapes():
    atoms_nx4, source = _inputs()
    bad_atoms = np.zeros((3, 6), dtype=np.float32)
    with pytest.raises(ValueError):
        ll.dyeDensityAV1_arr(bad_atoms, source, 20.0, 2.0, 3.5, 0.9)

    bad_source = np.zeros((2,), dtype=np.float32)
    with pytest.raises(ValueError):
        ll.dyeDensityAV1_arr(atoms_nx4, bad_source, 20.0, 2.0, 3.5, 0.9)
