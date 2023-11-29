import numpy as np
from numpy.testing import assert_allclose

from openseestools.basic import linspaceCoords3d


def test_single_axis():
    expected = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    x, _, _ = linspaceCoords3d(0, 0, 0, 4, 0, 0, 4)
    _, y, _ = linspaceCoords3d(0, 0, 0, 0, 4, 0, 4)
    _, _, z = linspaceCoords3d(0, 0, 0, 0, 0, 4, 4)
    assert_allclose(x, expected)
    assert_allclose(y, expected)
    assert_allclose(z, expected)


def test_equal_axis():
    expected = np.array([0, 0.5, 1.0, 1.5, 2.0])
    x, y, z = linspaceCoords3d(0, 0, 0, 2, 2, 2, 4)
    assert_allclose(x, expected)
    assert_allclose(y, expected)
    assert_allclose(z, expected)
