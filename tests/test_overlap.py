import numpy as np
import qutip
from qutip import ket

import krotov


def test_overlap():
    Qmagic = (1.0 / np.sqrt(2.0)) * qutip.Qobj(
        np.array(
            [[1, 0, 0, 1j], [0, 1j, 1, 0], [0, 1j, -1, 0], [1, 0, 0, -1j]],
            dtype=np.complex128,
        ),
        dims=[[2, 2], [2, 2]],
    )

    def ketbra(a, b):
        return ket(a) * ket(b).dag()

    rho_2 = ketbra('01', '10')

    expected = complex((Qmagic.dag() * rho_2).tr())
    res = krotov.second_order._overlap(Qmagic, rho_2)
    assert abs(res - expected) < 1e-14
