import qutip as qt
import numpy as np


def build_model():
    # read from file
    energies = np.loadtxt("qutip_input/hcn_tddft_excitation_energies.txt")
    dipoles = np.loadtxt("qutip_input/hcn_tddft_transition_dipoles.txt")

    # basis
    n = len(energies) + 1  # number of states including ground state
    g = qt.basis(n, 0)
    H0 = 0.0 * g * g.dag()
    muz = 0.0 * g * g.dag()
    for i in range(1, n):
        e = qt.basis(n, i)
        omega = energies[i - 1]
        mu12 = dipoles[i - 1, 2]  # z direction only
        H0 += omega * e * e.dag()

        sigmax = g * e.dag() + e * g.dag()
        dip = mu12 * sigmax
        muz += dip

    mu_ops = {"x": None, "y": None, "z": muz}

    return dict(H0=H0, mu_ops=mu_ops)
