QuTiP driver
============

The QuTiP driver embeds MaxwellLink in the rich modelling ecosystem provided by
`QuTiP <https://qutip.org/>`_. It can run preset TLS models or load user-defined
Hamiltonians, collapse operators, and initial states from an external Python
module.

Requirements
------------

- ``qutip`` (install via ``conda install -c conda-forge qutip``).
- Any custom module referenced through ``module=...`` must be importable on the
  driver side.

TLS preset example
------------------

.. code-block:: bash

   mxl_driver.py --model qutip --port 31415 \
     --param "preset=tls, fd_dmudt=false, \
              preset_kwargs=omega=0.242,mu12=187,orientation=2,pe_initial=1e-3, \
              gamma_relax=0.0, gamma_dephase=0.0"

Key parameters
--------------

.. list-table::
   :header-rows: 1

   * - Name
     - Description
   * - ``preset``
     - ``tls`` for the built-in two-level system or ``custom`` to load a user
       module.
   * - ``fd_dmudt``
     - When ``true`` the driver evaluates :math:`\mathrm{d}\mu/\mathrm{d}t`
       through finite differences instead of closed-form matrix elements.
   * - ``preset_kwargs``
     - Comma-separated overrides for the TLS preset (same parameters as the TLS
       driver, plus Lindblad rates ``gamma_relax`` and ``gamma_dephase``).
   * - ``module``
     - (``preset=custom``) Path to a Python file that defines
       ``build_model(**kwargs)`` returning a dict with ``H0``, ``mu_ops``,
       optional ``c_ops``, and an initial density matrix ``rho0``.
   * - ``kwargs``
     - Parameters forwarded to ``build_model`` when using the custom preset.

Reference custom module
-----------------------

The TLS regression test
``tests/test_qutip/test_meep_2d_socket_qutip_tls_relaxation.py`` launches the
driver in ``preset=custom`` mode. The helper ``tests/test_qutip/build_tls.py``
defines:

.. code-block:: python

   def build_model(**kwargs):
       return {
           "H0": qutip.Qobj(...),            # bare Hamiltonian
           "mu_ops": {"x": mux, "y": muy, "z": muz},
           "c_ops": [ ... ],                 # optional collapse operators
           "rho0": qutip.Qobj(...),          # initial density matrix or ket
       }

Any unused dipole components can be set to ``None``. Collapse operators are
optional—omit ``c_ops`` for unitary dynamics.

Diagnostics
-----------

- Returned ``extra`` payloads match the TLS driver's ``Pe`` and ``time_au`` keys
  to ease comparison.
- The test compares the QuTiP trajectory against the analytical golden-rule
  decay and demonstrates how to embed user-defined Lindblad terms.
