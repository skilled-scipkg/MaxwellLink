ASE driver
==========

The ASE driver embeds MaxwellLink in the `Atomic Simulation Environment
<https://wiki.fysik.dtu.dk/ase/>`_, letting you reuse any ASE-compatible
electronic structure calculator (Psi4, ORCA, DFTB+, …) while coupling to the EM
field through effective forces.

Requirements
------------

- ``ase`` (install via ``conda install -c conda-forge ase``).
- The desired calculator backends must be available on ``PATH`` or importable
  (e.g. Psi4 or ORCA Python bindings).

Example command
---------------

.. code-block:: bash

   mxl_driver.py --model ase --port 31415 \
     --param "atoms=${PWD}/tests/data/hcn.xyz, calculator=psi4, \
              calc_kwargs=method=b3lyp,basis=sto-3g, \
              charges=[1.0,-1.0,0.0], thermostat=MaxwellBoltzmann(300 K)"

Key parameters
--------------

- ``atoms`` – Path to an XYZ file or any format readable by ``ase.io.read``.
- ``calculator`` – Name of the ASE calculator to wrap (``psi4``, ``orca``,
  ``dftb`` …). MaxwellLink falls back to importing ``ase.calculators.<name>``.
- ``calc_kwargs`` – Comma-separated keyword arguments passed to the calculator
  constructor (strings, numbers, or lists—handled by the internal parser).
- ``charges`` – Optional per-atom charges. When omitted, the driver can
  recompute charges every step if the calculator exposes them.
- ``thermostat`` – Notation for ASE thermostat initialisation (e.g.
  ``MaxwellBoltzmann(300)``). Optional.

Operation
---------

- The driver integrates nuclear motion through ASE’s ``VelocityVerlet`` scheme.
- External electric fields are converted from atomic units to ASE’s eV/Å forces
  using the same conversion factors as the RT-TDDFT drivers.
- The ``additional_data_history`` contains serialized snapshots (positions,
  velocities, energies) returned by ``ASEModel._snapshot()``.

Validation
----------

``tests/test_ase/test_ase_psi4_bomd.py`` compares:

- MaxwellLink’s ``RTEhrenfestModel`` running with Psi4 forces.
- The ASE driver using Psi4 in Born–Oppenheimer mode.

The matching trajectories confirm that the ASE interface preserves the expected
forces and unit conversions.
