ASE driver
==========

The ASE driver embeds MaxwellLink in the `Atomic Simulation Environment
<https://wiki.fysik.dtu.dk/ase/>`_, enabling Born–Oppenheimer molecular dynamics
with any ASE-compatible calculator. The implementation is provided by
:class:`maxwelllink.mxl_drivers.python.models.ASEModel`.

Requirements
------------

- ``ase`` (install via ``conda install -c conda-forge ase``).
- The desired calculator backends (Psi4, ORCA, DFTB+, …) must be installed and
  discoverable by ASE.

Usage
-----

Socket mode
^^^^^^^^^^^

.. code-block:: bash

   mxl_driver --model ase --port 31415 \
     --param "atoms=${PWD}/tests/data/hcn.xyz, calculator=psi4, \
              calc_kwargs=method=b3lyp,basis=sto-3g, \
              charges=[1.0 -1.0 0.0], n_substeps=5, temperature_K=300, \
              recompute_charges=false"

Non-socket mode
^^^^^^^^^^^^^^^

.. code-block:: python

   mxl.Molecule(
       driver="ase",
       driver_kwargs={
           "atoms": "tests/data/hcn.xyz",
           "calculator": "psi4",
           "calc_kwargs": "method=b3lyp,basis=sto-3g",
           "charges": "[1.0 -1.0 0.0]",
           "n_substeps": 5,
           "temperature_K": 300.0,
       },
       # ...
   )


Parameters
----------

.. list-table::
   :header-rows: 1

   * - Name
     - Description
   * - ``atoms``
     - ASE ``Atoms`` instance or path to a structure file readable by
       ``ase.io.read`` (``.xyz``, ``.pdb`` …). Required.
   * - ``calculator``
     - Name of the ASE calculator to wrap (``psi4``, ``orca``, ``dftb`` …).
       Default: ``psi4``.
   * - ``calc_kwargs``
     - Comma-separated ``key=value`` pairs (or a dict) forwarded to the
       calculator constructor. Default: ``""``.
   * - ``charges``
     - Optional per-atom charges specified as a space-separated list in square
       brackets, e.g. ``"[0.3 -0.3 0.0]"``. Default: ``None``.
   * - ``recompute_charges``
     - When ``True`` the driver queries the wrapped calculator for charges at
       every step instead of using ``charges``. Default: ``False``.
   * - ``n_substeps``
     - Number of velocity-Verlet steps per MaxwellLink time step. Default: ``1``.
   * - ``temperature_K``
     - Initial temperature passed to ``MaxwellBoltzmannDistribution`` before
       propagation. Default: ``0.0``.
   * - ``verbose``
     - When ``True`` print calculator setup and integration diagnostics.
       Default: ``False``.
   * - ``checkpoint``
     - When ``True`` write positions/velocities to ``ase_checkpoint_id_<n>.npz``.
       Default: ``False``.
   * - ``restart``
     - When ``True`` try to restore the most recent checkpoint on start-up.
       Default: ``False``.

Returned data
-------------

- ``time_au`` – Simulation time in atomic units.
- ``temperature_K`` – Instantaneous temperature reported by ASE.

Notes
-----

- Provide either ``charges`` or set ``recompute_charges=true``; the driver
  raises an error if no charges are available.
- Calculator-specific options can be supplied via ``calc_kwargs=...`` or as
  additional ``key=value`` pairs in ``--param``; unrecognised tokens are
  forwarded to the calculator constructor.
