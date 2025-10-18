Single-Mode Cavity
==================

The :mod:`maxwelllink.em_solvers.single_mode_cavity` module implements a
lightweight cavity solver that replaces the full FDTD grid with a single damped
harmonic oscillator. It runs entirely in atomic units and can couple to the same
``Molecule`` abstractions used by the Meep backend, making it ideal for quick
prototyping, regression tests, or scenarios where a single optical mode
suffices.

Core API
--------

- :class:`maxwelllink.SingleModeSimulation` evolves the cavity field with the
  equation described in the class docstring. The constructor accepts:

  - ``dt_au`` – integration time step (must be positive).
  - ``frequency_au`` and ``damping_au`` – oscillator frequency :math:`\omega_0`
    and damping :math:`\kappa` in atomic units.
  - ``molecules`` – iterable of :class:`maxwelllink.Molecule` objects in either
    socket or embedded mode.
  - ``drive`` – constant float or callable returning the external drive term.
  - ``coupling_strength`` and ``coupling_axis`` – control how molecular dipoles
    feed back into the cavity. Axes accept ``0/1/2`` or ``"x"/"y"/"z"``.
  - ``record_history`` – when ``True`` stores time, field, velocity, drive, and
    molecular response traces on ``*.history`` attributes.

- :class:`maxwelllink.em_solvers.single_mode_cavity.MoleculeSingleModeWrapper`
  normalizes each molecule, refreshes its time units, and keeps an
  ``additional_data_history`` in sync with the cavity clock.

Socket integration
------------------

Socket-backed molecules must all share the same
:class:`maxwelllink.SocketHub`. Before each step the simulation calls
:meth:`maxwelllink.SocketHub.wait_until_bound` until every driver is connected.
Per-step payloads include the instantaneous electric field (aligned with the
selected axis) and accept JSON diagnostics via the driver's ``extra`` field.
See ``tests/test_tls/test_meep_2d_socket_tls1_relaxation.py`` for the same
socket handshake pattern used with Meep.

Embedded molecules
------------------

When ``mode="non-socket"``, MaxwellLink initializes the driver once via
:meth:`~maxwelllink.em_solvers.single_mode_cavity.MoleculeSingleModeWrapper.initialize_driver`.
Each step then calls :meth:`~maxwelllink.Molecule.propagate`,
collects the source amplitude with
:meth:`~maxwelllink.Molecule.calc_amp_vector`, and appends any diagnostics
exposed by the driver.

Example
-------

.. code-block:: python

   import maxwelllink as mxl

   tls = mxl.Molecule(
       driver="tls",
       driver_kwargs=dict(omega=0.242, mu12=187.0, orientation=2, pe_initial=1e-3),
   )

   sim = mxl.SingleModeSimulation(
       dt_au=0.05,
       frequency_au=0.242,
       damping_au=1e-3,
       molecules=[tls],
       drive=0.0,
       coupling_strength=1.0,
       record_history=True,
   )

   sim.run(steps=500)
   print(sim.field_history[-1], tls.additional_data_history[-1]["Pe"])

The example mirrors the snippet in ``docs/source/usage.rst`` and the notebook
``tutorials/notebook/single_mode_tls.ipynb``. History arrays let you compute
observables such as the excited-state population without re-running the driver.

Additional notes
----------------

- Both ``steps`` and ``until`` arguments are mutually exclusive when calling
  :meth:`~maxwelllink.SingleModeSimulation.run`.
- ``drive`` may be any callable taking the current time (a.u.) and returning a
  float. A constant value keeps the cavity coherently driven.
- Setting ``record_history=False`` skips all list allocations for minimal
  overhead when only the final field value is required.
