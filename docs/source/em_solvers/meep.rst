Meep FDTD Solver
================

The :mod:`maxwelllink.em_solvers.meep` backend couples MaxwellLink molecules to
`Meep <https://meep.readthedocs.io/>`_ (``pymeep``) simulations. It converts
between Meep's native unit system and atomic units, injects molecular sources
into the grid, and streams regularized field integrals back to the quantum
drivers.

Overview
--------

- :class:`maxwelllink.MeepSimulation` subclasses :class:`meep.Simulation` and
  automatically inserts the MaxwellLink coupling step. It enforces
  ``Courant = 0.5`` and rescales each molecule to the chosen
  ``time_units_fs`` (default ``0.1`` fs).
- :class:`maxwelllink.em_solvers.meep.MoleculeMeepWrapper` builds Meep
  ``CustomSource`` objects for the molecule's polarization kernel and caches
  them so identical kernels share sources.
- :class:`maxwelllink.em_solvers.meep.MeepUnits` handles conversions such as
  :meth:`~maxwelllink.em_solvers.meep.MeepUnits.efield_em_to_au` and
  :meth:`~maxwelllink.em_solvers.meep.MeepUnits.source_amp_au_to_em`, and emits a
  units summary via :meth:`~maxwelllink.em_solvers.meep.MeepUnits.units_helper`.
- Two step functions expose the coupling logic directly:
  :func:`maxwelllink.update_molecules` (socket + MPI aware) and
  :func:`maxwelllink.update_molecules_no_socket` (embedded drivers only).

Prerequisites
-------------

- ``pymeep`` (imported as ``meep``) must be installed; import failures raise a
  clear :class:`ImportError`.
- Optional ``mpi4py`` support is detected at runtime. When present, only Meep's
  master rank talks to the molecular drivers while amplitudes and field
  integrals are broadcast across ranks.
- Socket workflows require a :class:`maxwelllink.SocketHub` and external driver
  processes such as ``mxl_driver.py``. End-to-end examples live in
  ``tests/test_tls/test_meep_2d_socket_tls1_relaxation.py`` and related files.

Typical socket workflow
-----------------------

.. code-block:: python

   import meep as mp
   import maxwelllink as mxl
   from maxwelllink import sockets as mxs

   host, port = mxs.get_available_host_port()
   hub = mxl.SocketHub(host=host, port=port, timeout=10.0, latency=1e-5)
   print(f"SocketHub listening on {host}:{port}")

   molecule = mxl.Molecule(
       hub=hub,
       resolution=10,
       center=mp.Vector3(0, 0, 0),
       size=mp.Vector3(1, 1, 1),
       sigma=0.1,
       dimensions=2,
   )

   sim = mxl.MeepSimulation(
       cell_size=mp.Vector3(8, 8, 0),
       geometry=[],
       sources=[],
       boundary_layers=[mp.PML(3.0)],
       resolution=10,
       hub=hub,
       molecules=[molecule],
       time_units_fs=0.1,
   )

   # launch the molecular driver separately (see tests for subprocess helpers)
   sim.run(until=90)

At runtime the coupling step sends regularized field integrals (Meep units) to
the driver, receives source amplitudes in atomic units, converts them back to
Meep units, and accumulates them inside the shared ``CustomSource`` objects. The
driver's diagnostic payloads (``extra`` JSON blobs) are appended to
``molecule.additional_data_history``.

.. note::
   Legacy scripts (see ``tests/test_tls/test_meep_2d_socket_tls1_relaxation.py``)
   still use :class:`maxwelllink.SocketMolecule` in conjunction with
   :func:`maxwelllink.update_molecules`. Both interfaces remain supported; the
   :class:`maxwelllink.Molecule` + :class:`maxwelllink.MeepSimulation` path shown
   above removes the need to wire the coupling loop manually.

Embedded / non-socket molecules
-------------------------------

When molecules are embedded (``mode="non-socket"``), wrap them with
:class:`~maxwelllink.em_solvers.meep.MoleculeMeepWrapper` and pass them to
:func:`~maxwelllink.update_molecules_no_socket`. MaxwellLink initializes each
driver once via :meth:`~maxwelllink.em_solvers.meep.MoleculeMeepWrapper.initialize_driver`
and reuses the cached polarization sources on subsequent runs.

MPI considerations
------------------

:func:`~maxwelllink.update_molecules` detects ``mpi4py`` automatically. It keeps
the socket exchange on rank 0 (and Meep's master), broadcasts amplitudes to all
ranks, and pauses cleanly if a driver disconnects until
:meth:`maxwelllink.SocketHub.wait_until_bound` confirms reconnection. This logic
matches the regression tests under ``tests/test_tls/test_meep_*``.

Further reading
---------------

- :doc:`tutorials/index` summarises the available notebooks and how they map to
  automated regression tests.
- :doc:`tutorials/notebook/socket_tls_workflow` demonstrates the end-to-end TLS
  socket workflow with executable cells.
- :doc:`tutorials/notebook/rttddft_hcn` couples the same Meep setup to the
  Psi4-based RT-TDDFT driver.
