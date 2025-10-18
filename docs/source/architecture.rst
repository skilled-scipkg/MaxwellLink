Architecture
============

Architecture overview
---------------------

MaxwellLink separates electromagnetic propagation from molecular dynamics by
placing them in different processes (or even different nodes) and letting them
communicate through a socket protocol derived from the i-PI project:

1. The FDTD code (Meep at present) advances Maxwell's equations on its grid.
2. After each time step, MaxwellLink measures the regularised electric field at
   every coupled molecule and converts it to atomic units.
3. Those field vectors are sent to the driver processes through a ``SocketHub``
   barrier call.
4. Each driver propagates its molecular model for one FDTD step (possibly using
   sub-steps) and returns the time-derivative of the dipole moment and optional
   metadata.
5. The returned amplitudes are converted back to EM units and injected as Meep
   ``CustomSource`` objects before the next time step begins.

SocketHub
---------

``SocketHub`` (:mod:`maxwelllink.sockets`) manages the transport layer:

- Supports both TCP sockets (``host``/``port``) and UNIX domain sockets
  (``unixsocket``).
- Generates molecule IDs on demand or respects IDs provided by the FDTD script.
- Implements the ``NEEDINIT → INIT → READY/HAVEDATA`` handshake for each client.
- Detects dropped connections during sends or receives and pauses the EM solver
  until all expected drivers reconnect (see the reconnection loops exercised in
  ``tests/test_tls/test_meep_2d_socket_tls1_relaxation.py``).
- Exposes helpers such as :func:`maxwelllink.sockets.get_available_host_port`
  and :func:`maxwelllink.sockets.mpi_bcast_from_master` so the same configuration
  can be shared across MPI ranks.

Molecules
---------

Molecule wrapper
----------------

``maxwelllink.Molecule`` covers both socket-coupled and in-process models. Pass
``hub=SocketHub(...)`` to connect to an external driver or ``driver="tls"`` (and
``driver_kwargs``) to instantiate the model locally. Every molecule records
diagnostics in ``additional_data_history`` and maintains a hashable
``polarization_fingerprint`` so that Meep can reuse spatial kernels for
identical molecules (critical for the superradiance tests).

Each molecule stores the EM resolution (``resolution``), dimensions (1D/2D/3D),
Gaussian width (``sigma``), and centre/size of its coupling volume. These values
are used by :mod:`maxwelllink.em_solvers.meep` to build the correct spatial
envelope and interpolate fields.

EM solver integration
---------------------

``MeepSimulation`` derives from :class:`meep.Simulation` and automatically
inserts the appropriate coupling step when ``run`` is called:

- With a ``SocketHub`` it delegates to :func:`maxwelllink.em_solvers.meep.update_molecules`.
- Without a hub it falls back to
  :func:`maxwelllink.em_solvers.meep.update_molecules_no_socket`.
- Meep's Courant factor must currently be 0.5; the helper checks for this
  constraint on construction.

The module also defines a :class:`MeepUnits` adapter that converts between Meep
units and atomic units. Conversion factors were validated against the TLS
analytics checks in ``tests/test_tls/test_meep_2d_socket_tls1_relaxation.py`` and
``tests/test_tls/test_meep_2d_tlsmolecule_n_relaxation.py``.

When a full FDTD grid is unnecessary, :class:`maxwelllink.SingleModeSimulation`
approximates the field as a single damped harmonic oscillator evolving in atomic
units. It supports the same socket and non-socket molecule interfaces, making it
useful for rapid prototyping or unit tests without launching Meep.

MPI awareness
-------------

Meep simulations can be launched under MPI without extra code. Only the master
rank (rank 0) interacts with sockets; field integrals and returned amplitudes
are broadcast to the other ranks via ``mpi4py``. The helper methods
``am_master`` and ``mpi_bcast_from_master`` ensure that molecule IDs and hub
state remain in sync across ranks.

Resilience and checkpoints
--------------------------

Driver classes that inherit from :class:`maxwelllink.mxl_drivers.python.models.DummyModel`
support checkpointing. When ``checkpoint=true`` the driver writes state files
after each step; setting ``restart=true`` lets a reconnected driver resume from
disk. The hub blocks the EM solver inside ``wait_until_bound`` until all
expected molecules report back, so even long-lived RT-TDDFT simulations remain
consistent if a driver is restarted.
