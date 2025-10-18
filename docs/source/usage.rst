Usage Guide
===========

This guide walks through three ways to couple :class:`maxwelllink.Molecule` with
Meep using a single TLS molecule as the working example.

.. note::

   MaxwellLink still ships the legacy :class:`maxwelllink.SocketMolecule` API
   used throughout ``tests/test_tls``. The patterns below focus on the unified
   :class:`maxwelllink.Molecule` interface, which works for both socket and
   embedded drivers; swap in ``SocketMolecule`` plus
   :func:`maxwelllink.update_molecules` if you need full parity with the tests.

Single process (no sockets)
---------------------------

The simplest setup instantiates the TLS driver inside the Meep process. This
avoids socket traffic entirely and is ideal for prototyping or small-scale
benchmarks.

.. code-block:: python

   import meep as mp
   import maxwelllink as mxl

   tls = mxl.Molecule(
       driver="tls",
       resolution=10,
       center=mp.Vector3(0, 0, 0),
       size=mp.Vector3(1, 1, 1),
       sigma=0.1,
       dimensions=2,
       driver_kwargs=dict(
           omega=0.242,      # TLS frequency (a.u.)
           mu12=187.0,       # dipole moment (a.u.)
           orientation=2,    # Ez
           pe_initial=1e-3,  # initial excited population
       ),
   )

   sim = mxl.MeepSimulation(
       molecules=[tls],
       cell_size=mp.Vector3(8, 8, 0),
       geometry=[],
       sources=[],
       boundary_layers=[mp.PML(3.0)],
       resolution=10,
       time_units_fs=0.1,
   )

   sim.run(until=90)

Within the same interpreter, you can analyze the TLS diagnostics through
``tls.additional_data_history``.

Local multi-process run (UNIX socket)
-------------------------------------

When you want Meep and the molecular driver to run as separate processes on the
same machine, use a UNIX domain socket. The hub listens on ``/tmp/socketmxl_<name>``
and the driver connects with ``--unix``.

Meep script:

.. code-block:: python

   import meep as mp
   import maxwelllink as mxl

   hub = mxl.SocketHub(unixsocket="tls_demo", timeout=10.0, latency=1e-4)

   tls = mxl.Molecule(
       hub=hub,
       resolution=10,
       center=mp.Vector3(0, 0, 0),
       size=mp.Vector3(1, 1, 1),
       sigma=0.1,
       dimensions=2,
   )

   sim = mxl.MeepSimulation(
       hub=hub,
       molecules=[tls],
       cell_size=mp.Vector3(8, 8, 0),
       boundary_layers=[mp.PML(3.0)],
       resolution=10,
       time_units_fs=0.1,
   )

   sim.run(until=90)

Driver command (same laptop/workstation):

.. code-block:: bash

   mxl_driver.py --model tls --unix --address tls_demo \
     --param "omega=0.242, mu12=187, orientation=2, pe_initial=1e-3"

UNIX sockets avoid port collisions and usually require no firewall changes.
MaxwellLink waits for the driver to connect before advancing the simulation.

Distributed run (TCP socket)
----------------------------

For multi-node deployments (e.g., Meep on one node and ``mxl_driver`` on another),
use a TCP socket. Let the OS pick a free port to prevent clashes.

Meep script:

.. code-block:: python

   import meep as mp
   import maxwelllink as mxl
   from maxwelllink import sockets as mxs

   _, port = mxs.get_available_host_port()
   hub = mxl.SocketHub(host="", port=port, timeout=30.0, latency=1e-4)

   print(f"SocketHub listening on port {port}. Share the host/IP with the driver.")

   tls = mxl.Molecule(
       hub=hub,
       resolution=10,
       center=mp.Vector3(0, 0, 0),
       size=mp.Vector3(1, 1, 1),
       sigma=0.1,
       dimensions=2,
   )

   sim = mxl.MeepSimulation(
       hub=hub,
       molecules=[tls],
       cell_size=mp.Vector3(8, 8, 0),
       boundary_layers=[mp.PML(3.0)],
       resolution=10,
       time_units_fs=0.1,
   )

   sim.run(until=90)

Setting ``host=\"\"`` binds the hub to all interfaces (equivalent to ``0.0.0.0``).
Share the public hostname or IP of the Meep node together with ``port``.

Driver command (run on the remote node):

.. code-block:: bash

   mxl_driver.py --model tls --address <meep-hostname> --port <port> \
     --param "omega=0.242, mu12=187, orientation=2, pe_initial=1e-3"

Replace ``<meep-hostname>`` with the reachable address of the Meep node. Open
firewall ports if required by your cluster configuration.

Inspecting TLS output
---------------------

All three workflows expose the TLS diagnostics through the molecule wrapper.
After ``sim.run`` completes, the code below recovers the excited-state
population and the simulation clock in atomic units.

.. code-block:: python

   import numpy as np

   population = np.array([entry["Pe"] for entry in tls.additional_data_history])
   time_au = np.array([entry["time_au"] for entry in tls.additional_data_history])

Compare the trajectory to the analytical golden-rule result as shown in
``tests/test_tls/test_meep_2d_socket_tls1_relaxation.py`` to validate a new build.

MPI execution
-------------

MaxwellLink detects MPI automatically. Only rank 0 communicates with drivers
while field integrals and returned amplitudes are broadcast to worker ranks via
``mpi4py``. Launch a run with:

.. code-block:: bash

   mpirun -np 4 python run_tls.py

When using sockets under MPI, ensure the driver process is launched exactly once
(typically by rank 0).

Driver restarts
---------------

If a driver disconnects unexpectedly, the hub pauses the Meep time loop and
waits for the driver to reconnect. Enabling ``checkpoint=true`` and
``restart=true`` in the driver parameters lets long RT-TDDFT or Ehrenfest jobs
recover from transient failures without restarting the EM simulation.

Single-mode cavity emulator
---------------------------

For quick prototyping without launching Meep, use
:class:`maxwelllink.SingleModeSimulation`. It models the EM field as one damped
harmonic oscillator in atomic units and couples to the same ``Molecule``
objects. Example:

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
   )

   sim.run(steps=500)
   print(sim.field_history[-1], tls.additional_data_history[-1]["Pe"])
