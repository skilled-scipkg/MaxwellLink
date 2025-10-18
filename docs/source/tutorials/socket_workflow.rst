Socket TLS workflow
===================

.. important::

   This text walkthrough mirrors :doc:`notebook/socket_tls_workflow` and the
   regression test ``tests/test_tls/test_meep_2d_socket_tls1_relaxation.py``.
   Refer to the notebook for an interactive, executable version with captured
   outputs.

This tutorial mirrors ``tests/test_tls/test_meep_2d_socket_tls1_relaxation.py``
and ``tests/test_qutip/test_meep_2d_socket_qutip_tls_relaxation.py``. It shows
how to:

- Start a Meep simulation with a TCP ``SocketHub``.
- Launch a TLS or QuTiP driver that connects to the hub.
- Collect the returned population history and compare it to the analytical
  golden-rule decay.

1. Prepare the Meep script
--------------------------

Create ``run_tls.py`` with the following contents:

.. code-block:: python

   import numpy as np
   import meep as mp
   import maxwelllink as mxl
   from maxwelllink import sockets as mxs

   host, port = mxs.get_available_host_port()
   hub = mxl.SocketHub(host=host, port=port, timeout=10.0, latency=1e-5)
   print(f"Running TLS demo on {host}:{port}")

   molecule = mxl.Molecule(
       hub=hub,
       resolution=10,
       center=mp.Vector3(0, 0, 0),
       size=mp.Vector3(1, 1, 1),
       sigma=0.1,
       dimensions=2,
   )

   sim = mxl.MeepSimulation(
       hub=hub,
       molecules=[molecule],
       cell_size=mp.Vector3(8, 8, 0),
       boundary_layers=[mp.PML(3.0)],
       resolution=10,
       time_units_fs=0.1,
   )

   sim.run(until=90)

   population = np.array([entry["Pe"] for entry in molecule.additional_data_history])
   time_au = np.array([entry["time_au"] for entry in molecule.additional_data_history])
   print("Final excited population:", population[-1])

2. Launch the driver
--------------------

In a second terminal:

.. code-block:: bash

   mxl_driver.py --model tls --address 127.0.0.1 --port <port> \
     --param "omega=0.242, mu12=187, orientation=2, pe_initial=1e-4"

Replace ``<port>`` with the value printed by the Meep script. When launching the
driver from a different machine, substitute ``127.0.0.1`` with the public
hostname or IP address of the Meep node.

MaxwellLink pauses the Meep stepper until the driver connects. When the driver
exits, the EM simulation terminates cleanly. Replace ``--model tls`` with
``--model qutip --param "preset=custom, module=/path/to/build_tls.py"`` to use
the QuTiP implementation from the tests.

3. Compare with the analytical result
-------------------------------------

Append the following snippet to ``run_tls.py`` if you want to reproduce the test
assertions:

.. code-block:: python

   gamma = (0.1 ** 2) * (1.0 ** 2) / 2.0
   time_fs = time_au * 0.02418884254
   time_meep = time_fs / 0.1
   p_ref = population[0] * np.exp(-time_meep * gamma)
   rel_std = np.std(population - p_ref) / population[0]
   rel_max = np.max(np.abs(population - p_ref)) / population[0]
   print(f"std_dev={rel_std:.3e}, max_abs_diff={rel_max:.3e}")

The tolerance used in the automated test is ``3e-3`` for the standard deviation
and ``8e-3`` for the maximum deviation. Matching those numbers confirms that
your socket setup reproduces the reference behavior.
