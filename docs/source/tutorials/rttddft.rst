Psi4 RT-TDDFT tutorial
======================

.. important::

   This walkthrough complements :doc:`notebook/rttddft_hcn` and the regression
   test ``tests/test_psi4_rttddft/test_meep_2d_socket_rttddft.py``. Prefer the
   notebook for interactive execution and captured Psi4 output.

Follow this walkthrough to reproduce the vacuum HCN simulation from
``tests/test_psi4_rttddft/test_meep_2d_socket_rttddft.py``. It computes the
time-dependent dipole moment from Psi4 and compares it to a stored reference.

Prerequisites
-------------

- Install ``psi4`` in the same environment as MaxwellLink.
- Ensure the ``hcn.xyz`` sample geometry is available (``tests/data/hcn.xyz`` in
  the repository provides one with charge/multiplicity on the second line).

1. Start the Meep script
------------------------

``run_rttddft.py``:

.. code-block:: python

   import numpy as np
   import meep as mp
   import maxwelllink as mxl
   from maxwelllink import sockets as mxs

   host, port = mxs.get_available_host_port()
   hub = mxl.SocketHub(host=host, port=port, timeout=10.0, latency=1e-5)
   print(f"Listening on {host}:{port}")

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

   sim.run(until=5)

   mu_z = np.array([entry["mu_z_au"] for entry in molecule.additional_data_history])
   time_au = np.array([entry["time_au"] for entry in molecule.additional_data_history])
   np.savetxt("mu_z_vs_time.txt", np.c_[time_au, mu_z])

2. Launch the Psi4 driver
-------------------------

.. code-block:: bash

   mxl_driver.py --model rttddft --address 127.0.0.1 --port <port> \
     --param "molecule_xyz=${PWD}/tests/data/hcn.xyz, functional=SCF, \
              basis=sto-3g, delta_kick_au=1e-1, dt_rttddft_au=0.04, \
              electron_propagation=pc, threshold_pc=1e-6"

Again, substitute ``<port>`` with the value printed by the Meep script. Provide
the public hostname instead of ``127.0.0.1`` if the driver runs on a different
node.

The driver performs the initial SCF, optionally applies a delta-kick, and then
propagates for five FDTD steps. You can increase ``until`` in the script to
collect longer trajectories.

3. Validate against reference data
----------------------------------

The regression test bundles a reference file for the same setup. After running
the simulation, compare the output:

.. code-block:: bash

   python - <<'PY'
   import numpy as np
   ref = np.loadtxt("tests/test_psi4_rttddft/test_meep_2d_socket_rttddft_mu_z_au_ref.txt")
   calc = np.loadtxt("mu_z_vs_time.txt")
   assert np.allclose(calc[:, 0], ref[:, 0], atol=1e-8)
   assert np.allclose(calc[:, 1], ref[:, 1], atol=1e-8)
   print("RT-TDDFT trajectory matches the reference data.")
   PY

Adapting to different molecules is as simple as changing ``molecule_xyz`` and
the Psi4 options. When the FDTD time step is longer than the RT-TDDFT step, the
driver sub-steps automatically so that the EM solver always receives a single
response per Meep iteration.
