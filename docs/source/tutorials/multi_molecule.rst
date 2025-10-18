Multiple molecules without sockets
==================================

.. important::

   This page distils the workflow from
   :doc:`notebook/multi_molecule_superradiance` and the regression test
   ``tests/test_tls/test_meep_2d_tlsmolecule_n_relaxation.py``. Use the notebook
   when you want live plots or to tweak parameters interactively.

This tutorial is based on
``tests/test_tls/test_meep_2d_tlsmolecule_n_relaxation.py``. It shows how to
instantiate many TLS models directly inside the Meep process and how MaxwellLink
shares polarization fingerprints so that common spatial envelopes only need to
be built once.

1. Build the molecule list
--------------------------

.. code-block:: python

   import numpy as np
   import meep as mp
   import maxwelllink as mxl

   n_tls = 10
   molecules = []
   for _ in range(n_tls):
       molecules.append(
           mxl.Molecule(
               driver="tls",
               resolution=10,
               center=mp.Vector3(0, 0, 0),
               size=mp.Vector3(1, 1, 1),
               sigma=0.1,
               dimensions=2,
               driver_kwargs=dict(
                   omega=0.242,
                   mu12=187 / np.sqrt(n_tls),
                   orientation=2,
                   pe_initial=1e-2,
               ),
           )
       )

2. Run the simulation
----------------------

.. code-block:: python

   sim = mxl.MeepSimulation(
       molecules=molecules,
       cell_size=mp.Vector3(8, 8, 0),
       geometry=[],
       boundary_layers=[mp.PML(3.0)],
       resolution=10,
       time_units_fs=0.1,
   )
   sim.run(until=90)

   population = np.array(
       [np.real(entry["Pe"]) for entry in molecules[0].additional_data_history]
   )
   time_au = np.array(
       [np.real(entry["time_au"]) for entry in molecules[0].additional_data_history]
   )

3. Analyse collective decay
---------------------------

Because all molecules share the same position and dipole orientation, the decay
rate scales with the number of emitters (superradiance). The test compares the
trajectory to the analytical formula:

.. code-block:: python

   gamma_single = (0.1 ** 2) * (1.0 ** 2) / 2.0
   gamma_collective = gamma_single * n_tls
   time_meep = time_au * 0.02418884254 / 0.1
   p_ref = np.exp(-time_meep * gamma_collective) / (
       np.exp(-time_meep * gamma_collective) + (1.0 - population[0]) / population[0]
   )

   rel_std = np.std(population - p_ref) / population[0]
   rel_max = np.max(np.abs(population - p_ref)) / population[0]
   print(f"std_dev={rel_std:.3e}, max_abs_diff={rel_max:.3e}")

The same analysis applies when molecules are distributed in space—adjust the
``center`` and ``size`` parameters as needed. Running without sockets is ideal
for quick prototyping or when the EM solver and molecular model share the same
compute resource.
