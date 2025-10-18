ASE vs Ehrenfest cross-check
============================

.. important::

   This comparison mirrors :doc:`notebook/ase_vs_ehrenfest` and the regression
   test ``tests/test_ase/test_ase_psi4_bomd.py``. Move to the notebook for an
   interactive environment with inline diagnostics.

This tutorial reproduces ``tests/test_ase/test_ase_psi4_bomd.py`` to compare two
MaxwellLink drivers:

1. :class:`maxwelllink.RTEhrenfestModel` running with Psi4 forces.
2. :class:`maxwelllink.ASEModel` wrapping ASE’s Psi4 calculator.

Both drivers are executed in-process (no sockets) to highlight how MaxwellLink
can orchestrate purely molecular dynamics calculations before coupling them to
FDTD simulations.

Setup
-----

.. code-block:: python

   import numpy as np
   import maxwelllink as mxl

   xyz_path = "tests/data/hcn.xyz"

   rt = mxl.RTEhrenfestModel(
       engine="psi4",
       molecule_xyz=xyz_path,
       functional="b3lyp",
       basis="sto-3g",
       dt_rttddft_au=10.0,
       delta_kick_au=0.0,
       memory="2GB",
       verbose=False,
       remove_permanent_dipole=False,
       n_fock_per_nuc=1,
       n_elec_per_fock=1,
       homo_to_lumo=False,
       force_type="bo",
       partial_charges=[1.0, -1.0, 0.0],
   )
   rt.initialize(dt_new=10.0, molecule_id=0)

   ase_model = mxl.ASEModel(
       atoms=xyz_path,
       calculator="psi4",
       calc_kwargs="method=b3lyp, basis=sto-3g",
       charges="[1.0 -1.0 0.0]",
   )
   ase_model.initialize(dt_new=10.0, molecule_id=0)

Propagation
-----------

.. code-block:: python

   n_steps = 10
   field = np.array([1e-2, 1e-2, 0.0])
   rt_traj = []
   ase_traj = [ase_model._snapshot()["positions"] * 1.8897259886]

   for _ in range(n_steps):
       rt.propagate(effective_efield_vec=field)
       ase_model.propagate(effective_efield_vec=field)
       rt_traj.append(rt.traj_R[-1])
       snap = ase_model._snapshot()
       ase_traj.append(snap["positions"] * 1.8897259886)

   bonds_rt = [np.linalg.norm(R[0] - R[1]) for R in rt_traj]
   bonds_ase = [np.linalg.norm(R[0] - R[1]) for R in ase_traj]
   print("max bond deviation =", np.max(np.abs(np.array(bonds_rt) - np.array(bonds_ase))))

Interpretation
--------------

- With ``force_type="bo"`` both implementations agree exactly, confirming that
  the ASE wrapper applies the same electric-field forces as the Ehrenfest
  driver.
- Set ``force_type="ehrenfest"`` and ``charges="[0.0 0.0 0.0]"`` to switch the
  comparison to Ehrenfest forces, mirroring the second half of the regression
  test.

By validating molecular dynamics in isolation you can confidently couple the
same drivers to Meep and investigate field–matter interactions under more
realistic electromagnetic environments.
