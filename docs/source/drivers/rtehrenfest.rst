RT-TDDFT-Ehrenfest driver
=========================

The RT-TDDFT-Ehrenfest driver extends the RT-TDDFT model with classical nuclear
propagation using Ehrenfest or Born–Oppenheimer forces. It is implemented on top
of the same Psi4 backend and shares most parameters with
:doc:`rttddft <rttddft>`, adding controls for nuclear time stepping, stochastic
thermostats, and partial charges.

Requirements
------------

- ``psi4`` (same as the RT-TDDFT driver).
- Atomic masses, partial charges, or thermostat settings are optional and can be
  provided through parameters.

Example command
---------------

.. code-block:: bash

   mxl_driver.py --model rtehrenfest --port 31415 \
     --param "molecule_xyz=${PWD}/tests/data/hcn.xyz, functional=b3lyp, \
              basis=sto-3g, dt_rttddft_au=0.04, force_type=ehrenfest, \
              n_fock_per_nuc=1, n_elec_per_fock=10, \
              partial_charges=[1.0,-1.0,0.0], checkpoint=false, restart=false"

Additional parameters
---------------------

- ``force_type`` – ``ehrenfest`` for mean-field forces or ``bo`` to use
  Born–Oppenheimer gradients (see the ASE comparison tests).
- ``n_fock_per_nuc`` / ``n_elec_per_fock`` – Control the nested propagation loop
  (how many electronic steps per nuclear move).
- ``mass_amu`` – Override Psi4’s default atomic masses.
- ``friction_gamma_au`` / ``temperature_K`` / ``rng_seed`` – Enable Langevin
  thermostats.
- ``partial_charges`` – Optional per-atom charges for coupling to external
  fields.
- ``fix_nuclei_indices`` – Freeze selected nuclei.

Validation
----------

``tests/test_ase/test_ase_psi4_bomd.py`` uses the driver in two configurations
and compares the resulting trajectories to ASE/PSI4 BOMD runs:

- With ``force_type="bo"`` and a constant external field to match ASE’s BOMD
  driver.
- With ``force_type="ehrenfest"`` and zero external field, validating the use of
  Ehrenfest gradients.

These tests demonstrate how MaxwellLink’s Ehrenfest implementation reproduces
standard ASE dynamics, giving confidence when coupling to FDTD simulations.
