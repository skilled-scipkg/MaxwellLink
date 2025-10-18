RT-TDDFT driver
===============

The RT-TDDFT driver couples MaxwellLink to the Psi4 quantum chemistry package
for real-time time-dependent density functional theory propagation. It supports
delta-kick excitations, automatic sub-stepping when ``dt_rttddft_au`` is smaller
than the FDTD step, and optional checkpoint/restart workflows.

Requirements
------------

- ``psi4`` must be available in the driver environment.
- The molecule is provided through an XYZ file whose second line contains the
  ``charge multiplicity`` descriptor (e.g. ``0 1``).

Example command
---------------

.. code-block:: bash

   mxl_driver.py --model rttddft --port 31415 \
     --param "molecule_xyz=${PWD}/tests/data/hcn.xyz, functional=SCF, \
              basis=sto-3g, delta_kick_au=1e-1, delta_kick_direction=xyz, \
              dt_rttddft_au=0.04, electron_propagation=pc, threshold_pc=1e-6, \
              checkpoint=false, restart=false, memory=2GB, num_threads=1"

Key parameters
--------------

- ``molecule_xyz`` – Path to the geometry file. Relative paths are resolved on
  the driver side.
- ``functional`` / ``basis`` – Any labels understood by Psi4 (``PBE``,
  ``B3LYP``, ``cc-pVDZ`` …).
- ``dt_rttddft_au`` – Electronic time step (a.u.). The driver automatically
  chooses an integer number of sub-steps per FDTD step.
- ``delta_kick_au`` / ``delta_kick_direction`` – Optional impulsive perturbation
  used for linear-response spectra.
- ``electron_propagation`` – ``etrs`` (enforced time-reversal symmetry) or ``pc``
  (predictor-corrector). ``threshold_pc`` controls the predictor-corrector
  convergence.
- ``remove_permanent_dipole`` – Subtract the permanent dipole from the coupling
  term.
- ``checkpoint`` / ``restart`` – Enable recovery when the driver restarts.

Returned diagnostics
--------------------

The driver pushes several diagnostics back to the EM process (see
``tests/test_psi4_rttddft/test_meep_2d_socket_rttddft.py``):

- ``time_au`` – Simulation clock in atomic units.
- ``mu_x_au``, ``mu_y_au``, ``mu_z_au`` – Time-dependent dipole moments.
- ``Etot_au`` – Total electronic energy (when ``verbose=true``).

These quantities can be saved to disk to generate spectra or to cross-check
against standalone Psi4 runs.
