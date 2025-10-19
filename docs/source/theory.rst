Theory
======

MaxwellLink provides a practical pathway to run self-consistent light-matter
simulations without forcing researchers to rewrite their preferred solvers.
The framework links a finite-difference time-domain (FDTD) electromagnetics
engine to one or many external molecular drivers through an intentionally thin
socket interface inspired by i-PI :cite:`XXX`. This page summarizes the central
ideas in plain language so that new users can reason about what the code is
doing under the hood before diving into model- or solver-specific details.

Overview of the Coupled Problem
-------------------------------

Self-consistent light-matter simulations propagate the electromagnetic (EM) field
and the microscopic matter response on the same time grid. In MaxwellLink:

- The EM field evolves on a Yee grid using Maxwell's equations with material
  response encoded by classical dielectric functions for coarse-grained media.
- Microscopic molecules are handled by external drivers that evaluate how their
  dipoles respond to the local EM field.
- A socket hub exchanges the regularized electric fields and induced dipole
  currents between the EM and molecular subsystems at every FDTD step.

As long as both sides agree on the time step and units, the entire simulation
remains synchronized and energy flows self-consistently between the field and
the molecules.

Classical Maxwell Solver
------------------------

The EM side solves the standard source-coupled Maxwell equations in
dimensionless FDTD units :cite:`XXX`:

.. math::

   \partial_t \mathbf{D} = \nabla \times \mathbf{H} - \mathbf{J}_\text{mol},
   \qquad
   \partial_t \mathbf{B} = -\nabla \times \mathbf{E}.

Material dispersion or absorption is encoded in the constitutive relations,
typically through permittivity models available in the chosen EM solver (Meep in
the initial release). The role of the molecular drivers is to deliver the
microscopic current density :math:`\mathbf{J}_\text{mol}` that sits on the right
hand side.

Regularized Electric Fields
---------------------------

Placing molecular dipoles on a discrete FDTD grid can lead to noisy, sometimes
unstable, self-fields if one simply samples the raw electric field at a single
grid point. MaxwellLink avoids that issue by using spatial kernel functions to
compute a smoothed, or *regularized*, electric field for each molecule:

.. math::

   \widetilde{\mathbf{E}}_m(t)
     = \int_{\Omega_m} \mathbf{E}(\mathbf{r}, t)\,
       \boldsymbol{\kappa}_m(\mathbf{r})\, \mathrm{d}\mathbf{r}.

The kernel :math:`\boldsymbol{\kappa}_m` is localized around the molecular
position and normalized so that the integral picks out the averaged field within
roughly one-tenth of the wavelength. Each driver receives only the three
components of :math:`\widetilde{\mathbf{E}}_m`, which keeps the data exchange
minimal while filtering the singular self-interaction that would otherwise
effectively double-count the molecule's emission :cite:`XXX`.

Molecular Response Models
-------------------------

MaxwellLink intentionally pushes all microscopic physics into independent
drivers. They advance their internal state, compute updated dipoles, and report
the time derivative of that dipole back to the EM solver. Three driver families
ship with the toolkit:

Model systems
   Lightweight quantum models—such as a two-level system or arbitrary
   QuTiP-defined Hamiltonians—propagate a von Neumann equation with optional
   Lindblad terms for dissipation. The driver evaluates
   :math:`\langle \hat{\boldsymbol\mu} \rangle` and its time derivative so that
   the EM solver can form :math:`\mathbf{J}_\text{mol}`.

Real-time electronic structure
   Psi4-based real-time time-dependent density functional theory (RT-TDDFT)
   integrates the electronic density matrix with enforced time-reversal symmetry
   propagators :cite:`XXX`. Sub-stepping can be applied internally when the
   electronic time step is smaller than the EM step.

Classical and ab initio molecular mechanics
   ASE and LAMMPS interfaces evaluate Born–Oppenheimer forces, update atomic
   positions, and compute the classical molecular dipole. The drivers return the
   dipole velocity :math:`\dot{\boldsymbol\mu}` required by Maxwell's equations.

Because every driver implements the same socket contract—read the regularized
field, advance itself, and send the dipole derivative—users can mix heterogeneous
levels of theory within a single simulation. For example, a plasmonic geometry
might simultaneously host a two-level system, a RT-TDDFT molecule, and a
classical water cluster.

Socket-Based Synchronization
----------------------------

The socket hub (``SocketHub``) maintains a registry of ``SocketMolecule`` objects,
each associated with one molecular driver process. At every EM time step:

1. Each ``SocketMolecule`` computes its regularized field and pushes it to the hub.
2. The hub forwards the field to the corresponding driver and awaits a response.
3. Drivers perform a *trial* propagation, send back their dipole derivatives, and
   acknowledge completion.
4. Once the hub has gathered all responses, it commits the step, hands the dipole
   derivatives to the EM solver, and instructs the drivers to finalize the new state.

This two-phase commit pattern keeps the coupled simulation robust against transient
driver failures or restarts. If a driver disconnects, the hub pauses the global
step, waits for reconnection, and replays the pending field data so that the
replacement driver rejoins without desynchronizing its peers :cite:`XXX`.

Units and Scaling
-----------------

EM solvers such as Meep operate in unit systems where
:math:`c = \epsilon_0 = \mu_0 = 1`, while quantum chemistry packages typically use
atomic units. MaxwellLink performs all conversions internally. Users only need to
specify the mapping between the FDTD time unit and real time (for example,
choosing a femtosecond scale), after which the framework consistently converts:

- EM time steps into atomic units for the drivers.
- Regularized electric fields into atomic units before they leave the EM solver.
- Dipole derivatives back into FDTD units for reinsertion into Maxwell's equations.

These conversions ensure you can combine disparate codes without manually tracking
unit systems or scaling factors.

Key Assumptions and Limitations
-------------------------------

The current design of MaxwellLink rests on a few practical assumptions:

- EM fields are treated classically; quantum EM effects (e.g., photon statistics)
  are not included :cite:`XXX`.
- Local matter is partitioned into coarse material response (handled by the EM
  solver) plus point-like microscopic molecules represented by drivers.
- Spatial kernel functions are static; phenomena requiring dynamic charge
  distributions or X-ray wavelengths may need extended kernels or finer grids.
- Exchange or charge-transfer interactions between explicit molecules and the
  classical background are not captured; scenarios such as plasmon-driven
  catalysis may require hybrid quantum-classical treatments beyond the current
  protocol.

Understanding these approximations helps determine whether MaxwellLink already
meets your needs or if additional method development is required for a given
application.

Where to Go Next
----------------

- The :doc:`architecture` page dives deeper into the socket implementation and
  runtime packaging.
- :doc:`drivers/index` documents all available drivers and their configuration.
- :doc:`tutorials/index` walks through complete examples, ranging from simple
  two-level systems to heterogeneous RT-TDDFT and molecular dynamics cases.
