Introduction
============

MaxwellLink provides a flexible software bridge between finite-difference
time-domain (FDTD) electromagnetic solvers and quantum or molecular dynamics
engines. The toolkit originated in the TEL Research Group at the University of
Delaware and is designed for self-consistent simulations in which light evolves
inside an FDTD grid while molecules respond through first-principles or model
drivers.

.. image:: ../../media/workflow.png
   :alt: MaxwellLink workflow diagram
   :align: center
   :scale: 25

Key capabilities
----------------

- Couple an external FDTD solver (currently Meep) to one or many molecular
  drivers running in the same process or across networked nodes.
- Exchange fields and dipole-response data through an i-PI inspired socket
  protocol that supports both TCP and UNIX domain sockets.
- Mix heterogeneous molecular theories within a single simulation, ranging from
  two-level systems and QuTiP model Hamiltonians to Psi4-based RT-TDDFT,
  Ehrenfest dynamics, and ASE-powered molecular mechanics.
- Tolerate driver pauses or restarts—the ``SocketHub`` automatically waits for
  reconnections before resuming a simulation step.
- Operate under MPI: only the master rank speaks to the drivers while field
  data and source amplitudes are broadcast to worker ranks.

Included driver families
------------------------

- **Model systems** – a lightweight TLS driver and a QuTiP interface for custom
  Hamiltonians with optional Lindblad terms.
- **First-principles quantum dynamics** – Psi4-backed RT-TDDFT and
  RT-TDDFT-Ehrenfest drivers capable of sub-stepping when the electronic time
  step is shorter than the FDTD step.
- **Classical and mixed quantum-classical dynamics** – ASE drivers that wrap any
  ASE calculator (including Psi4, ORCA, DFTB+) and a LAMMPS plugin using
  ``fix mxl``.

Learning path
-------------

1. :doc:`installation` walks through the recommended conda environment and
   optional driver dependencies.
2. :doc:`usage` shows how to launch a Meep simulation, connect a TLS driver, and
   inspect the returned dipole data.
3. :doc:`architecture` explains the socket protocol, molecular wrappers, and how
   MaxwellLink keeps the EM solver and drivers synchronized.
4. :doc:`drivers/index` documents all shipped molecular drivers, their CLI
   parameters, and expected outputs.
5. :doc:`tutorials/index` provides task-focused walkthroughs, from reproducing
   the TLS relaxation tests to running Psi4 RT-TDDFT with delta-kick spectra.
