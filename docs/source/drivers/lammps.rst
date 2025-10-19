LAMMPS driver
=============

MaxwellLink ships a C++ ``fix mxl`` for LAMMPS so classical molecular dynamics
jobs can exchange electric fields and dipole currents with the EM solver. The
sources live in :mod:`maxwelllink.mxl_drivers.lammps` alongside a convenience
installer.

Requirements
------------

- POSIX-like environment (the fix uses BSD sockets).
- A LAMMPS build with atom IDs enabled and non-LJ units (such as ``units real`` or
  ``metal``).

Usage
-----

Socket preparation
^^^^^^^^^^^^^^^^^^

On the MaxwellLink side create a socket hub, for example in Meep:

.. code-block:: python

   hub = mxl.SocketHub(host="0.0.0.0", port=31415, timeout=60.0)

LAMMPS build helper
^^^^^^^^^^^^^^^^^^^

Use the bundled script to build a LAMMPS executable that contains the fix:

.. code-block:: bash

   mxl_install_lammps

The script clones LAMMPS, copies ``fix_maxwelllink.cpp``/``.h`` into ``src/MISC``,
and compiles an ``lmp_mxl`` binary placed on ``PATH``. Advanced users may copy
the fix into an existing source tree and rebuild manually.

LAMMPS input
^^^^^^^^^^^^

Add the fix to the group of atoms to be coupled:

.. code-block:: lammps

   fix 1 all mxl localhost 31415

During the run the fix connects to the socket hub, receives the electric field,
applies :math:`\mathbf{F}_i = q_i \mathbf{E}`, and returns
:math:`\mathrm{d}\boldsymbol{\mu}/\mathrm{d}t = \sum_i q_i \mathbf{v}_i`.
LAMMPS halts automatically if the hub closes the connection.


Parameters
----------

.. list-table::
   :header-rows: 1

   * - Name
     - Description
   * - ``host``
     - Hostname or IP address of the MaxwellLink process (third argument in the
       ``fix`` command). Required.
   * - ``port``
     - TCP port exposed by :class:`maxwelllink.SocketHub`. Must lie in ``(1024,
       65536]`` for TCP mode. Required.
   * - ``unix``
     - Optional flag (append ``unix`` to the ``fix`` command) to connect via a
       Unix-domain socket at ``/tmp/socketmxl_<host>``.

Returned data
-------------


- ``t_fs`` – Time in femtoseconds.
- ``temp_K`` – Temperature of the MD system in Kelvin.
- ``pe_au`` – Potential energy in atomic units (not implemented, always zero).
- ``ke_au`` – Kinetic energy in atomic units.
- ``dmudt_au`` – Dipole time derivative vector in atomic units.

Notes
-----

- Ensure every atom has a charge (``q``) defined; otherwise the external force is
  zero.
- The fix forces neighbor-list rebuilds every step, matching the requirement of
  exchanging fields at the FDTD cadence.
- See :mod:`maxwelllink.mxl_drivers.lammps.install` for scripting helpers that
  manage download and compilation.
