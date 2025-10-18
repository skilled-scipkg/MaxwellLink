LAMMPS driver
=============

MaxwellLink includes a helper to build a LAMMPS executable with ``fix mxl``,
allowing classical molecular dynamics simulations to exchange data with the EM
solver.

Installation
------------

.. code-block:: bash

   mxl_install_lammps

The command downloads LAMMPS, injects ``fix_maxwelllink.cpp`` / ``fix_maxwelllink.h``,
and compiles a binary named ``lmp_mxl`` that is added to your ``PATH``. Advanced
users can copy the fix into an existing LAMMPS source tree (``src/MISC``) and
rebuild manually instead.

LAMMPS input
------------

Add the MaxwellLink fix to the atom group you want to couple:

.. code-block:: lammps

   fix mxl all mxl localhost 31415

- ``localhost`` – Hostname or IP address of the machine running the FDTD code.
- ``31415`` – Port number supplied to the ``SocketHub`` in the EM script.

Workflow
--------

1. Start your MaxwellLink-enabled Meep script and initialize the hub with
   ``hub = mxl.SocketHub(host="0.0.0.0", port=31415, timeout=60.0)``.
2. Launch LAMMPS with the modified input file. The fix connects to the hub,
   completes the initialization handshake, and begins exchanging fields and
   dipole derivatives each time step.
3. When LAMMPS terminates the hub removes the client; the Meep simulation blocks
   until every connected driver finishes (mirroring the behavior for Python
   drivers).

See :mod:`maxwelllink.mxl_drivers.lammps` for API helpers that wrap the fix and
provide scripted install routines.
