Installation
============

The recommended way to install MaxwellLink is through a conda environment that
provides Meep (with MPI support) alongside Python 3.9+.

Prerequisites
-------------

- A recent Python (``>=3.9``). MaxwellLink is tested with CPython.
- A working MPI stack (e.g. MPICH or OpenMPI) whenever you plan to run Meep or
  MaxwellLink under ``mpirun``.
- C/Fortran toolchains supplied by your platform if you plan to build drivers
  such as LAMMPS from source.

Create an environment
---------------------

.. code-block:: bash

   # Create a fresh environment and install Meep w/ MPI support
   CONDA_ENV="mxl"
   conda create -n "$CONDA_ENV" -c conda-forge pymeep="*=mpi_mpich_*" python=3.10
   conda activate "$CONDA_ENV"

Install MaxwellLink from source
-------------------------------

.. code-block:: bash

   git clone https://github.com/TaoELi/MaxwellLink.git
   cd MaxwellLink
   pip install .

This installs the Python package together with the ``mxl_driver`` console entry
point used to launch molecular drivers.

Optional driver dependencies
----------------------------

Install any drivers that you want to use from the command line. Each driver can
be pulled into the same conda environment.

.. code-block:: bash

   # Model Hamiltonians via QuTiP
   conda install -n "$CONDA_ENV" -c conda-forge qutip

   # Psi4 for RT-TDDFT and Ehrenfest drivers
   conda install -n "$CONDA_ENV" -c conda-forge psi4

   # Atomic Simulation Environment (ASE) for classical/MM drivers
   conda install -n "$CONDA_ENV" -c conda-forge ase

   # LAMMPS driver with fix mxl (installs a custom binary lmp_mxl)
   mxl_install_lammps

The TLS model ships with MaxwellLink and does not require extra packages. The
LAMMPS helper downloads, patches, and builds a LAMMPS executable that contains
``fix mxl``; alternatively copy the provided ``fix_maxwelllink.*`` sources into
your existing LAMMPS build and recompile.

Verify the installation
-----------------------

After installing the desired drivers, run the core regression tests to confirm
that the coupling between Meep and the drivers works in your environment. Test
selection relies on the markers defined in ``pyproject.toml`` (``core``,
``optional``, ``slow``), so prefer ``-m`` over name-based ``-k`` filtering.

.. code-block:: bash

   # Run the fast tests (TLS socket/no-socket paths)
   pytest -m core -v

   # Run the optional Psi4 and QuTiP tests when the dependencies are available
   pytest -m "optional or slow" -v

If the optional tests are skipped, confirm that the corresponding packages were
installed or adjust your ``PYTHONPATH`` so that the drivers can import them.
