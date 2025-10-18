TLS driver
==========

The TLS driver implements a closed two-level system with a single transition
frequency and dipole moment. It is lightweight, ships with MaxwellLink, and is
used extensively throughout the regression tests to validate the EM coupling.

Requirements
------------

- No third-party packages are required; the driver is installed with
  MaxwellLink.
- Optional ``scipy`` (already a transitive dependency) is used to build the
  short-time propagator.

Command-line usage
------------------

.. code-block:: bash

   mxl_driver.py --model tls --port 31415 \
     --param "omega=0.242, mu12=187, orientation=2, pe_initial=1e-4, \
              checkpoint=false, restart=false, verbose=false"

Parameters passed through ``--param``:

.. list-table::
   :header-rows: 1

   * - Name
     - Description
     - Default
   * - ``omega``
     - Transition frequency (a.u.), equals 1.0 in Meep units when
       ``time_units_fs=0.1``
     - ``0.242``
   * - ``mu12``
     - Dipole moment (a.u.); scales the emitted source amplitude
     - ``187``
   * - ``orientation``
     - Dipole orientation: ``0`` (``Ex``), ``1`` (``Ey``), ``2`` (``Ez``)
     - ``2``
   * - ``pe_initial``
     - Initial excited-state population
     - ``0.0``
   * - ``checkpoint``
     - Enable writing ``tls_checkpoint_id_<n>.npz`` after each step
     - ``false``
   * - ``restart``
     - Resume from an existing checkpoint instead of the initial state
     - ``false``
   * - ``verbose``
     - Print field values and density-matrix diagnostics each step
     - ``false``

Outputs
-------

- The driver reports the excited-state population ``Pe`` and current simulation
  time in the ``extra`` payload. When coupled to a molecule you can access the
  history through ``molecule.additional_data_history``.
- The analytical decay comparison used in
  ``tests/test_tls/test_meep_2d_socket_tls1_relaxation.py`` is a good validation
  scenario for new environments.

In-process usage
----------------

The same model can run without sockets by instantiating it via
``Molecule(driver="tls", driver_kwargs={...})``. This path is useful for unit
tests or when the EM solver and molecular model run in the same MPI task (see
``tests/test_tls/test_meep_2d_tlsmolecule_n_relaxation.py`` for an example).
