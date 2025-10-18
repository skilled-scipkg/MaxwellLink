Tutorials
=========

MaxwellLink's official tutorials are Jupyter notebooks rendered through
``nbsphinx``. Each notebook mirrors a regression test under ``tests/`` so you
can experiment interactively and then compare results against the automated
checks.

.. toctree::
   :maxdepth: 1

   notebook/socket_tls_workflow.ipynb
   notebook/multi_molecule_superradiance.ipynb
   notebook/rttddft_hcn.ipynb
   notebook/ase_vs_ehrenfest.ipynb
   notebook/single_mode_tls.ipynb

Notebook ↔ test mapping
-----------------------

- :doc:`notebook/socket_tls_workflow` → ``tests/test_tls/test_meep_2d_socket_tls1_relaxation.py``
- :doc:`notebook/multi_molecule_superradiance` → ``tests/test_tls/test_meep_2d_tlsmolecule_n_relaxation.py``
- :doc:`notebook/rttddft_hcn` → ``tests/test_psi4_rttddft/test_meep_2d_socket_rttddft.py``
- :doc:`notebook/ase_vs_ehrenfest` → ``tests/test_ase/test_ase_psi4_bomd.py``
- :doc:`notebook/single_mode_tls` → ``tests/test_tls/test_meep_2d_tlsmolecule_1_relaxation.py``

Text walkthroughs
-----------------

The repository still ships text-only guides derived from the same regression
tests. They focus on command sequences and code snippets and are useful when a
notebook runtime is unavailable.

.. seealso::

   ``tutorials/socket_workflow.rst`` · ``tutorials/multi_molecule.rst`` · ``tutorials/rttddft.rst`` · ``tutorials/ase_ehrenfest.rst``
