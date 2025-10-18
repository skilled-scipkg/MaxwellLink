Documentation
======================================

MaxwellLink couples finite-difference time-domain (FDTD) electromagnetics with
quantum and molecular dynamics packages through a lightweight socket interface.

.. image:: ../../media/workflow.png
   :alt: MaxwellLink workflow diagram
   :align: center
   :scale: 25

Use these docs to install the toolkit, build your first coupled simulation, and
explore the available molecular drivers.


.. toctree::
   :maxdepth: 1
   :caption: Get Started

   introduction
   installation
   usage

.. toctree::
   :maxdepth: 1
   :caption: Learn More

   architecture
   em_solvers/index
   drivers/index

.. toctree::
   :maxdepth: 1
   :caption: Tutorials

   tutorials/index
   tutorials/notebook/socket_tls_workflow.ipynb
   tutorials/notebook/multi_molecule_superradiance.ipynb
   tutorials/notebook/rttddft_hcn.ipynb
   tutorials/notebook/ase_vs_ehrenfest.ipynb
   tutorials/notebook/single_mode_tls.ipynb

.. toctree::
   :maxdepth: 1
   :caption: Reference

   api/maxwelllink
