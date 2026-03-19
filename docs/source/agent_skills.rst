Agent Skills
============

Agent Skills for **MaxwellLink** provide a simple way to get started with this package. Using natural language prompts, users can easily create input files and run jobs on both local machines and HPC systems. Users can then inspect the input files and modify them as needed for more customized simulations.


Autonomous light-matter simulations with command line interfaces (CLIs)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A very elegant way to use MaxwellLink for autonomous light-matter simulations is to use the CLIs.

After installation, create a clean working directory and run ``mxl init`` within this directory to set up the **MaxwellLink** environment:

.. code-block:: bash

   mkdir myproject
   cd myproject
   mxl init

Then, the package knowledge base (including the agent skills) will be automatically loaded and the agent will be ready to assist with your simulations. You can then use any of the supported agent providers (e.g., ``codex``, ``claude``, or ``gemini`` CLIs, IDE extensions, or desktop applications) to provide natural language prompts for setting up and running MaxwellLink simulations.

For persistent HPC defaults across projects, install your profile once:

.. code-block:: bash

   mxl hpc set path/to/HPC_PROFILE.json

This writes ``~/.maxwelllink/HPC_PROFILE.json``. On systems where ``sbatch`` is available, ``mxl init`` links local ``HPC_PROFILE.json`` to that persistent profile.

For example, with the ``codex`` CLI, you can run the following command to start chatting with the agent: 

.. code-block:: bash

   codex 

Finally, after the simulations, you can remove the **MaxwellLink** runtime knowledge base by running:

.. code-block:: bash

   mxl clean


Prerequisites
~~~~~~~~~~~~~~

This feature requires cloning the MaxwellLink GitHub repository to access the skills folder. First, either on your local machine or on an HPC cluster, clone the repository and install **MaxwellLink**:

.. code-block:: bash

   git clone https://github.com/TaoELi/MaxwellLink.git
   cd MaxwellLink
   pip install .

Then, follow :doc:`installation` to install the third-party EM solver (MEEP FDTD) and molecular drivers.

After installation, you can use MaxwellLink's agent skills to run light-matter simulations automatically.


.. note::

   An agent provider, such as ``codex``, ``claude``, or ``gemini`` [command line interfaces (CLIs), IDE extensions, or desktop applications] needs to be installed and configured in your machine. 


Autonomous light-matter simulations on local machines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Please watch the following `walkthrough video <https://www.youtube.com/watch?v=ttAvvNByMLg>`_ for an introduction to using MaxwellLink's agent skills with VS Code and the Codex extension:

.. youtube:: ttAvvNByMLg

The aboe video tutorial contains the following steps:

- Open VS Code -> ``File`` -> ``Open Folder...`` -> select ``path/to/MaxwellLink``.
- Install/enable the **Codex** extension (from Marketplace). Make sure the extension has access to the workspace.
- Open the Codex chat panel (usually the side activity bar) and provide your prompt. The agent will load the relevant skills at ``skills/`` and try to accomplish your request.
- When prompted, let the agent run the suggested terminal commands in VS Code's integrated terminal; it will create ``projects/YYYY-MM-DD-<name>/`` and create input files for **MaxwellLink**.

The above video tutorial uses the following input prompt:

.. code-block:: text

   In my local machine, run an initially weakly excited two-level system coupled to 2d vacuum using meep fdtd and plot the excited-state population dynamics


Autonomous light-matter simulations on Anvil HPC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Using `Purdue Anvil HPC system <https://www.rcac.purdue.edu/knowledge/anvil>`_ (available via the `NSF ACCESS program <https://www.access-ci.org/>`_) as an example, we can also run MaxwellLink light-matter simulations.

Please watch the following `demo video <https://www.youtube.com/watch?v=c0PVxcvriDM>`_ for an introduction to using MaxwellLink's agent skills on HPC systems:

.. youtube:: c0PVxcvriDM

Different from the local machine setup, we need to **configure VS Code to connect to the remote HPC** system via SSH. Please refer to the official VS Code documentation for detailed instructions on setting up SSH connections:

`https://code.visualstudio.com/docs/remote/ssh <https://code.visualstudio.com/docs/remote/ssh>`_

After setting up the SSH connection, follow these steps:

- Open VS Code -> ``File`` -> ``Open Folder...`` -> select ``path/to/MaxwellLink`` on the remote HPC system.
- Install/enable the **Codex** extension (from Marketplace). Make sure the extension has access to the workspace.
- Open the Codex chat panel (usually the side activity bar), switch mode to **Agent (full access)** (bottom panel setting). This will allow to run commands without approval for each command.
- Provide your prompt. The agent will load the relevant skills at ``skills/`` and try to accomplish your request.
- When prompted, let the agent run the suggested terminal commands in VS Code's integrated terminal; it will create ``projects/YYYY-MM-DD-<name>/`` and create input files for **MaxwellLink**.

The above video tutorial uses the following input prompt:

.. code-block:: text

   In this Anvil HPC system, run a slurm job of an initially weakly excited two-level system coupled to 2d vacuum using meep fdtd and then plot the excited-state population dynamics after the slurm job


Customizing HPC settings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use the root-level ``HPC_PROFILE.json`` schema for cluster defaults, and install a persistent user profile with:

.. code-block:: bash

   mxl hpc set path/to/HPC_PROFILE.json

This updates ``~/.maxwelllink/HPC_PROFILE.json`` and is reused in future ``mxl init`` workspaces.



Vibe simulations
--------------------------------------

Users can run ``vibe simulations``, i.e., using natural language prompts to set up and run
**MaxwellLink** simulations, with popular agent providers.

.. note::

   This feature relies on installed ``codex``, ``claude``, or ``gemini`` command line interfaces (CLIs), IDE extensions, or desktop applications. 

Using ``codex`` for an example, install either `VS Code IDE + Codex extension <https://developers.openai.com/codex/ide>`_ or `Codex CLI <https://developers.openai.com/codex/cli/>`_.

Then, open Codex at the root directory of the **MaxwellLink** repository (requiring installing **MaxwellLink** from source) and start to chat:

.. code-block:: text

   In my local machine, run an initially weakly excited two-level system coupled to 2d vacuum using meep fdtd and plot the excited-state population dynamics

It supports ``vibe simulations`` on both local machines and HPC clusters. See :doc:`agent_skills` for more details.

.. note::

   When running ``vibe simulations`` on a local machine (**not on HPC**), the agent typically runs in a sandboxed environment. This may conflict with the MPI environment used by `Meep <https://meep.readthedocs.io/en/latest/>`_, leading to failed simulations.
   To resolve this issue, consider installing the serial version of `Meep <https://meep.readthedocs.io/en/latest/>`_ (i.e., without MPI support) instead:

   .. code-block:: bash

      conda install -c conda-forge pymeep="*=nompi_*"

Testing and Validation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, ``pytest`` does not run unit tests for ``vibe simulations`` that require agent providers. To run these tests, first install ``codex`` CLI (the other agent providers are not supported in the unit tests), 
then make sure you have logged in with your OpenAI API key or account credentials. 

Run the unit tests with the ``agent`` marker:

.. code-block:: bash

   # run agent tests, default run on local machines
   pytest -m agent -v
   # run agent tests on local machines
   pytest -m agent --codex-prompts local -v
   # run agent tests on hpc machines
   pytest -m agent --codex-prompts hpc -v
   # run agent tests by skipping the passed tests in the last run
   pytest -m agent -v --codex-resume

.. note::

   A **passed** simulation here means that the agent was able to generate a working simulation script and get some results. It does **not guarantee** that the simulation results are **physically accurate**.


Running a batch of prompts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Users can also use this unit test feature to automatically run a set of Codex prompts (each line corresponds to one independent prompt) stored in a text file:

.. code-block:: bash

   # run agent tests, default run on local machines
   pytest -m agent -v --codex-prompts-file path/to/your/prompts.txt

The output is stored at ``tests/test_agents/your_prompt_file_name/`` by default.
