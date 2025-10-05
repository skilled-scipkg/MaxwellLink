# Python drivers connecting with FDTD engines via MaxwellLink

A few python drivers for the coupled Maxwell-molecular dynamics simulations using the external FDTD engine are provided here.

**All parameters in python drivers are in atomic units**. MaxwellLink converts the FDTD units to atomic units when sending data to Python drivers. Similarly, Python drivers send data to MaxwellLink in atomic units, and MaxwellLink converts the data to FDTD units for propagating Maxwell's equations.

## **Electronic two-level system (tls)** 

A sample bash input for setting up the **tls** model is as follows: 

```bash
mxl_driver.py --model tls --port 31415 --param "omega=0.242, mu12=187, orientation=2, pe_initial=0.01" --verbose
```

1. The **--port** number should match that in **SocketHub**; see [../README.md](../README.md) for the initialization of **SocketHub** in MaxwellLink. 

2. The **--param** input provides the initialization of the Python class of the two-level system: [models/tls_model.py](./models/tls_model.py). 

3. **orientation=2** means setting the dipole vector to orient along the **z** (0-x, 1-y, 2-z) direction. 

4. **pe_initial=0.01** sets the initial electronic excited-state population.

5. If **checkpoint=true**, after each step, the necessary checkpoint data will be written to disk.

6. If **restart=true** and **checkpoint=true**, the driver code can restart from the checkpoint file in disk and resume the simulation. This setting is necessary when many  drivers are connected to the FDTD engine at the same time. If one driver is terminated in one machine (by different reasons), all the other drivers will pause the simulation and wait for the restart of this driver. 

7. If the driver and the FDTD code are running in different HPC nodes, please also set **--address <FDTD_CODE_IP_OR_DNS>** so that this driver can connect with FDTD.

## **An arbitrary model Quantum Hamiltonian supported by QuTiP (qutip)** 

The [QuTiP](https://qutip.org/) Python package can efficiently simulate a wide range of model Quantum systems with efficient support of Lindblian dissipation and quantum master solvers. 

We can propagate self-consistent Maxwell--Schr\"odinger equations for **any model Hamiltonian supported by QuTiP**.

A sample bash input for setting up the **qutip** model to simulate **a single TLS** is as follows: 

```bash
mxl_driver.py --model qutip --port 31415 --param "preset=tls, fd_dmudt=True, preset_kwargs=omega=0.242,mu12=187,orientation=2,pe_initial=1e-3,gamma_relax=0.0, gamma_dephase=0.0" 
```

1. **preset**: The available model Hamiltonians supported by this driver (**tls** or **custom**).

2. **fd_dmudt**: Whether using finite difference to evaluate dmu/dt, the quantity returned by the driver for constructing current densities in FDTD (default **false**). If set as **true**, analytical calculations will be used to evaluate dmu/dt, which can be expensive for large systems.

3. **preset_kwargs**: Additional kwargs (**omega=0.242,mu12=187,orientation=2,pe_initial=1e-3,gamma_relax=0.0, gamma_dephase=0.0**) to define the TLS as above. Here, **gamma_relax** defines the diagonal relaxation, and **gamma_dephase** defines the off-diagonal dephasing of the TLS. If these two values are provided, a Lindblan term will be constructed.


A sample bash input for setting up the **qutip** model to simulate **an arbitrary model quantum system** is as follows:
```bash
mxl_driver.py --model qutip --port 31415 --param "preset=custom, module=/path/to/spec.py, fd_dmudt=True, kwargs=omega=0.242,mu12=187,orientation=2,pe_initial=1e-3" 
```

1. **preset**: The available model Hamiltonians supported by this driver (**tls** or **custom**).

2. **module**: Path of the user-defined Python file. This file **must** define `build_model(**kwargs)` that returns a dict with keys:
```python
def build_model(**kwargs):
    return {
        "H0": qutip.Qobj (NxN),
        "mu_ops": {"x":Qobj|None, "y":Qobj|None, "z":Qobj|None},
        "c_ops": [Qobj, ...],
        "rho0":  Qobj (ket or density matrix),
}
```
Optional fields may be omitted; defaults: no c_ops, rho0=ground state if not provided.

For example, to define a TLS, we can provide a file **spec.py** as follows:

```python
# file spec.py
import numpy as np
import qutip as qt

def build_model(omega=0.242, mu12=187, orientation=2, pe_initial=0.0,
                gamma_relax=0.0, gamma_dephase=0.0):
    """
    Simple 2-level model preset like TLSModel, but using QuTiP objects.
    H0 = |e><e| * omega
    mu = mu12 * (|g><e| + |e><g|) along chosen axis
    Lindblad: relaxation (sigma_-) at rate gamma_relax, pure dephasing at gamma_dephase

    + **`omega`** (float): Transition frequency in atomic units (a.u.). Default is 0.242 a.u.
    + **`mu12`** (float): Dipole moment in atomic units (a.u.). Default is 187 a.u.
    + **`orientation`** (int): Orientation of the dipole moment, can be 0 (x), 1 (y), or 2 (z). Default is 2 (z).
    + **`pe`** (float): Initial population in the excited state. Default is 0.0.
    + **`gamma_relax`** (float): Relaxation rate (a.u.). Default is 0.0.
    + **`gamma_dephase`** (float): Pure dephasing rate (a.u.). Default is 0.0.
    """
    # basis
    g = qt.basis(2, 0)
    e = qt.basis(2, 1)
    H0 = omega * e * e.dag()

    # dipole operator along chosen axis; others set to None
    sigmax = g*e.dag() + e*g.dag()
    mux = muy = muz = None
    dip = mu12 * sigmax
    if orientation == 0:
        mux = dip
    elif orientation == 1:
        muy = dip
    else:
        muz = dip
    mu_ops = {"x": mux, "y": muy, "z": muz}

    # collapse operators
    c_ops = []
    if gamma_relax > 0.0:
        c_ops.append(np.sqrt(gamma_relax) * (g*e.dag())) 
    if gamma_dephase > 0.0:
        c_ops.append(np.sqrt(gamma_dephase) * (e*e.dag() - g*g.dag()))

    # initial state (density matrix form of a coherent state)
    pe = pe_initial
    rho0 = (1.0 - pe) * (g*g.dag()) + pe * (e*e.dag()) + np.sqrt(pe*(1.0-pe)) * (g*e.dag() + e*g.dag())

    return dict(H0=H0, mu_ops=mu_ops, c_ops=c_ops, rho0=rho0)
```

3. **kwargs**: Additional kwargs (**omega=0.242,mu12=187,orientation=2,pe_initial=1e-3**) to be parsed into the `build_model` function above.

4. **fd_dmudt**: Whether using finite difference to evaluate dmu/dt, the quantity returned by the driver for constructing current densities in FDTD (default **false**). If set as **true**, analytical calculations will be used to evaluate dmu/dt, which can be expensive for large systems.





## **Real-time time-dependent density functional theory (rttddft)** 

A sample bash input for setting up the **rttddft** model is as follows: 

```bash
mxl_driver.py --model rttddft --port 31415 --param "molecule_xyz=PATH_TO_MaxwellLink/tests/data/hcn.xyz, checkpoint=false, restart=false" --verbose
```

Similar as the **tls**, here the **--param** input provides the initialization of the Python class of the **rttddft** system: [models/rttddft_model.py](./models/rttddft_model.py). Currently, the Psi4 electronic structure engine is used to provide necessary quantities for propagating RT-TDDFT. 

We can freely change the molecular geometry, level of functional, basis, time step in RT-TDDFT, required memory, and number of CPU threads:
```bash
--param "molecule_xyz=PATH_TO_MaxwellLink/tests/data/hcn.xyz, functional="SCF",basis="sto-3g", dt_rttddft_au=0.04, memory="2GB",num_threads=1, checkpoint=false, restart=false"
```

A few tips when running RT-TDDFT with MaxwellLink:

1. The second line in the xyz file should be "charge multiplicity" (such as **0 1**). 

2. The time step in RT-TDDFT **dt_rttddft_au** (in atomic units) should be the same or smaller than the FDTD time step. If the RT-TDDFT time step is smaller than the FDTD time step, **dt_rttddft_au** will be adjusted automatically so an integer number of RT-TDDFT steps are propagated per FDTD step.

3. If **checkpoint** is **true**, after each FDTD step (or a few RT-TDDFT steps), the necessary checkpoint data in RT-TDDFT will be written to disk.

4. If **restart=true** and **checkpoint=true**, the driver code can restart from the checkpoint file in disk and resume the simulation. This setting is necessary when many  drivers are connected to MaxwellLink at the same time. If one driver is terminated in one machine (by different reasons), all the other drivers will pause the simulation and wait for the restart of this driver. 

5. If the driver and the FDTD code are running in different HPC nodes, please also set **--address <FDTD_CODE_IP_OR_DNS>** so that this driver can connect with FDTD.



## **Write your own driver** 

In [models/](./models/), each driver class inherits the **DummyModel** class. To implement a new Python driver, please add a new class file in [models/](./models/) and provide the implemenation of the class methods in **DummyModel** class.