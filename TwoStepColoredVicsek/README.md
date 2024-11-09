
# TwoStepColoredVicsek

This folder contains modified *c++* code for [hoomd-blue](https://github.com/glotzerlab/hoomd-blue) that introduces a new integration method for the equations of motion described [here](https://github.com/crpackard/Hoomd-ActiveTurningParticles/blob/main/README.md).

The locally saved [hoomd-blue](https://github.com/crpackard/Hoomd-ActiveTurningParticles/blob/main/hoomd-v2.8.2.tar.gz) package is automatically re-configured and built by the script [setup.sh](https://github.com/crpackard/Hoomd-ActiveTurningParticles/blob/main/setup.sh).

## File Descriptions

The *Python* file [integrate.py](https://github.com/crpackard/Hoomd-ActiveTurningParticles/blob/main/TwoStepColoredVicsek/integrate.py) has been modified to include the following additional integration method:

```python
class vicsek(_integration_method):

    def __init__(self, group, nlist=None, alpha=1.0, beta=1.0, v0=1.0, delta=0.0, bias=0.0, tau=1.0, sym=2, seed=0):

        # initialize base class
        _integration_method.__init__(self);

        ...
```

An example of how this class is called in order to run simulations is given in [run_simulation.py](https://github.com/crpackard/Hoomd-ActiveTurningParticles/blob/main/run_simulation.py).

When the ```vicsek(_integration_method)``` method is called, a call is internally made to the newly introduced *c++* module [TwoStepVicsek.cc](https://github.com/crpackard/Hoomd-ActiveTurningParticles/blob/main/TwoStepColoredVicsek/TwoStepVicsek.cc).
- This file contains the explicity implementation of the equations of motion.
- The name of the file derives from Vicsek-like manner in which the equations of motion are numerically solved; the orientations of particles are updated in step 1, and the particle positions are updated in step 2.

Finally, I note that the file [module-md.c](https://github.com/crpackard/Hoomd-ActiveTurningParticles/blob/main/TwoStepColoredVicsek/module-md.cc) has been modified so that ```cmake``` includes the *TwoStepVicsek* integration methods in the build.
