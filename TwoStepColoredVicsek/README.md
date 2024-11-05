
# File Descriptions:

**HOOMD-Blue** simulations are run using *python* scripts which access modules in the compiled */build/hoomd* directory; these modules then access the *c++* code at the heart of **HOOMD-Blue** which run the simulations. 

In order to write a custom Vicsek model [integration method](https://hoomd-blue.readthedocs.io/en/stable/module-md-integrate.html), a *vicsek* class was added to */md/integrate.py* - this class handles the passing of simulation parameters from a python script to the *c++* code. 

Four *TwoStepVicsek* *c++* files were then written to define the Vicsek equations of motion for the integration method. Variable names and general coding conventions followed those used by the existing **HOOMD-Blue** integration methods.

In order to compile these modified files with the existing **HOOMD-Blue** build, two *cmake* files in */md* were modified to include the custom integration method.
