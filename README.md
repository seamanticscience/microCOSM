# microCOSM
_microCOSM_ (Carbon/Climate Ocean Science Model) is a lightweight _in silico_ Planet that can be run from the command line or alternatively, and wonderfully conveniently, from a python or Jupyter Notebook environment. 

## Command line:
First, you must compile the exec for your particular destination.
```
>>make model OPTIONDEFS=""
```

There are currently two critical compile-time options to set through the OPTIONDEFS argument to _make_ or directly in the _Makefile_:
1. `WRITEOUTFILE`, which produces a text file for each unique model run. See `microCOSM000001.dat` for an example.
2. `FIXATMPCO2`, which controls whether the atmospheric box acts as an infinite carbon reservoir, or whether the atmosphere and ocean are coupled, and air-sea CO_2 fluxes can change the atmospheric pCO_2 in response to ocean carbon content changes.

Simply add `-DWRITEOUTFILE` and/or `-DFIXATMPCO2` to the `OPTIONDEFS` argument to _make_ or in in the _Makefile_.

For the command line version make sure to use `WRITEOUTFILE` to capture your output. You can see in _microcosm_model.F90_ how to assign the input values, with some examples from a "cold start", and some values from a 100kyr spinup simulation...this run took a minute or two to complete, so feel free to tinker! Execute as:
```
>>./microCOSM
```
## Python:

Python fans can use the _run_microCOSM.ipynb_ notebook to find out how to compile, access, and run the model from there. 

If you are having difficulty (i.e. non-zero exit status reported in the notebook), you can compile the python model from the command line
```
>>make pymodel OPTIONDEFS=""
```
which will provide lots of output from f2py for troubleshooting.

Advanced users may be interested in running _microCOSM_ with the [_pandarallel_ package](https://github.com/nalepae/pandarallel) "parallel_apply" function to do efficient ensemble/parameter space explorations - see [this notebook](https://github.com/seamanticscience/Lauderdale_etal_2020_PNAS/blob/master/boxmodel.ipynb) (and the _utils.py_ file) for the strategy, or ask me how!

## More info:

The carbon system in _microCOSM_ is solved using routiines from the _SolveSAPHE_ package (v1.0.1) of [Munhoven (2013)](https://doi.org/10.5194/gmd-6-1367-2013), particularly _mod_precision_, _mod_chemconst_, _mod_chemspeciation_, and _mod_phsolvers_. The only modification I made was to promote all `PRIVATE` functions and variables to public in order to work with `f2py` (sorry!). See those files for their license arrangements.

A progenitor of this model was originally developed by [Mick Follows](http://mick.mit.edu/about/). 

Any questions, comments or feature requests, please get in contact! I hope to develop some more interesting avenues (dual numbers, northern hemisphere land box, atmospheric dynamics, climate feedbacks, etc) and maybe come up with a family of versions based around other simple geometries (e.g. the old Toggweiler configurations).

### Troubleshooting tips:
1. Because f2py is finicky about allocatable output arrays, it is ideal to aim for 1000 output timesteps (i.e. 10kyrs at 10yr output, 100kyrs at 100yr output, etc.) You can, of course change this manually by altering `outstepmax` in _comdeck.h_ and _microCOSM_model.F90_.
2. Importing the model package was causing python environment crashes due to memory address segfaults. Turns out `f2py` was a bit confused about what version of OSX I was using and was cross compiling to Mavericks instead of Catalina. I was able to solve this by issuing:
```
>>export MACOSX_DEPLOYMENT_TARGET=10.15 
```
before compiling.
