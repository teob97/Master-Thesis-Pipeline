# Implementation and simulation of the pointing reconstruction model for the LSPE/Strip telescope.

This repository, together with [PRMaps](https://github.com/teob97/PRMaps.jl), contains the code used for the simulations in my [thesis](Thesis.pdf). Basically, it is a Julia port of some functions of the following Python libraries:
- [PySM](https://github.com/galsci/pysm): which allows to simulate a map of the sky to observe.
- [FGBuster](https://github.com/fgbuster/fgbuster): which allows you to perform a component separation analysis using frequencies associated with the following surveys: LSPE/Strip, LSPE/SWIPE, Planck, Quijote.

## System requirements

Once the Julia environment has been successfully activated, it is necessary to configure a Python environment containing PySM (v.3.3.2) and FGBuster (v2.0.0). It will then be necessary to set PyCall.jl correctly by indicating the path to the correct Python environment (see the documentation [here](https://github.com/JuliaPy/PyCall.jl#specifying-the-python-version)).

## Simple usage example

