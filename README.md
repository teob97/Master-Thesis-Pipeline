# Implementation and simulation of the pointing reconstruction model for the LSPE/Strip telescope.

This repository, together with [PRMaps](https://github.com/teob97/PRMaps.jl), contains the code used for the simulations in my [thesis](Thesis.pdf). Basically, it is a Julia port of some functions of the following Python libraries:
- [PySM](https://github.com/galsci/pysm): which allows to simulate a map of the sky to observe.
- [FGBuster](https://github.com/fgbuster/fgbuster): which allows you to perform a component separation analysis using frequencies associated with the following surveys: LSPE/Strip, LSPE/SWIPE, Planck, Quijote.

## System requirements

Once the Julia environment has been successfully activated, it is necessary to configure a Python environment containing PySM (v.3.3.2) and FGBuster (v2.0.0). It will then be necessary to set PyCall.jl correctly by indicating the path to the correct Python environment (see the documentation [here](https://github.com/JuliaPy/PyCall.jl#specifying-the-python-version)).

## Simple usage example

Below I report a simple example containing a component separation analysis implemented by introducing a systematic pointing error in the 43 GHz LSPE/Strip channel.

First it is necessary to choose which instruments you have to simulate; all of the aviable instruments are stored into `.pkl` file into the [instruments](instruments) directory:

```
using Pandas

swipe = read_pickle("../instruments/lspe_swipe_instrument.pkl")
strip = read_pickle("../instruments/lspe_strip_instrument.pkl")

instruments = concat([strip, swipe])
```

Then you have to decide the simulation setup for the LSPE/Strip telescope:

```
# Choose a Healpix map resolution
nside = 128

# Choose the sky model to observe
sky_model = "c1s0d0"

# Choose the starting day of the simulation and the number of days to simulate 
t_start = DateTime(2022, 1, 1, 12, 0, 0)
days = 10
    
setup = PRMaps.Setup(
    sampling_freq_Hz = 50.0,
    total_time_s = 24. * 3600. * days
)

# Choose the detector to simulate (in this case H0)
cam_ang = Sl.CameraAngles()
```