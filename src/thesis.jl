module thesis

using PyCall
using PRMaps
using Healpix
import Stripeline as Sl

export makeSky

# This function is automatically called when module is loaded
# I define here all the python function that I need
function __init__()
    py"""
    import pysm3
    import pysm3.units as u
    import healpy as hp

    def pysm_sky_IQU(nside, frequency):
        sky = pysm3.Sky(nside=nside, preset_strings=["c1","d0","s0"])
        emission = sky.get_emission(frequency * u.GHz)
        emission = emission.to(u.uK_CMB, equivalencies=u.cmb_equivalencies(frequency * u.GHz))
        return pysm3.apply_smoothing_and_coord_transform(emission, rot=hp.Rotator(coord=("G", "C")))
    """
end

Base.@kwdef struct Instrument
    name :: String = ""
    frequency :: Float64 = 0.0
    noisePerPixel :: Float64 = 0.0
end

# Produce HealpixPolarizedMap in μK_CMB
function get_sky(nside, frequency)
    
    map = Healpix.PolarizedHealpixMap{Float64, RingOrder}(nside)
    py_map = py"pysm_sky_IQU"(nside, frequency)

    map.i.pixels = py_map[1,:]
    map.q.pixels = py_map[2,:]
    map.u.pixels = py_map[3,:]

    return map
end

function get_emissions()
    # Passare un vettore di Instrument e per ognuno calcolare la mappa totale alla data frequenza.
    nothing
end


end # module thesis

