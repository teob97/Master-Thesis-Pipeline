module thesis

using PyCall
using PRMaps
using Healpix
using Random, Distributions
import Stripeline as Sl

export Instrument
export get_single_map, get_foreground_maps, get_observations, get_white_noise
export add_white_noise!

# This function is automatically called when module is loaded
# I define here all the python function that I need
function __init__()
    py"""
    import pysm3
    import pysm3.units as u
    import healpy as hp

    def pysm_sky_IQU(frequency, nside):
        sky = pysm3.Sky(nside=nside, preset_strings=["c1","d0","s0"])
        emission = sky.get_emission(frequency * u.GHz)
        #emission = emission.to(u.uK_CMB, equivalencies=u.cmb_equivalencies(frequency * u.GHz))
        return pysm3.apply_smoothing_and_coord_transform(emission, rot=hp.Rotator(coord=("G", "C")))
    """
end

Base.@kwdef struct Instrument
    name :: String = ""
    frequency :: Float64 = 0.0
    noisePerPixel :: Float64 = 0.0
end

# Produce PolarizedHealpixMap in μK_CMB
function get_single_map(frequency, nside)
    
    map = Healpix.PolarizedHealpixMap{Float64, RingOrder}(nside)
    py_map = py"pysm_sky_IQU"(frequency, nside)

    map.i.pixels = py_map[1,:]
    map.q.pixels = py_map[2,:]
    map.u.pixels = py_map[3,:]

    return map
end

# Calcola il cielo alle varie frequenze per ogni strumento
function get_foreground_maps(instruments::Array{Instrument}, nside)
    maps = Healpix.PolarizedHealpixMap[]
    for i in instruments
        push!(maps, get_single_map(i.frequency, nside))
    end
    return maps
end

# Lancia la scanning strategy per  i vari strumenti
# cam_ang è sempre lo stesso nell'ipotesi che i vari strumenti osservino gli stessi pixel nel cielo
function get_observations(
    cam_ang :: Sl.CameraAngles, 
    signals :: Vector{PolarizedHealpixMap}, 
    setup :: PRMaps.Setup
)
    observations = Healpix.PolarizedHealpixMap[]
    for signal in signals
        push!(observations, PRMaps.makeIdealMapIQU(cam_ang, signal, setup)[1])
    end
    return observations 
end

function get_white_noise(nside, sigma)

    map = PolarizedHealpixMap{Float64, RingOrder}(nside)

    noise_i = rand(Normal(0.0, sigma), nside*nside*12)
    noise_q = rand(Normal(0.0, sigma), nside*nside*12)
    noise_u = rand(Normal(0.0, sigma), nside*nside*12)

    map.i.pixels = noise_i
    map.q.pixels = noise_q
    map.u.pixels = noise_u
    
    return map
end

function add_white_noise!(signal, instrument)

    pixel_size = Healpix.nside2pixarea(signal.i.resolution.nside)
    sigma = instrument.noisePerPixel * pixel_size * (180/π)^2
    noise = get_white_noise(signal.i.resolution.nside, sigma)

    signal.i = signal.i + noise.i
    signal.q = signal.q + noise.q
    signal.u = signal.u + noise.u

end

end # module thesis

