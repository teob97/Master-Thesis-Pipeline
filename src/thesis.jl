module thesis

using PyCall
using PRMaps
using Healpix
using Random, Distributions
import Stripeline as Sl

export Instrument
export get_single_map, get_foreground_maps, get_observations, get_white_noise
export add_white_noise!
export map2vec, fgbuster_basic_comp_sep

# Define python functions ---------------------------------------------------------------------------

# This function is automatically called when module is loaded
function __init__()
    py"""
    import pysm3
    import pysm3.units as u
    
    import pandas as pd
    import healpy as hp
    
    import fgbuster

    def pysm_sky_IQU(frequency, nside):
        sky = pysm3.Sky(nside=nside, preset_strings=["c1","d0","s0"])
        emission = sky.get_emission(frequency * u.GHz)
        emission = emission.to(u.uK_CMB, equivalencies=u.cmb_equivalencies(frequency * u.GHz))
        return pysm3.apply_smoothing_and_coord_transform(emission, rot=hp.Rotator(coord=("G", "C")))

    def get_df(instruments):
        buffer = []
        for i in instruments:
            buffer.append([i.frequency])       
        df = pd.DataFrame(buffer, columns = ['frequency'])        
        return df

    def fgbuster_pipeline(data, instruments):
        instrument = get_df(instruments)
        # Componenti settate per le mappe in polarizzazione [NON SICURO SE LA FREQ DI RIFERIMENTO SIA GIUSTA]
        components = [fgbuster.CMB(), fgbuster.Dust(353.), fgbuster.Synchrotron(23.)]
        return fgbuster.basic_comp_sep(components, instrument, data[:,1:])

    """
end

#------------------------------------------------------------------------------------------------------

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

# Genera rumero bianco per una singola PolarizedHealpixMap
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
    sigma = instrument.noisePerPixel * pixel_size * (180/π)^(2)
    noise = get_white_noise(signal.i.resolution.nside, sigma)

    signal.i = signal.i + noise.i
    signal.q = signal.q + noise.q
    signal.u = signal.u + noise.u

end

# Converte una lista di PolarizedHealpixMap in un vettore con row-major order,
# conpatibile con Python ed fgbuster
function map2vec(maps)
    v = zeros(maps[1].i.resolution.numOfPixels, 3, length(maps))
    for (index, m) in enumerate(maps)
        v[:,1,index] = m.i.pixels
        v[:,2,index] = m.q.pixels
        v[:,3,index] = m.u.pixels
    end
    v[isnan.(v)] .= -1.6375e+30 # TROVARE UN MODO MIGLIORE USANDO IN PYTHON healpy.UNSEEN
    return PyReverseDims(v)
end

function fgbuster_basic_comp_sep(maps, instruments)
    data = map2vec(maps)
    return py"fgbuster_pipeline"(data, instruments)    
end

end # module thesis

