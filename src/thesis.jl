module thesis

using PyCall
using PRMaps
using Healpix
using StatsPlots
using StructArrays
using Random, Distributions
using Logging
import Stripeline as Sl

export Instrument, run_simulation, run_simulation_with_error, get_corrplot, get_map_and_hist
export get_single_map, get_observations, get_white_noise
export add_white_noise!
export map2vec, fgbuster_basic_comp_sep

export get_foreground_maps, get_noise_maps

# Define python functions ---------------------------------------------------------------------------

# This function is automatically called when module is loaded
function __init__()
    py"""
    import pysm3
    import pysm3.units as u
    import fgbuster

    import pandas as pd
    import healpy as hp
    import numpy as np

    #def pysm_sky_IQU(frequency, nside):
    #    sky = pysm3.Sky(nside=nside, preset_strings=["c1","d0","s0"], output_unit=u.uK_RJ)
    #    emission = sky.get_emission(frequency * u.GHz)
    #    return pysm3.apply_smoothing_and_coord_transform(emission, rot=hp.Rotator(coord=("G", "C")))

    def get_freq_maps(instruments, sky_model, nside, unit):
        emissions = fgbuster.get_observation(instruments, sky_model, noise = False, nside = nside, unit = unit)
        for i,emission in enumerate(emissions):
            emissions[i] = pysm3.apply_smoothing_and_coord_transform(emission, rot=hp.Rotator(coord=("G", "C")))
        return emissions

    def get_noise_maps(instruments, nside, unit):
        return fgbuster.get_noise_realization(nside, instruments, unit)

    #def get_df(instruments):
    #    buffer = []
    #    for i in instruments:
    #        buffer.append([i.frequency])       
    #    df = pd.DataFrame(buffer, columns = ['frequency'])        
    #    return df

    def fgbuster_pipeline(instruments, data):
        components = [fgbuster.CMB(), fgbuster.Dust(353.), fgbuster.Synchrotron(23.)] 
        return fgbuster.basic_comp_sep(components, instruments, data[:,1:])

    def get_mvDistibution(result):
        return np.random.multivariate_normal(result['x'], result['Sigma'], 10000)
    
    """
end

#------------------------------------------------------------------------------------------------------

#= Base.@kwdef struct Instrument
    name :: String = ""
    frequency :: Float64 = 0.0
    noisePerPixel :: Float64 = 0.0
end =#

#--------------------------------------------------------------------------------------------------------

# Converte una lista di PolarizedHealpixMap in un vettore con row-major order,
# conpatibile con Python ed fgbuster
function map2vec(maps)
    v = zeros(maps[1].i.resolution.numOfPixels, 3, size(maps,1))
    for (index, m) in enumerate(maps)
        v[:,1,index] = m.i.pixels
        v[:,2,index] = m.q.pixels
        v[:,3,index] = m.u.pixels
    end
    v[isnan.(v)] .= -1.6375e+30 # TROVARE UN MODO MIGLIORE USANDO IN PYTHON healpy.UNSEEN
    return PyReverseDims(v)
end

function vec2map(vec)
    maps = Healpix.PolarizedHealpixMap[]
    for indx in axes(vec,1)
        map = PolarizedHealpixMap{Float64, RingOrder}(64)
        map.i.pixels = vec[indx,1,:]
        map.q.pixels = vec[indx,2,:]
        map.u.pixels = vec[indx,3,:]
        push!(maps, map)
    end
    return maps
end

#--------------------------------------------------------------------------------------------------------

# Calcola il cielo alle varie frequenze per ogni strumento e restituisce un PolarizedHealpixMap
function get_foreground_maps(
    instruments, 
    sky_model::String, 
    nside::Int; 
    unit::String = "uK_CMB"
)
    results = py"get_freq_maps"(instruments, sky_model, nside, unit)
    maps = vec2map(results)
    return maps
end

function get_noise_maps(
    instruments,
    nside::Int;
    unit::String = "uK_CMB" 
)
    results = py"get_noise_maps"(instruments, nside, unit)
    maps = vec2map(results)
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

function get_observations_with_error(
    cam_ang :: Sl.CameraAngles,
    tel_ang :: Sl.TelescopeAngles,
    signals :: Vector{PolarizedHealpixMap}, 
    setup :: PRMaps.Setup
)
    observations = Healpix.PolarizedHealpixMap[]
    for signal in signals
        push!(observations, PRMaps.makeErroredMapIQU(cam_ang, tel_ang, signal, setup)[1])
    end
    return observations 
end

# Genera rumero bianco per una singola PolarizedHealpixMap
#= function get_white_noise(nside, sigma)

    map = PolarizedHealpixMap{Float64, RingOrder}(nside)

    noise_i = rand(Normal(0.0, sigma), nside*nside*12)
    noise_q = rand(Normal(0.0, sigma), nside*nside*12)
    noise_u = rand(Normal(0.0, sigma), nside*nside*12)

    map.i.pixels = noise_i
    map.q.pixels = noise_q
    map.u.pixels = noise_u
    
    return map
end

function add_white_noise!(signal, instrument, setup)

    nside = signal.i.resolution.nside

    # Calcolo grandezza pixel in deg^2
    pixel_size = Healpix.nside2pixarea(nside) * (180.0/π)^(2)
    # Calcolo numero medio di volte che un pixel viene colpito
    average_hits = (setup.total_time_s * setup.sampling_freq_Hz) / (nside * nside * 12)
    # Divido errore pr la radice del numero medio di volte che pixel è colpito
    sigma = (instrument.noisePerPixel * pixel_size) / sqrt(average_hits)
    
    noise = get_white_noise(nside, sigma)

    signal.i = signal.i + noise.i
    signal.q = signal.q + noise.q
    signal.u = signal.u + noise.u

end =#



function fgbuster_basic_comp_sep(instruments, maps)
    data = map2vec(maps)
    return py"fgbuster_pipeline"(instruments, data)    
end


# Simulation without poining errors
function run_simulation(
    instruments,
    cam_ang :: Sl.CameraAngles,
    setup :: PRMaps.Setup,
    sky_model :: String,
    nside :: Int,
)
    @info "Generating sky signal and noise"
    signals = get_foreground_maps(instruments, sky_model, nside)
    noise = get_noise_maps(instruments, nside)

    @info "Simulating telescope scanning strategy [LSPE/Strip]"
    observations = get_observations(cam_ang, signals, setup)

    @info "Adding white noise based on the instruments sensitivity"
    for indx in axes(signals, 1)
        observations[indx].i.pixels += noise[indx].i.pixels
        observations[indx].q.pixels += noise[indx].q.pixels 
        observations[indx].u.pixels += noise[indx].u.pixels 
    end

    @info "Run component separation using fgbuster"
    return fgbuster_basic_comp_sep(instruments, observations)

end

function run_simulation_with_error(
    instruments :: StructArrays.StructArray,
    cam_ang :: Sl.CameraAngles,
    tel_ang :: Sl.TelescopeAngles,
    setup :: PRMaps.Setup,
    nside :: Int,
)
    # Separate LSPE/Strip experiment from the others
    if !any(instruments.name .== "LSPE/Strip")
        error("Invalid simulation: no Instrument named LSPE/Strip found")
    end
    # Remove LSPE/Strip from the instruments vector and store them in another vector 
    strip = filter(z -> z.name == "LSPE/Strip", instruments)
    filter!(z -> z.name != "LSPE/Strip", instruments)
    
    @info "Generating sky signal using pysm3"
    signals = get_foreground_maps(instruments, nside)
    signals_strip = get_foreground_maps(strip, nside)

    @info "Simulating telescope scanning strategy [LSPE/Strip WITH pointing error]"
    observations = get_observations(cam_ang, signals, setup)
    observations_strip = get_observations_with_error(cam_ang, tel_ang, signals_strip, setup)

    # Append LSPE/Strip result to the other results 
    StructArrays.append!!(instruments, strip)
    append!(observations, observations_strip)

    @info "Adding white noise based on the instruments sensitivity"
    for indx in size(signals, 1)
        add_white_noise!(observations[indx], instruments[indx], setup) 
    end

    @info "Run component separation using fgbuster"
    return fgbuster_basic_comp_sep(observations, instruments)  

end

function get_corrplot(result)
    sampling = py"get_mvDistibution"(result)
    return StatsPlots.corrplot(sampling, label = result["params"], size = (1500,700), fillcolor =:thermal, bottom_margin = 6Plots.mm, left_margin = 6Plots.mm)
end

function get_map_and_hist(result, stokes_param::String, nside::Int)

    if size(result["s"], 2) == 1
        if stokes_param == "I"
            stokes_indx = 1
        else
            error("Wrong stokes_param: it must be I")
        end
    elseif size(result["s"], 2) == 2
        if stokes_param == "Q"
            stokes_indx = 1
        elseif stokes_param == "U"
            stokes_indx = 2
        else
            error("Wrong stokes_param: it must be Q or U")
        end
    elseif size(result["s"], 2) == 3
        if stokes_param == "I"
            stokes_indx = 1
        elseif stokes_param == "Q"
            stokes_indx = 2
        elseif stokes_param == "U"
            stokes_indx = 3
        else
            error("Wrong stokes_param: it must be one between I,Q,U")
        end
    else
        error("Wrong result size")
    end

    result["s"][result["s"] .== -1.6375e+30] .= NaN

    cmb = HealpixMap{Float64, RingOrder}(nside)
    dust = HealpixMap{Float64, RingOrder}(nside)
    synchrotron = HealpixMap{Float64, RingOrder}(nside)

    cmb.pixels = result["s"][1,stokes_indx,:]
    dust.pixels = result["s"][2,stokes_indx,:]
    synchrotron.pixels = result["s"][3,stokes_indx,:]

    p1 = Plots.plot(cmb, title = "CMB_"*stokes_param, show = false)
    p2 = Plots.plot(dust, title = "Dust_"*stokes_param, show = false)
    p3 = Plots.plot(synchrotron, title = "Synchrotron_"*stokes_param, show = false)

    h1 = Plots.histogram(cmb[isfinite.(cmb)], normalize = true, show = false, legend = false)
    h2 = Plots.histogram(dust[isfinite.(dust)], normalize = true, show = false, legend = false)
    h3 = Plots.histogram(synchrotron[isfinite.(synchrotron)], normalize = true, show = false, legend = false)

    return Plots.plot(p1,p2,p3, h1,h2,h3, layout = (2,3), size = (1500, 500))
end

end # module thesis

