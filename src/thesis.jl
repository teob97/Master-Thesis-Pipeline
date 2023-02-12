module thesis

using PyCall
using PRMaps
using Healpix
using StatsPlots
using Pandas
using Random, Distributions
using Dates
import Stripeline as Sl

export get_foreground_maps, get_noise_maps, get_observations
export fgbuster_basic_comp_sep
export run_fgbuster, run_fgbuster_with_error, get_corrplot, get_map_and_hist

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

    def get_freq_maps(instruments, sky_model, nside, unit):
        emissions = fgbuster.get_observation(instruments, sky_model, noise = False, nside = nside, unit = unit)
        for i,emission in enumerate(emissions):
            emissions[i] = pysm3.apply_smoothing_and_coord_transform(emission, rot=hp.Rotator(coord=("G", "C")))
        return emissions

    def get_noise_maps(instruments, nside, unit):
        return fgbuster.get_noise_realization(nside, instruments, unit)

    def fgbuster_pipeline_c1s0d0(instruments, data):
        components = [fgbuster.CMB(), fgbuster.Dust(353.), fgbuster.Synchrotron(23.)] 
        return fgbuster.basic_comp_sep(components, instruments, data[:,1:])

    def fgbuster_pipeline_c1s3d0(instruments, data):
        components = [fgbuster.CMB(), fgbuster.Dust(353.), fgbuster.Synchrotron(23., running=None, nu_pivot=23.)] 
        return fgbuster.basic_comp_sep(components, instruments, data[:,1:])

    def get_mvDistibution(result):
        return np.random.multivariate_normal(result['x'], result['Sigma'], 10000)
    
    """
end

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

function vec2map(vec, nside)
    maps = Healpix.PolarizedHealpixMap[]
    for indx in axes(vec,1)
        map = PolarizedHealpixMap{Float64, RingOrder}(nside)
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
    maps = vec2map(results, nside)
    return maps
end

function get_noise_maps(
    instruments,
    nside::Int;
    unit::String = "uK_CMB" 
)
    results = py"get_noise_maps"(instruments, nside, unit)
    maps = vec2map(results, nside)
    return maps
end

function fgbuster_basic_comp_sep_c1s0d0(instruments, maps)
    data = map2vec(maps)
    return py"fgbuster_pipeline_c1s0d0"(instruments, data)    
end

function fgbuster_basic_comp_sep_c1s3d0(instruments, maps)
    data = map2vec(maps)
    return py"fgbuster_pipeline_c1s3d0"(instruments, data)    
end


# Simulation without poining errors
function run_fgbuster(
    instruments,
    cam_ang :: Sl.CameraAngles,
    setup :: PRMaps.Setup,
    sky_model :: String,
    nside :: Int,
    t_start :: Dates.DateTime
)
    # "Generating sky signal and noise"
    signals = get_foreground_maps(instruments, sky_model, nside)
    noise = get_noise_maps(instruments, nside)

    # "Simulating telescope scanning strategy [LSPE/Strip]"
    observations, _ = PRMaps.makeIdealMapsIQU(cam_ang, signals, setup, t_start)

    # "Adding white noise based on the instruments sensitivity"
    for indx in axes(observations, 1)
        observations[indx].i.pixels += noise[indx].i.pixels
        observations[indx].q.pixels += noise[indx].q.pixels 
        observations[indx].u.pixels += noise[indx].u.pixels 
    end

    # "Run component separation using fgbuster"
    if sky_model == "c1s0d0"
        return fgbuster_basic_comp_sep_c1s0d0(instruments, observations)
    elseif sky_model == "c1s3d0"
        return fgbuster_basic_comp_sep_c1s3d0(instruments, observations)
    else
        error("Invalid sky model")
    end

end

function run_fgbuster_with_error(
    instruments,
    cam_ang :: Sl.CameraAngles,
    tel_ang :: Sl.TelescopeAngles,
    setup :: PRMaps.Setup,
    sky_model :: String,
    nside :: Int,
    t_start :: Dates.DateTime
)
    # Separate LSPE/Strip experiment from the others
    if !any(instruments.instrument .== "LSPE/Strip")
        error("Invalid simulation: no Instrument named LSPE/Strip found")
    end
    # Remove LSPE/Strip from the instruments vector and store them in another vector 
    strip = query(instruments, :(instrument=="LSPE/Strip"))
    instruments = query(instruments, :(instrument!="LSPE/Strip"))
    
    # "Generating sky signal"
    signals = get_foreground_maps(instruments, sky_model, nside)
    signals_strip = get_foreground_maps(strip, sky_model, nside)

    # "Simulating telescope scanning strategy [LSPE/Strip WITH pointing error]"
    observations, _ = makeIdealMapsIQU(cam_ang, signals, setup, t_start)
    observations_strip, _ = makeErroredMapsIQU(cam_ang, tel_ang, signals_strip, setup, t_start)

    # Append LSPE/Strip result to the other results 
    instruments = Pandas.concat([instruments, strip])
    append!(observations, observations_strip)

    # "Adding white noise based on the instruments sensitivity"
    noise = get_noise_maps(instruments, nside)

    for indx in axes(observations, 1)
        observations[indx].i.pixels += noise[indx].i.pixels
        observations[indx].q.pixels += noise[indx].q.pixels 
        observations[indx].u.pixels += noise[indx].u.pixels 
    end

    # "Run component separation using fgbuster"
    if sky_model == "c1s0d0"
        return fgbuster_basic_comp_sep_c1s0d0(instruments, observations)
    elseif sky_model == "c1s3d0"
        return fgbuster_basic_comp_sep_c1s3d0(instruments, observations)
    else
        error("Invalid sky model")
    end
    
end

function get_corrplot(result)
    sampling = py"get_mvDistibution"(result)
    return StatsPlots.corrplot(sampling, label = result["params"], fillcolor =:thermal)
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

    return (p1,p2,p3),(h1,h2,h3)
end


end # module thesis

