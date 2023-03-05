only pointing accuracy, no observed signal

t_start = DateTime(2022, 1, 1, 12, 0, 0)
tel_ang = Sl.TelescopeAngles(wheel2ang_0_rad = deg2rad(1.0/60))
days = [5, 50, 100, 200, 400, 600]
setups = [PRMaps.Setup(sampling_freq_Hz = 50., total_time_s = 24. * 3600. * i) for i in days]    
result = [pointing_error_arcmin(Sl.CameraAngles(), tel_ang, setup, t_start) for setup in setups]