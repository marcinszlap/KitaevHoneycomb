# DMRG settings container
struct DMRGSettings
    sweeps
    init_linkdims::Int
    outputlevel::Int
end

function make_sweeps(dmrgcfg)
    ns = dmrgcfg["sweeps"]
    s = Sweeps(ns)
    setmaxdim!(s, dmrgcfg["maxdim"]...)
    setcutoff!(s, dmrgcfg["cutoff"])
    return s
end

function dmrg_settings_from_cfg(dmrgcfg::Dict)
    return DMRGSettings(
        make_sweeps(dmrgcfg),
        dmrgcfg["init_linkdims"],
        dmrgcfg["outputlevel"],
    )
end

function DMRGcalc(os::OpSum, psi0, sites::AbstractVector; sweeps, outputlevel::Int=1)
    H = MPO(os, sites)
    obs = TruncErrObserver()

    t0 = time_ns()
    energy, psi = dmrg(H, psi0, sweeps; outputlevel=outputlevel, observer=obs)
    dmrg_time_s = (time_ns() - t0) / 1e9

    return energy, psi, obs.last_truncerr, dmrg_time_s
end

function EnergyGridSweep(
    theta_range::Vector{Float64},
    phi_range::Vector{Float64},
    FixedSite::Int,
    SpecialSite::Int,
    os::OpSum,
    state,
    sites::AbstractVector;
    dmrgset::DMRGSettings
)
    if SpecialSite != 0
        println("Grid sweep: rotating site $FixedSite and special site $SpecialSite")
    else
        println("Grid sweep: rotating site $FixedSite (no special site)")
    end

    E = Matrix{Float64}(undef, length(theta_range), length(phi_range))

    # Warm start MPS across grid points (huge speedup)
    psi = random_mps(sites, state; linkdims=dmrgset.init_linkdims)
    err_last = NaN
    dmrg_time_last = NaN

    for (it, theta) in enumerate(theta_range)
        for (ip, phi) in enumerate(phi_range)
            mgJx = sin(theta) * cos(phi)
            mgJy = sin(theta) * sin(phi)
            mgJz = cos(theta)

            ops = deepcopy(os)
            ops += mgJx, "Sx", FixedSite
            ops += mgJy, "Sy", FixedSite
            ops += mgJz, "Sz", FixedSite

            if SpecialSite != 0
                # fixed special orientation
                ops += sin(pi/2)*cos(0), "Sx", SpecialSite
                ops += sin(pi/2)*sin(0), "Sy", SpecialSite
                ops += cos(pi/2),        "Sz", SpecialSite
            end

            energy, psi, err_last, dmrg_time_s =
                DMRGcalc(ops, psi, sites; sweeps=dmrgset.sweeps, outputlevel=dmrgset.outputlevel)

            E[it, ip] = energy
            dmrg_time_last = dmrg_time_s
        end
    end

    return E, err_last, dmrg_time_last
end

function EnergySweep(
    sweep::Symbol,
    angle_range::Vector{Float64},
    fixed_angle::Float64,
    FixedSite::Int,
    SpecialSite::Int,
    os::OpSum,
    state,
    sites::AbstractVector;
    dmrgset::DMRGSettings
)
    if !(sweep in (:theta, :phi))
        psi0 = random_mps(sites, state; linkdims=dmrgset.init_linkdims)
        energy, _, err_last, dmrg_time_s =
            DMRGcalc(os, psi0, sites; sweeps=dmrgset.sweeps, outputlevel=dmrgset.outputlevel)
        return [energy], err_last, dmrg_time_s
    end

    if SpecialSite != 0
        println("Classical Spin initialization on rotating site $FixedSite and $SpecialSite")
    end

    data = Vector{Float64}(undef, length(angle_range))
    theta_vals = sweep == :theta ? collect(angle_range) : fill(fixed_angle, length(angle_range))
    phi_vals   = sweep == :phi   ? collect(angle_range) : fill(fixed_angle, length(angle_range))

    psi = random_mps(sites, state; linkdims=dmrgset.init_linkdims)
    err_last = NaN
    dmrg_time_last = NaN

    for i in eachindex(angle_range)
        theta = float(theta_vals[i])
        phi   = float(phi_vals[i])

        mgJx = sin(theta) * cos(phi)
        mgJy = sin(theta) * sin(phi)
        mgJz = cos(theta)

        ops = deepcopy(os)
        ops += mgJx, "Sx", FixedSite
        ops += mgJy, "Sy", FixedSite
        ops += mgJz, "Sz", FixedSite

        if SpecialSite != 0
            ops += sin(pi/2)*cos(0), "Sx", SpecialSite
            ops += sin(pi/2)*sin(0), "Sy", SpecialSite
            ops += cos(pi/2),        "Sz", SpecialSite
        end

        energy, psi, err_last, dmrg_time_s =
            DMRGcalc(ops, psi, sites; sweeps=dmrgset.sweeps, outputlevel=dmrgset.outputlevel)

        data[i] = energy
        dmrg_time_last = dmrg_time_s
    end

    println("End of calculations for rotating site $FixedSite and $SpecialSite")
    return data, err_last, dmrg_time_last
end
