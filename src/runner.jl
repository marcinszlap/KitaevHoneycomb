function run_one_job(cfg::Dict; config_path_for_plot::AbstractString="")
    Nx        = cfg["system"]["Nx"]
    Ny        = cfg["system"]["Ny"]
    spin_type = cfg["system"]["spin_type"]

    Jx = cfg["model"]["Jx"]
    Jy = cfg["model"]["Jy"]
    Jz = cfg["model"]["Jz"]

    yperiodic = cfg["lattice"]["yperiodic"]

    fixed_site        = cfg["sweep"]["fixed_site"]
    angle_name        = Symbol(cfg["sweep"]["angle_name"])
    angle_start       = parse_numeric(cfg["sweep"]["angle_start"])
    angle_stop        = parse_numeric(cfg["sweep"]["angle_stop"])
    angle_points      = cfg["sweep"]["angle_points"]
    other_angle       = parse_numeric(cfg["sweep"]["other_angle"])
    special_sites_cfg = cfg["sweep"]["special_sites"]
    sweep_mode        = lowercase(get(cfg["sweep"], "sweep_mode", "1d"))

    outdir  = cfg["io"]["outfile_dir"]
    prefix  = cfg["io"]["outfile_prefix"]

    dmrgset = dmrg_settings_from_cfg(cfg["DMRG"])

    Nsites      = 2 * Ny * Nx
    angle_range = collect(range(angle_start, stop=angle_stop, length=angle_points))

    size_dir = joinpath(outdir, "$(Nx)x$(Ny)")
    mkpath(size_dir)

    # ---------- lattice plot (portable) ----------
    if has_lattice_plot(size_dir)
        println("Lattice plot found in $size_dir")
    else
        project_root = normpath(joinpath(@__DIR__, ".."))
        plot_script  = joinpath(project_root, "lattice_plot.jl")

        plot_cfg_path = joinpath(size_dir, "_plot_config.toml")
        plot_cfg = Dict{String,Any}(
            "Nx" => Nx,
            "Ny" => Ny,
            "system"    => Dict("Nx" => Nx, "Ny" => Ny),
            "lattice"   => Dict("yperiodic" => yperiodic),
            "yperiodic" => yperiodic,
        )
        open(plot_cfg_path, "w") do io
            TOML.print(io, plot_cfg)
        end

        run(`julia --project=$project_root $plot_script --config $plot_cfg_path --outdir $size_dir`)
    end

    # ---------- lattice & Hamiltonian ----------
    my_lattice = HoneyCombLattice(Nx, Ny, yperiodic=yperiodic)

    os_base = OpSum()
    for b in my_lattice
        if b.bond_type == "xx"
            os_base += Jx, "Sx", b.site1, "Sx", b.site2
        elseif b.bond_type == "yy"
            os_base += Jy, "Sy", b.site1, "Sy", b.site2
        elseif b.bond_type == "zz"
            os_base += Jz, "Sz", b.site1, "Sz", b.site2
        end
    end

    sites = siteinds(spin_type, Nsites)
    state = [isodd(n) ? "Dn" : "Up" for n in 1:Nsites]

    # ---------- benchmark / no-sweep ----------
    if !(angle_name in (:theta, :phi)) && sweep_mode != "grid"
        psi0 = random_mps(sites, state; linkdims=dmrgset.init_linkdims)
        energy, _, last_truncerr, dmrg_time_s =
            DMRGcalc(os_base, psi0, sites; sweeps=dmrgset.sweeps, outputlevel=dmrgset.outputlevel)

        outfile_name = "$(prefix)_bench_$(Nx)x$(Ny).dat"
        OutFile = joinpath(size_dir, outfile_name)
        open(OutFile, "w") do io
            println(io, "Nx\tNy\tNsites\tEnergy\tLastTruncErr\tDMRGTime_s")
            @printf(io, "%d\t%d\t%d\t%.12f\t%.3e\t%.3f\n", Nx, Ny, Nsites, energy, last_truncerr, dmrg_time_s)
        end

        return (Nx=Nx, Ny=Ny, Nsites=Nsites,
                energy=energy, last_truncerr=last_truncerr,
                dmrg_time_s=dmrg_time_s, outfile=OutFile)
    end

    # ---------- special sites list ----------
    PossibleSites = normalize_special_sites(special_sites_cfg, fixed_site, Nsites)

    # ============================================================
    # GRID MODE: compute E(theta, phi) and early return
    # ============================================================
    if sweep_mode == "grid"
        # Parse grid ranges *here* (fail early if missing keys)
        theta_start  = parse_numeric(cfg["sweep"]["theta_start"])
        theta_stop   = parse_numeric(cfg["sweep"]["theta_stop"])
        theta_points = cfg["sweep"]["theta_points"]

        phi_start    = parse_numeric(cfg["sweep"]["phi_start"])
        phi_stop     = parse_numeric(cfg["sweep"]["phi_stop"])
        phi_points   = cfg["sweep"]["phi_points"]

        theta_range = collect(range(theta_start, stop=theta_stop, length=theta_points))
        phi_range   = collect(range(phi_start,   stop=phi_stop,   length=phi_points))

        outfile_name = "$(prefix)_$(Nx)x$(Ny)_grid_fix$(fixed_site).dat"
        OutFile = joinpath(size_dir, outfile_name)

        Emats    = Vector{Matrix{Float64}}(undef, length(PossibleSites))
        lasterrs = Vector{Float64}(undef, length(PossibleSites))
        times    = Vector{Float64}(undef, length(PossibleSites))

        @threads for idx in eachindex(PossibleSites)
            s = PossibleSites[idx]
            E, err_last, t_last = EnergyGridSweep(theta_range, phi_range,
                                                  fixed_site, s,
                                                  os_base, state, sites;
                                                  dmrgset=dmrgset)
            Emats[idx]    = E
            lasterrs[idx] = err_last
            times[idx]    = t_last
        end

        open(OutFile, "w") do io
            println(io, "# grid sweep: rows=theta, cols=phi")
            println(io, "# Nx=$Nx Ny=$Ny fixed_site=$fixed_site")
            println(io, "# theta: start=$(theta_range[1]) stop=$(theta_range[end]) points=$(length(theta_range))")
            println(io, "# phi:   start=$(phi_range[1]) stop=$(phi_range[end]) points=$(length(phi_range))")
            println(io, "# special_sites=$(PossibleSites)")
            println(io)

            for (idx, s) in enumerate(PossibleSites)
                println(io, "# --- SpecialSite = $s ---")
                @printf(io, "theta\\phi")
                for phi in phi_range
                    @printf(io, "\t%g", phi)
                end
                println(io)

                E = Emats[idx]
                for (it, theta) in enumerate(theta_range)
                    @printf(io, "%g", theta)
                    for ip in eachindex(phi_range)
                        @printf(io, "\t%.10f", E[it, ip])
                    end
                    println(io)
                end
                println(io)
            end
        end

        println("Grid data written to $OutFile")

        return (Nx=Nx, Ny=Ny, Nsites=Nsites,
                energy=Emats[end][end, end],
                last_truncerr=lasterrs[end],
                dmrg_time_s=times[end],
                outfile=OutFile)
    end

    # ============================================================
    # 1D MODE: compute E(angle) for each special site
    # ============================================================
    outfile_name = "$(prefix)_$(Nx)x$(Ny)_$(angle_name)_fix$(fixed_site).dat"
    OutFile = joinpath(size_dir, outfile_name)

    results  = Vector{Vector{Float64}}(undef, length(PossibleSites))
    lasterrs = Vector{Float64}(undef, length(PossibleSites))
    times    = Vector{Float64}(undef, length(PossibleSites))

    @threads for idx in eachindex(PossibleSites)
        s = PossibleSites[idx]
        data, err_last, t_last = EnergySweep(angle_name, angle_range, other_angle,
                                             fixed_site, s,
                                             os_base, state, sites;
                                             dmrgset=dmrgset)
        results[idx]  = data
        lasterrs[idx] = err_last
        times[idx]    = t_last
    end

    open(OutFile, "w") do io
        @printf(io, "SpecialSite\\%s", string(angle_name))
        for angle in angle_range
            @printf(io, "\t%g", angle)
        end
        println(io)

        for (idx, s) in enumerate(PossibleSites)
            @printf(io, "%d", s)
            for e in results[idx]
                @printf(io, "\t%.8f", e)
            end
            println(io)
        end
    end

    println("Data written to $OutFile")

    return (Nx=Nx, Ny=Ny, Nsites=Nsites,
            energy=results[end][end],
            last_truncerr=lasterrs[end],
            dmrg_time_s=times[end],
            outfile=OutFile)
end


function run_job_dict(cfg::Dict; summary_path::AbstractString="results/benchmark.tsv",
                      job_id::Int=0, config_path_for_plot::AbstractString="")
    t0 = time_ns()
    ok = true
    errtxt = ""

    Nx = Ny = Nsites = -1
    energy = NaN
    last_truncerr = NaN
    outfile = ""
    dmrg_time_s = NaN

    try
        res = run_one_job(cfg; config_path_for_plot=config_path_for_plot)
        Nx, Ny, Nsites = res.Nx, res.Ny, res.Nsites
        energy, last_truncerr = res.energy, res.last_truncerr
        outfile, dmrg_time_s = res.outfile, res.dmrg_time_s
    catch e
        ok = false
        errtxt = sprint(showerror, e, catch_backtrace())
    end

    elapsed_s = (time_ns() - t0) / 1e9
    timestamp = Dates.format(now(), dateformat"yyyy-mm-dd HH:MM:SS")

    mkpath(dirname(summary_path))
    newfile = !isfile(summary_path)

    open(summary_path, "a") do io
        if newfile
            println(io, "timestamp\tjob_id\tok\tNx\tNy\tNsites\tdmrg_time_s\tenergy\tlast_truncerr\toutfile\telapsed_s\terror")
        end
        println(io, join((
            timestamp, string(job_id), string(ok),
            string(Nx), string(Ny), string(Nsites),
            @sprintf("%.3f", dmrg_time_s),
            @sprintf("%.12f", energy),
            @sprintf("%.3e", last_truncerr),
            outfile,
            @sprintf("%.3f", elapsed_s),
            replace(errtxt, '\n' => ' ')
        ), '\t'))
    end

    return ok
end

function run_batch(config_path::AbstractString; summary_path="results/benchmark.tsv")
    cfg = TOML.parsefile(config_path)

    defaults = Dict(
        "DMRG"    => cfg["DMRG"],
        "model"   => cfg["model"],
        "lattice" => cfg["lattice"],
        "system"  => cfg["system"],
        "sweep"   => cfg["sweep"],
        "io"      => cfg["io"],
    )

    jobs = get(cfg, "jobs", Any[Dict{String,Any}()])

    for (idx, job_any) in enumerate(jobs)
        job = Dict{String,Any}(job_any)
        eff = deepmerge(defaults, job_to_nested(job))
        println("\n=== Job $idx ===")
        run_job_dict(eff; summary_path=summary_path, job_id=idx, config_path_for_plot=config_path)
    end
end
