function parse_numeric(x)
    x isa Number && return float(x)
    x isa String || error("Invalid numeric field: $x")
    s = lowercase(strip(x))

    s == "pi"   && return pi
    s == "2pi"  && return 2pi
    s == "pi/2" && return pi/2
    s == "pi/4" && return pi/4

    ex = Meta.parse(s)
    return float(eval(:(let pi = $pi; $ex end)))
end

function normalize_special_sites(x, fixed_site::Int, Nsites::Int)::Vector{Int}
    if x isa String
        s = lowercase(strip(x))
        if s == "none"
            return Int[0]  # single run, no special site
        elseif s == "all"
            return setdiff(collect(1:Nsites), [fixed_site])
        else
            error("special_sites must be 'none', 'all', or an integer array; got '$x'")
        end
    elseif x isa AbstractVector
        v = unique(Int.(x))
        any(v .< 1)      && error("special_sites contains site < 1")
        any(v .> Nsites) && error("special_sites contains site > Nsites")
        fixed_site in v  && error("special_sites contains fixed_site=$fixed_site (would double-field one site)")
        return v
    else
        error("special_sites must be String or Array; got $(typeof(x))")
    end
end

function deepmerge(a::Dict, b::Dict)
    out = copy(a)
    for (k, v) in b
        if haskey(out, k) && out[k] isa Dict && v isa Dict
            out[k] = deepmerge(out[k], v)
        else
            out[k] = v
        end
    end
    return out
end

function job_to_nested(job::Dict)
    out = Dict{String,Any}()

    sys = Dict{String,Any}()
    for k in ("Nx","Ny","spin_type")
        haskey(job,k) && (sys[k] = job[k])
    end
    !isempty(sys) && (out["system"] = sys)

    lat = Dict{String,Any}()
    haskey(job,"yperiodic") && (lat["yperiodic"] = job["yperiodic"])
    !isempty(lat) && (out["lattice"] = lat)

    sw = Dict{String,Any}()
    for k in (
        "fixed_site",
        "angle_name", "angle_start", "angle_stop", "angle_points", "other_angle",
        "special_sites",
        "sweep_mode",
        "theta_start", "theta_stop", "theta_points",
        "phi_start", "phi_stop", "phi_points",
    )
        haskey(job, k) && (sw[k] = job[k])
    end
    !isempty(sw) && (out["sweep"] = sw)

    mdl = Dict{String,Any}()
    for k in ("Jx","Jy","Jz")
        haskey(job,k) && (mdl[k] = job[k])
    end
    !isempty(mdl) && (out["model"] = mdl)

    return out
end
