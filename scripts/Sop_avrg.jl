using ITensors, ITensorMPS
using Printf
using TOML
using Base.Threads
using LinearAlgebra
using DelimitedFiles
using Dates
using CUDA
using Random

CUDA.allowscalar(false)
function save_spin_component(path::AbstractString,
                             vals_all::Vector{Vector{Float64}},
                             mg_sites::AbstractVector{<:Integer},
                             Nsites::Int)
    ncfg = length(mg_sites)
    @assert length(vals_all) == ncfg

        data = Matrix{Float64}(undef, Nsites, 1 + ncfg)
    data[:, 1] .= collect(1:Nsites)
    for (k, _) in enumerate(mg_sites)
        @assert length(vals_all[k]) == Nsites
        data[:, 1 + k] .= vals_all[k]
    end

    open(path, "w") do io
        println(io, "# i\t", join(mg_sites, "\t"))
        writedlm(io, data, '\t')
    end
end

struct LatticeBond
    site1::Int
    site2::Int
    param1::Float64
    param2::Float64
    param3::Float64
    param4::Float64
    bond_type::String
end

const Lattice = Vector{LatticeBond}

mutable struct QualityObserver <: AbstractObserver
    energy_tol::Float64
    last_energy::Float64
    spin_tol::Float64

    last_sx::Vector{Float64}
    last_sy::Vector{Float64}
    last_sz::Vector{Float64}

    converged_at_sweep::Int
    extra_sweeps::Int

    QualityObserver(; energy_tol=0.0, spin_tol=0.0, extra_sweeps=2) =
        new(energy_tol, 1e9, spin_tol,
            Float64[], Float64[], Float64[],
            0, extra_sweeps)
end

function ITensorMPS.checkdone!(o::QualityObserver; kwargs...)
    sw = kwargs[:sweep]
    energy = kwargs[:energy]
    psi = kwargs[:psi]

    sx = expect(psi, "Sx")
    sy = expect(psi, "Sy")
    sz = expect(psi, "Sz")

    energy_conv = abs(energy - o.last_energy)

    spin_conv = Inf
    if !isempty(o.last_sx)
        spin_conv =
            maximum(abs.(sx .- o.last_sx)) +
            maximum(abs.(sy .- o.last_sy)) +
            maximum(abs.(sz .- o.last_sz))
    end

    energy_good = energy_conv < o.energy_tol
    spin_good = spin_conv < o.spin_tol

    converged_now = energy_good && spin_good
    println("Sweep $sw | ?E=$energy_conv | ?Spin=$spin_conv")

    o.last_energy = energy
    o.last_sx = copy(sx)
    o.last_sy = copy(sy)
    o.last_sz = copy(sz)

    if converged_now && o.converged_at_sweep == 0
        o.converged_at_sweep = sw
    end

    return false
end

function site_index_snake(x::Int, y::Int, sub::Int, Nx::Int, Ny::Int)
    @assert 1 <= x <= Nx
    @assert 1 <= y <= Ny
    @assert sub == 1 || sub == 2

    col_offset = 2 * Ny * (x - 1)

    if isodd(x)
        yeff = y
        subeff = sub
    else
        yeff = Ny - y + 1
        subeff = 3 - sub   # 1->2, 2->1
    end

    return col_offset + 2 * (yeff - 1) + subeff
end

function site_index_snake_reverse(x, y, sub, Nx, Ny)
    Nsites = 2 * Nx * Ny
    i = site_index_snake(x, y, sub, Nx, Ny)
    return Nsites - i + 1
end

function energy_variance(H::MPO, psi::MPS)
    E = inner(psi', H, psi)
    H2 = inner(H, psi, H, psi)
    return real(H2 - E^2)
end

function HoneyCombLattice(Nx::Int, Ny::Int; yperiodic=false)::Lattice
    Nsites = 2*Ny*Nx
    if yperiodic
        Nbonds = 3*Nx*Ny - Ny
    else
        Nbonds = 3*Nx*Ny - Ny - Nx
    end
    @show(Nbonds)
    @show(Nsites)

    latt = LatticeBond[]

    A(x,y) = site_index_snake_reverse(x,y,1,Nx,Ny)
    B(x,y) = site_index_snake_reverse(x,y,2,Nx,Ny)
    for x in 1:Nx
        for y in 1:Ny
            push!(latt, LatticeBond(A(x,y), B(x,y), 0.0, 0.0, 0.0, 0.0, "xx"))

            if x < Nx
                push!(latt, LatticeBond(B(x,y), A(x+1,y), 0.0, 0.0, 0.0, 0.0, "yy"))
            end

            if y < Ny
                push!(latt, LatticeBond(B(x,y), A(x,y+1), 0.0, 0.0, 0.0, 0.0, "zz"))
            elseif yperiodic
                push!(latt, LatticeBond(B(x,y), A(x,1), 0.0, 0.0, 0.0, 0.0, "zz"))
            end
        end
    end
    return latt
end
let
    ####### CONFIGURATION #######
    Nx = 5
    Ny = 3
    Nsites = 2*Ny*Nx

    Jx = -1.0
    Jy = -1.0
    Jz = -1.0

    theta = pi/2
    phi = 0.0

    second_spin = false
    s2_Pos = 2

    #[i for i in 1:Nsites if i != s2_Pos]
    mg_Site = 15
    yperiodic = false

    sweeps = Sweeps(40)
    maxdim!(sweeps,
        50,100,150,200,300,400,
        500,700,900,1200,1600,2000,
        2500,3000,3000,3000,3000,3000,
        3000,3000,3000,3000,3000,3000,
        3000,3000,3000,3000,3000,3000,
        3000,3000,3000,3000,3000,3000,
        3000,3000,3000,3000
    )
    cutoff!(sweeps, 1e-8)
    noise!(sweeps,
        1e-5,1e-5,1e-6,1e-6,1e-7,1e-7,
        1e-8,1e-8,1e-9,1e-9,1e-10,1e-10,
        0.0
    )    
    # observer checking for quality of energy and spin measurements
    # energy tolerance, spin tolerance
    obs = QualityObserver(energy_tol=1E-7, spin_tol=1E-5, extra_sweeps=2)

    ####### DATA #######
    sx_all = Vector{Vector{Float64}}()
    sy_all = Vector{Vector{Float64}}()
    sz_all = Vector{Vector{Float64}}()
    E_all  = Float64[]
    ####### LATTICE DEFINITION ########
    my_lattice = HoneyCombLattice(Nx, Ny, yperiodic=yperiodic)

    Nbonds = 3*Nx*Ny - Ny - Nx
    for b in my_lattice
        println("$(b.bond_type): $(b.site1) <-> $(b.site2)")
    end

    for x in 1:Nx
    println("column x = $x")
        for y in 1:Ny
            a = site_index_snake(x,y,1,Nx,Ny)
            b = site_index_snake(x,y,2,Nx,Ny)
            println("cell ($x,$y): A=$a  B=$b")
        end
        println()
    end

    for t in ["xx","yy","zz"]
        ds = [abs(b.site1 - b.site2) for b in my_lattice if b.bond_type == t]
        println("\n$t")
        println("count = ", length(ds))
        println("min   = ", minimum(ds))
        println("max   = ", maximum(ds))
        
    end

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
    
    ####### ADDING MAGNETIC SPIN #######
    mgJx = sin(theta) * cos(phi)
    mgJy = sin(theta) * sin(phi)
    mgJz = cos(theta)
    
    os = deepcopy(os_base)
    os += mgJx, "Sx", mg_Site
    os += mgJy, "Sy", mg_Site
    os += mgJz, "Sz", mg_Site
    if second_spin
        os += mgJx, "Sx", s2_Pos
        os += mgJy, "Sy", s2_Pos
        os += mgJz, "Sz", s2_Pos
    end
    sites = siteinds("S=1/2", Nsites)
    ######## DMRG RUNS INICIALIZATION ########

    H = MPO(os, sites)
    H = cu(H)

    # ===== GROUND STATE =====
    obs0 = QualityObserver(energy_tol=1e-7, spin_tol=1e-5, extra_sweeps=2)

    psi0_init = random_mps(sites; linkdims=10)
    psi0_init = cu(psi0_init)

    E0, psi0 = dmrg(H, psi0_init, sweeps; observer=obs0)

    var0 = energy_variance(H, psi0)

    println("\n=== GROUND STATE ===")
    println("E0      = ", E0)
    println("var0    = ", var0)
    println("maxdim0 = ", maxlinkdim(psi0))

    # ===== FIRST ORTHOGONAL STATE =====
    obs1 = QualityObserver(energy_tol=1e-7, spin_tol=1e-5, extra_sweeps=2)

    psi1_init = random_mps(sites; linkdims=10)
    psi1_init = cu(psi1_init)

    weight = 20.0
    E1, psi1 = dmrg(H, [psi0], psi1_init, sweeps; weight=weight, observer=obs1)

    var1 = energy_variance(H, psi1)

    println("\n=== FIRST ORTHOGONAL STATE ===")
    println("E1      = ", E1)
    println("var1    = ", var1)
    println("maxdim1 = ", maxlinkdim(psi1))
    
    println("gap01   = ", E1 - E0)
    println("|<0|1>| = ", abs(inner(psi0, psi1)))
    push!(energies, energy)
    push!(variances, var)
    push!(psis, psi)

    if E1 < E0
        println("Orthogonal state has lower energy, swapping labels.")
        E0, E1 = E1, E0
        psi0, psi1 = psi1, psi0
        var0, var1 = var1, var0
    end
    energy = E0
    psi = psi0

    push!(E_all, energy)
    push!(sx_all, expect(psi, "Sx"))
    push!(sy_all, expect(psi, "Sy"))
    push!(sz_all, expect(psi, "Sz"))

    script_dir = @__DIR__

    outdir = joinpath(
        "results",
        "$(Nx)x$(Ny)_PBC=$(yperiodic)_Jx=$(Jx)_Jy=$(Jy)_Jz=$(Jz)",
        "Savrg_theta=$(theta)_phi=$(phi)_$(mg_Site)"
    )
    mkpath(outdir)

    println("PWD = ", pwd())
    println("Script dir = ", script_dir)
    println("Writing to = ", outdir)
    flush(stdout)

    mg_sites = [mg_Site]

    if second_spin
        save_spin_component(joinpath(outdir, "Sx_2spins_fixed$(s2_Pos).dat"), sx_all, mg_sites, Nsites)
        save_spin_component(joinpath(outdir, "Sy_2spins_fixed$(s2_Pos).dat"), sy_all, mg_sites, Nsites)
        save_spin_component(joinpath(outdir, "Sz_2spins_fixed$(s2_Pos).dat"), sz_all, mg_sites, Nsites)
    else
        save_spin_component(joinpath(outdir, "Sx.dat"), sx_all, mg_sites, Nsites)
        save_spin_component(joinpath(outdir, "Sy.dat"), sy_all, mg_sites, Nsites)
        save_spin_component(joinpath(outdir, "Sz.dat"), sz_all, mg_sites, Nsites)
    end
end