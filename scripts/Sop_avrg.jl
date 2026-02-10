using ITensors, ITensorMPS
using Printf
using TOML
using Base.Threads
using LinearAlgebra
using DelimitedFiles
using Dates

function save_spin_component(path::AbstractString,
                             vals_all::Vector{Vector{Float64}},
                             mg_sites::Vector{Int},
                             Nsites::Int)
    ncfg = length(mg_sites)
    @assert length(vals_all) == ncfg

    # macierz: [i  vals(i; s=mg_sites[1])  vals(i; s=mg_sites[2]) ...]
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

function HoneyCombLattice(Nx::Int, Ny::Int; yperiodic=false)::Lattice
    Nsites = 2*Ny*Nx
    if yperiodic
        Nbonds = 3*Nx*Ny - Ny
    else
        Nbonds = 3*Nx*Ny - Ny - Nx
    end
    @show(Nbonds)
    @show(Nsites)

    latt = Lattice(undef, Nbonds)
    b = 0
    for i in 1:(Nsites-1)
        if mod(i,2) != 0
            latt[b+=1] = LatticeBond(i, i+1, 0,0,0,0, "xx")
            if yperiodic && (mod(i, 2*Ny) == 1)
                latt[b+=1] = LatticeBond(i, i+2*Ny-1, 0,0,0,0, "zz")
            end
        else
            if i <= 2*Ny*(Nx-1)
                latt[b+=1] = LatticeBond(i, i+2*Ny-1, 0,0,0,0, "yy")
            end
            if mod(i, 2*Ny) != 0
                latt[b+=1] = LatticeBond(i, i+1, 0,0,0,0, "zz")
            end
        end
    end
    return latt
end
let
    ####### CONFIGURATION #######
    Nx = 3
    Ny = 5
    Nsites = 2*Ny*Nx

    Jx = -1.0
    Jy = -1.0
    Jz = -1.0

    theta = 0.0
    phi = 0.0

    second_spin = false
    s2_Pos = 2

    #[i for i in 1:Nsites if i != s2_Pos]
    mg_Sites = [3, 16, 31, 38, 48, 80, 95]
    yperiodic = true

    sweeps = Sweeps(5)
    maxdim!(sweeps, 100,200,300,400,500)
    cutoff!(sweeps, 1e-10)

    ####### DATA #######
    sx_all = Vector{Vector{Float64}}()
    sy_all = Vector{Vector{Float64}}()
    sz_all = Vector{Vector{Float64}}()
    E_all  = Float64[]
    ####### LATTICE DEFINITION ########
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

    sites = siteinds("S=1/2", Nsites)
    state = [isodd(n) ? "Dn" : "Up" for n in 1:Nsites]
    psi0 = MPS(sites, state)

    ####### ADDING MAGNETIC SPIN #######
    mgJx = sin(theta) * cos(phi)
    mgJy = sin(theta) * sin(phi)
    mgJz = cos(theta)
    for s in mg_Sites
        os = deepcopy(os_base)
        os += mgJx, "Sx", s
        os += mgJy, "Sy", s
        os += mgJz, "Sz", s
        if second_spin
            os += mgJx, "Sx", s2_Pos
            os += mgJy, "Sy", s2_Pos
            os += mgJz, "Sz", s2_Pos
        end
        H = MPO(os, sites)
        energy, psi  =  dmrg(H, psi0, sweeps)
        
        push!(E_all, energy)
        push!(sx_all, expect(psi, "Sx"))
        push!(sy_all, expect(psi, "Sy"))
        push!(sz_all, expect(psi, "Sz"))
    end

    script_dir = @__DIR__                               # .../marcins/KitaevHoneycomb/scripts
    marcins_dir = normpath(joinpath(script_dir, "..", ".."))  # .../marcins

    outdir = joinpath("results", "$(Nx)x$(Ny)_PBC=$(yperiodic)", "Savrg")
    mkpath(outdir)
    println("PWD = ", pwd())
    println("Script dir = ", script_dir)
    println("Writing to = ", outdir)
    flush(stdout)

    if second_spin
        save_spin_component(joinpath(outdir, "Sx_2spins_fixed$(s2_Pos).dat"), sx_all, mg_Sites, Nsites)
        save_spin_component(joinpath(outdir, "Sy_2spins_fixed$(s2_Pos).dat"), sy_all, mg_Sites, Nsites)
        save_spin_component(joinpath(outdir, "Sz_2spins_fixed$(s2_Pos).dat"), sz_all, mg_Sites, Nsites)
    else
        save_spin_component(joinpath(outdir, "Sx.dat"), sx_all, mg_Sites, Nsites)
        save_spin_component(joinpath(outdir, "Sy.dat"), sy_all, mg_Sites, Nsites)
        save_spin_component(joinpath(outdir, "Sz.dat"), sz_all, mg_Sites, Nsites)
    end
end