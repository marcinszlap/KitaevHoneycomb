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
    Nbonds = 3*Nx*Ny - Ny - Nx
    @show(Nbonds)
    @show(Nsites)

    latt = Lattice(undef, Nbonds)
    b = 0
    for i in 1:(Nsites-1)
        if mod(i,2) != 0
            latt[b+=1] = LatticeBond(i, i+1, 0,0,0,0, "xx")
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

function has_lattice_plot(dir::AbstractString)
    isdir(dir) || return false
    for f in readdir(dir; join=false)
        lf = lowercase(f)
        if endswith(lf, ".pdf") || endswith(lf, ".svg")
            return true
        end
    end
    return false
end
