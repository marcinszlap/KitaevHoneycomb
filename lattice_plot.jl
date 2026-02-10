using TOML
using Printf

# --------------------------
# Same bond structure as your HoneyCombLattice
# --------------------------
struct LatticeBond
    site1::Int
    site2::Int
    bond_type::String
end

function honeycomb_bonds(Nx::Int, Ny::Int; yperiodic::Bool=false)
    Nsites = 2 * Nx * Ny
    Nbonds = 3 * Nx * Ny - Ny - Nx  # matches your count for open boundaries
    bonds = LatticeBond[]
    sizehint!(bonds, Nbonds)

    b = 0
    for i in 1:(Nsites-1)
        if mod(i, 2) != 0
            push!(bonds, LatticeBond(i, i + 1, "xx"))
            b += 1
        else
            if i <= 2 * Ny * (Nx - 1)
                push!(bonds, LatticeBond(i, i + 2 * Ny - 1, "yy"))
                b += 1
            end
            if mod(i, 2 * Ny) != 0
                push!(bonds, LatticeBond(i, i + 1, "zz"))
                b += 1
            end
        end
    end

    # NOTE: yperiodic is not implemented in your original bond loop either,
    # so we keep it unused here for consistency.
    return Nsites, bonds
end

# --------------------------
# Coordinates consistent with your indexing
#
# For each column x = 1..Nx and row r = 1..Ny:
#   A(x,r) index = (x-1)*2Ny + (2r-1)  (odd)
#   B(x,r) index = A+1                 (even)
#
# Place A and B on a standard honeycomb geometry so that:
#   A-B           (xx) has length 1
#   B(x,r)-A(x,r+1) (zz) has length 1
#   B(x,r)-A(x+1,r) (yy) has length 1
# --------------------------
function honeycomb_positions(Nx::Int, Ny::Int)
    Nsites = 2 * Nx * Ny
    pos = Vector{Tuple{Float64,Float64}}(undef, Nsites)
    placed = falses(Nsites)

    a = 1.0

    # Desired directed bond vectors (site1 -> site2)
    vx = (sqrt(3) * a,  a)   # for "xx"
    vz = (0.0,          2a)  # for "zz"
    vy = (sqrt(3) * a, -a)   # for "yy"

    # --- Rebuild the same bond list as in HoneyCombLattice ---
    # store as tuples: (s1, s2, bond_type)
    bonds = Vector{Tuple{Int,Int,String}}()
    for i in 1:(Nsites - 1)
        if mod(i, 2) != 0
            push!(bonds, (i, i + 1, "xx"))
        else
            if i <= 2 * Ny * (Nx - 1)
                push!(bonds, (i, i + 2 * Ny - 1, "yy"))
            end
            if mod(i, 2 * Ny) != 0
                push!(bonds, (i, i + 1, "zz"))
            end
        end
    end

    # --- Build adjacency with vectors ---
    adj = [Vector{Tuple{Int,Tuple{Float64,Float64}}}() for _ in 1:Nsites]

    bondvec(bt::String) = bt == "xx" ? vx : bt == "yy" ? vy : vz

    for (s1, s2, bt) in bonds
        v = bondvec(bt)
        push!(adj[s1], (s2, v))
        push!(adj[s2], (s1, (-v[1], -v[2])))
    end

    # --- BFS placement ---
    pos[1] = (0.0, 0.0)
    placed[1] = true
    queue = [1]

    while !isempty(queue)
        s = popfirst!(queue)
        xs, ys = pos[s]

        for (t, v) in adj[s]
            if !placed[t]
                pos[t] = (xs + v[1], ys + v[2])
                placed[t] = true
                push!(queue, t)
            end
        end
    end

    # If anything remained unplaced (shouldn't happen), anchor it
    for i in 1:Nsites
        if !placed[i]
            pos[i] = (0.0, 0.0)
        end
    end

    return pos
end




# --------------------------
# Simple SVG writer (no external dependencies)
# --------------------------
function write_svg_lattice(path::AbstractString, pos, bonds;
                           show_labels::Bool=true, show_bond_types::Bool=true)

    xs = [p[1] for p in pos]
    ys = [p[2] for p in pos]

    xmin, xmax = minimum(xs), maximum(xs)
    ymin, ymax = minimum(ys), maximum(ys)

    # --------------------
    # Canvas size scales with lattice size
    # --------------------
    margin = 60.0

    spanx = xmax - xmin
    spany = ymax - ymin

    # pixels per unit length (tune once, works for all Nx,Ny)
    ppu = 80.0

    W = spanx * ppu + 2margin
    H = spany * ppu + 2margin

    function X(x)
        return margin + (x - xmin) * ppu
    end
    function Y(y)
        # flip y for SVG
        return H - (margin + (y - ymin) * ppu)
    end

    open(path, "w") do io
        println(io, """<?xml version="1.0" encoding="UTF-8"?>""")
        @printf(io, """<svg xmlns="http://www.w3.org/2000/svg" width="%.0f" height="%.0f" viewBox="0 0 %.0f %.0f">\n""",
                W, H, W, H)
        println(io, """<rect x="0" y="0" width="100%" height="100%" fill="white"/>""")

        # bond style per type (still black, but different stroke styles)
        function bond_style(bt::String)
            if bt == "xx"
                return """stroke="blue" stroke-width="4" stroke-linecap="round" """
            elseif bt == "yy"
                return """stroke="red" stroke-width="4" stroke-linecap="round" """
            else # "zz"
                return """stroke="green" stroke-width="4" stroke-linecap="round" """
            end
        end


        # draw bonds
        for b in bonds
            x1, y1 = pos[b.site1]
            x2, y2 = pos[b.site2]
            @printf(io, """<line x1="%.3f" y1="%.3f" x2="%.3f" y2="%.3f" %s/>\n""",
                    X(x1), Y(y1), X(x2), Y(y2), bond_style(b.bond_type))
        end

        # draw sites
        r = 14.0
        for i in eachindex(pos)
            x, y = pos[i]
            @printf(io, """<circle cx="%.3f" cy="%.3f" r="%.3f" fill="white" stroke="black" stroke-width="3"/>\n""",
                    X(x), Y(y), r)
        end

        # site labels (optional)
        if show_labels
            for i in eachindex(pos)
                x, y = pos[i]
                @printf(io,
                    """<text x="%.3f" y="%.3f"
                        font-family="monospace"
                        font-size="15"
                        font-weight="bold"
                        text-anchor="middle"
                        dominant-baseline="middle"
                        fill="black">%d</text>\n""",
                    X(x), Y(y), i)
            end
        end
                # --------------------
        # Legend
        # --------------------
        legend_x = 60
        legend_y = 80
        dy = 35

        legend = [
            ("x bond", "blue"),
            ("y bond", "red"),
            ("z bond", "green"),
        ]

        for (i, (label, color)) in enumerate(legend)
            y = legend_y + (i - 1) * dy

            # line sample
            @printf(io,
                """<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"
                    stroke="%s" stroke-width="4" stroke-linecap="round"/>""",
                legend_x, y, legend_x + 40, y, color)

            # text
            @printf(io,
                """<text x="%.1f" y="%.1f" font-family="monospace"
                    font-size="18" dominant-baseline="middle">%s</text>\n""",
                legend_x + 55, y, label)
        end

        println(io, "</svg>")
    end
end

# --------------------------
# Optional PDF output via CairoMakie if installed
# (If not installed, it will just skip PDF.)
# --------------------------
function try_write_pdf_lattice(path::AbstractString, pos, bonds)
    try
        @eval begin
            using CairoMakie
        end
    catch
        println("CairoMakie not available -> skipping PDF output (SVG still produced).")
        return
    end

    # compute bounds
    xs = [p[1] for p in pos]
    ys = [p[2] for p in pos]
    xmin, xmax = minimum(xs), maximum(xs)
    ymin, ymax = minimum(ys), maximum(ys)

    fig = CairoMakie.Figure(size = (900, 700))
    ax = CairoMakie.Axis(fig[1,1], aspect = CairoMakie.DataAspect())
    CairoMakie.hidedecorations!(ax)
    CairoMakie.hidespines!(ax)

    # bonds
    for b in bonds
        x1, y1 = pos[b.site1]
        x2, y2 = pos[b.site2]
        CairoMakie.lines!(ax, [x1, x2], [y1, y2])
        mx, my = (x1+x2)/2, (y1+y2)/2
        CairoMakie.text!(ax, b.bond_type, position=(mx, my), align=(:center, :bottom))
    end

    # sites + labels
    for i in eachindex(pos)
        x, y = pos[i]
        CairoMakie.scatter!(ax, [x], [y], markersize=16)
        CairoMakie.text!(ax, string(i), position=(x, y), align=(:center, :top))
    end

    CairoMakie.xlims!(ax, xmin - 0.5, xmax + 0.5)
    CairoMakie.ylims!(ax, ymin - 0.5, ymax + 0.5)

    CairoMakie.save(path, fig)
end

# --------------------------
# Public entry point: draw from config and write into outdir
# --------------------------
function draw_lattice_from_config(config_path::AbstractString, outdir::AbstractString;
                                  svg_name::AbstractString="lattice.svg",
                                  pdf_name::AbstractString="lattice.pdf")
    cfg = TOML.parsefile(config_path)

    Nx = cfg["system"]["Nx"]
    Ny = cfg["system"]["Ny"]
    yperiodic = get(cfg["lattice"], "yperiodic", false)


    Nsites, bonds = honeycomb_bonds(Nx, Ny; yperiodic=yperiodic)
    pos = honeycomb_positions(Nx, Ny)

    isdir(outdir) || mkpath(outdir)

    svg_name = "lattice$(Nx)x$(Ny).svg"
    pdf_name = "lattice$(Nx)x$(Ny).pdf"

    svg_path = joinpath(outdir, svg_name)
    pdf_path = joinpath(outdir, pdf_name)


    println("Drawing honeycomb lattice Nx=$Nx Ny=$Ny (Nsites=$Nsites) -> $svg_path")
    write_svg_lattice(svg_path, pos, bonds; show_labels=true, show_bond_types=true)

    # PDF is optional
    try_write_pdf_lattice(pdf_path, pos, bonds)
end

# --------------------------
# CLI usage (optional)
# julia draw_lattice.jl config.toml outdir
# --------------------------
# --------------------------
# CLI usage:
# julia lattice_plot.jl --outdir <DIR> [--config <config.toml>]
# --------------------------
function _main()
    outdir = nothing
    config_path = "config.toml"

    i = 1
    while i <= length(ARGS)
        arg = ARGS[i]
        if arg == "--outdir"
            i += 1
            i <= length(ARGS) || error("Missing value after --outdir")
            outdir = ARGS[i]
        elseif arg == "--config"
            i += 1
            i <= length(ARGS) || error("Missing value after --config")
            config_path = ARGS[i]
        else
            error("Unknown argument: $arg")
        end
        i += 1
    end

    outdir === nothing && error("Usage: julia lattice_plot.jl --outdir <DIR> [--config <config.toml>]")

    draw_lattice_from_config(config_path, outdir)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    _main()
end