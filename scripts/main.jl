push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using KitaevDMRG

project_root = normpath(joinpath(@__DIR__, ".."))
config_path  = joinpath(project_root, "config.toml")
summary_path = joinpath(project_root, "results", "benchmark.tsv")

KitaevDMRG.run_batch(config_path; summary_path=summary_path)
