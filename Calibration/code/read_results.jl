#==========================================================================================
 An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
 ------------------------------------------------------------------------------------------
 File:    read_results.jl
 Purpose: Read a stored calibration without loading the model. Shared by
          calibration_table.jl and robustness.jl.
 ------------------------------------------------------------------------------------------
 Inputs:  <dir>/final_calibration.txt, <dir>/calibration_results.toml
 Outputs: none
==========================================================================================#

import TOML

const _NUM = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"

"""Parameters and data/model moments from final_calibration.txt in `dir`."""
function read_final_calibration(dir::AbstractString)
    path = joinpath(dir, "final_calibration.txt")
    isfile(path) || error("Missing $path; run the calibration first.")
    params = Dict{String,Float64}()
    moments = Dict{String,Tuple{Float64,Float64}}()   # label => (data, model)
    section = :none
    for line in eachline(path)
        t = strip(line)
        t == "PARAMETERS" && (section = :params; continue)
        startswith(t, "DATA AND MODEL MOMENTS") && (section = :moments; continue)
        t == "NOTES" && (section = :none; continue)
        if section == :params
            m = match(Regex("^\\s*([A-Za-z][A-Za-z0-9_]*)\\s+(" * _NUM.pattern * ")\\s+"), line)
            m === nothing || (params[m.captures[1]] = parse(Float64, m.captures[2]))
        elseif section == :moments
            # label, then data, model, model-data, role (label may contain spaces and digits)
            m = match(Regex("^(.*?)\\s+(" * _NUM.pattern * ")\\s+(" * _NUM.pattern * ")\\s+(" *
                            _NUM.pattern * ")\\s+\\S"), line)
            m === nothing || (moments[strip(m.captures[1])] =
                (parse(Float64, m.captures[2]), parse(Float64, m.captures[3])))
        end
    end
    isempty(params) && error("No parameters found in $path")
    (; params, moments)
end

"""calibration_results.toml in `dir`."""
read_results_toml(dir::AbstractString) = TOML.parsefile(joinpath(dir, "calibration_results.toml"))
