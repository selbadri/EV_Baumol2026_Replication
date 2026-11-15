#==========================================================================================
 An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
 ------------------------------------------------------------------------------------------
 File:    run_full_calibration.jl
 Purpose: Full calibration. Iterate the service-kappa/research-productivity block and the
          demand block until the calibrated parameters and the capital-price deflator ln D
          stop moving, then export the final report and fitted paths.
 ------------------------------------------------------------------------------------------
 Usage:   julia --project=. code/run_full_calibration.jl \
              ["k_m,k_rs_start,k_ns_start"] ["eta_m,eta_rs,sigma"] [max_sweeps] [tolerance]
          All arguments are optional; the defaults start from the published calibration.
 ------------------------------------------------------------------------------------------
 Inputs:  calibrate_service_kappas.jl, calibrate_demand.jl and dump_fit_paths.jl, each run in a
          child process; data files listed in paths.jl
 Outputs: output/calibration/sweep_XX_kappa_zeta.toml (+ _evaluations/)   kappa and zeta per sweep
          output/calibration/sweep_XX_demand.toml                          demand per sweep
          output/calibration/calibration_history.csv                       parameter history
          output/calibration/calibration_results.toml                      final parameters, status
          output/calibration/final_calibration.txt, fit_paths.csv          via dump_fit_paths.jl
 ------------------------------------------------------------------------------------------
 Notes:
 Automated outer fixed point between the kappa-zeta and demand calibration
 drivers. kappa_m is calibrated with the goods research productivity against the
 goods LP trend and patent-stock paths; kappa_rs and kappa_ns are calibrated to the
 slopes of the two service LP trends around the fully converged zeta loop.

 Capital prices enter the model in final-good units, through D(t), the model's chained
 VA deflator in final-good units (PK_DEFL in smm_ge_unified.jl). D is an equilibrium
 object: the demand block returns it at each sweep, the next sweep uses it, and its
 change enters the convergence check like a parameter.

 Each block runs in a fresh Julia process because smm_ge_unified.jl defines
 configuration-dependent constants at load time. Base.julia_cmd() reuses the
 exact Julia executable running this script, so `julia` need not be on PATH.
==========================================================================================#

import Pkg
import TOML
using Printf

include(joinpath(@__DIR__, "paths.jl"))
const ROOT = CAL_PATHS.root
Pkg.activate(ROOT)

parse_triplet(text::String, label::String) = begin
    values = parse.(Float64, split(text, ','))
    length(values) == 3 || error("$label must contain exactly three comma-separated numbers")
    all(isfinite, values) || error("$label must contain finite values")
    values
end

triplet_text(values) = join((@sprintf("%.17g", x) for x in values), ',')
ordered(section) = [Float64(section["m"]), Float64(section["rs"]), Float64(section["ns"])]

const INITIAL_KAPPA = parse_triplet(length(ARGS) >= 1 ? ARGS[1] : "2.537006999597075,2.1558318293414533,4.299089047389644", "initial kappa")
const INITIAL_DEMAND = parse_triplet(length(ARGS) >= 2 ? ARGS[2] : "0.7715027684980625,1.3597778226071489,0.005", "initial demand")
const MAX_SWEEPS = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 10
const TOLERANCE = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 1e-5
MAX_SWEEPS >= 1 || error("max_sweeps must be at least one")
TOLERANCE > 0 || error("tolerance must be positive")

const OUTPUT_DIR = CAL_PATHS.out_calibration
mkpath(OUTPUT_DIR)
const HISTORY_CSV = joinpath(OUTPUT_DIR, "calibration_history.csv")
# Capital-price deflator ln D at 1950, any INPUT_KNOTS and 2023. Without knots it starts at the
# published fixed point, with knots at zero, unless PK_DEFL is set.
const PK_NODES = let s = strip(get(ENV, "INPUT_KNOTS", ""))
    vcat(1950, isempty(s) ? Int[] : sort(parse.(Int, strip.(split(s, ",")))), 2023)
end
haskey(ENV, "PK_DEFL") || (ENV["PK_DEFL"] = length(PK_NODES) == 2 ? "0,0.06380282626824207" :
                                            join(zeros(length(PK_NODES)), ","))
const SUMMARY_TOML = joinpath(OUTPUT_DIR, "calibration_results.toml")

function run_block(script_name::String, args::Vector{String}, result_path::String)
    isfile(result_path) && rm(result_path)
    script = joinpath(@__DIR__, script_name)
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$ROOT $script $args`
    run(addenv(cmd, "CALIBRATION_RESULT_PATH" => result_path))
    isfile(result_path) || error("$script_name finished without writing $result_path")
    TOML.parsefile(result_path)
end

function save_summary(; converged, completed_sweeps, max_change, kappa, demand, zeta0, gzeta,
                      g2zeta, kappa_status, zeta_fit, demand_fit, knots=Float64[], pk_defl=Float64[])
    result = Dict(
        "status" => Dict(
            "converged" => converged,
            "completed_sweeps" => completed_sweeps,
            "maximum_parameter_change" => max_change,
            "tolerance" => TOLERANCE,
        ),
        "kappa" => Dict("m" => kappa[1], "rs" => kappa[2], "ns" => kappa[3]),
        "kappa_status" => kappa_status,
        "demand" => Dict("eta_m" => demand[1], "eta_rs" => demand[2], "sigma" => demand[3]),
        "zeta0" => Dict("m" => zeta0[1], "rs" => zeta0[2], "ns" => zeta0[3]),
        "gzeta" => Dict("m" => gzeta[1], "rs" => gzeta[2], "ns" => gzeta[3]),
        "g2zeta" => Dict("rs" => g2zeta[1], "ns" => g2zeta[2]),
        "zeta_knots" => isempty(knots) ? Dict{String,Any}() :
            Dict("ln_rs" => knots[1:end÷2], "ln_ns" => knots[end÷2+1:end]),
        "zeta_fit" => zeta_fit,
        "demand_fit" => demand_fit,
        "pk_defl" => Dict("nodes" => PK_NODES, "ln_defl" => pk_defl),
        "machine_markups" => get(ENV, "MUX_MODE", "quad"),
    )
    open(SUMMARY_TOML, "w") do io
        TOML.print(io, result)
    end
end

function main()
    kappa = copy(INITIAL_KAPPA)
    demand = copy(INITIAL_DEMAND)
    previous_zeta = nothing
    zeta0 = fill(NaN, 3)
    gzeta = fill(NaN, 3)
    g2zeta = parse.(Float64, split(get(ENV, "GZETA2", "0,0"), ","))
    previous_g2zeta = nothing
    knots = Float64[]
    previous_knots = Float64[]
    # Optional warm start for the first zeta block, "z0_m,z0_rs,z0_ns;g_m,g_rs,g_ns"; later sweeps
    # start from the previous sweep. It does not enter the convergence check.
    # Default: the published calibration (the services knots start from smm_ge_unified.jl's
    # ZETA_KNOTS_BASELINE).
    zw = get(ENV, "ZETA_WARM", "0.32079803472140844,0.46441626451242646,0.15450408408081595;0.008620983183703385,0.0009462573575271577,-0.014334022071997441")
    zeta_warm = isempty(zw) ? nothing :
        vcat(parse.(Float64, split(split(zw, ";")[1], ",")), parse.(Float64, split(split(zw, ";")[2], ",")))
    zeta_warm === nothing || length(zeta_warm) == 6 || error("ZETA_WARM must hold two triples")
    max_change = Inf
    converged = false
    pk_defl = parse.(Float64, split(ENV["PK_DEFL"], ","))
    zeta_fit = Dict{String,Any}()
    kappa_status = Dict{String,Any}()
    demand_fit = Dict{String,Any}()

    open(HISTORY_CSV, "w") do io
        println(io, "sweep,kappa_m,kappa_rs,kappa_ns,eta_m,eta_rs,sigma,zeta0_m,zeta0_rs,zeta0_ns,gzeta_m,gzeta_rs,gzeta_ns,max_change,goods_lp_rmse,goods_patent_rmse,scientist_rmse,ri_lp_residual,nri_lp_residual,demand_rmse,g2zeta_rs,g2zeta_ns")
    end

    println("="^88)
    println("Automated demand-innovation calibration")
    println("  initial kappa  = ", triplet_text(kappa), " (kappa_m is a starting value)")
    println("  initial demand = ", triplet_text(demand))
    println("  max sweeps     = ", MAX_SWEEPS)
    println("  tolerance      = ", TOLERANCE)
    println("  output         = ", OUTPUT_DIR)
    println("="^88)

    completed_sweeps = 0
    for sweep in 1:MAX_SWEEPS
        completed_sweeps = sweep
        println("\n", "#"^88)
        println("OUTER SWEEP $sweep: SERVICE-KAPPA + ZETA BLOCK")
        println("#"^88)

        previous_kappa = copy(kappa)
        zeta_file = joinpath(OUTPUT_DIR, @sprintf("sweep_%02d_kappa_zeta.toml", sweep))
        zeta_args = [triplet_text([kappa[1]]), triplet_text(kappa[2:3]),
                     triplet_text(demand)]
        zstart = previous_zeta !== nothing ? previous_zeta : zeta_warm
        if zstart !== nothing
            append!(zeta_args, [triplet_text(zstart[1:3]), triplet_text(zstart[4:6])])
        end
        zeta_result = run_block(
            "calibrate_service_kappas.jl",
            zeta_args,
            zeta_file,
        )
        kappa = ordered(zeta_result["kappa"])
        kappa_status = zeta_result["status"]
        zeta0 = ordered(zeta_result["zeta0"])
        gzeta = ordered(zeta_result["gzeta"])
        g2zeta = haskey(zeta_result, "g2zeta") ? [Float64(zeta_result["g2zeta"]["rs"]), Float64(zeta_result["g2zeta"]["ns"])] : [0.0, 0.0]
        # Child processes (demand block, next sweep, fit export) solve at this curvature.
        ENV["GZETA2"] = join(string.(g2zeta), ",")
        kn = get(zeta_result, "zeta_knots", Dict{String,Any}())
        knots = haskey(kn, "ln_rs") ? vcat(Float64.(kn["ln_rs"]), Float64.(kn["ln_ns"])) : Float64[]
        isempty(knots) || (ENV["ZETA_KNOTS"] = triplet_text(knots[1:end÷2]) * ";" * triplet_text(knots[end÷2+1:end]))
        zeta_fit = zeta_result["fit"]

        println("\n", "#"^88)
        println("OUTER SWEEP $sweep: DEMAND BLOCK")
        println("#"^88)

        demand_file = joinpath(OUTPUT_DIR, @sprintf("sweep_%02d_demand.toml", sweep))
        demand_result = run_block(
            "calibrate_demand.jl",
            [triplet_text(kappa), triplet_text(zeta0), triplet_text(gzeta),
             triplet_text(demand)],
            demand_file,
        )
        new_demand_section = demand_result["demand"]
        new_demand = [Float64(new_demand_section["eta_m"]),
                      Float64(new_demand_section["eta_rs"]),
                      Float64(new_demand_section["sigma"])]
        demand_fit = demand_result["fit"]

        # Capital-price deflator for the next sweep.
        new_pk = Float64.(demand_result["pk_defl"]["ln_defl"])
        pk_change = maximum(abs.(new_pk .- pk_defl))
        pk_defl = new_pk
        ENV["PK_DEFL"] = join(string.(pk_defl), ",")
        @printf("\nCapital-price deflator ln D at %s: %s (change %.3e)\n", join(PK_NODES, ", "), ENV["PK_DEFL"], pk_change)

        demand_change = maximum(abs.(new_demand .- demand))
        kappa_change = maximum(abs.(kappa .- previous_kappa))
        zeta_change = previous_zeta === nothing ? Inf :
            max(maximum(abs.(vcat(zeta0, gzeta) .- previous_zeta)), maximum(abs.(g2zeta .- previous_g2zeta)),
                isempty(knots) ? 0.0 : maximum(abs.(knots .- previous_knots)))
        max_change = max(demand_change, kappa_change, zeta_change, pk_change)

        @printf("\nOuter sweep %d: demand change = %.3e, kappa change = %.3e, zeta change = %.3e, maximum = %.3e\n",
            sweep, demand_change, kappa_change, zeta_change, max_change)

        open(HISTORY_CSV, "a") do io
            values = vcat([sweep], kappa, new_demand, zeta0, gzeta, [max_change,
                Float64(zeta_fit["goods_lp_rmse"]),
                Float64(zeta_fit["goods_patent_rmse"]),
                Float64(zeta_fit["scientist_rmse"]),
                Float64(zeta_fit["lp_residual_rs"]),
                Float64(zeta_fit["lp_residual_ns"]),
                Float64(demand_fit["overall_rmse"])], g2zeta)
            println(io, join(values, ','))
        end

        demand = new_demand
        previous_zeta = vcat(zeta0, gzeta)
        previous_g2zeta = copy(g2zeta)
        previous_knots = copy(knots)
        converged = sweep > 1 && max_change < TOLERANCE
        save_summary(; converged, completed_sweeps, max_change, kappa, demand, zeta0, gzeta,
                     g2zeta, kappa_status, zeta_fit, demand_fit, knots, pk_defl)
        converged && break
    end

    println("\n", "#"^88)
    println("FINAL FIT-PATH EXPORT")
    println("#"^88)
    dump_script = joinpath(@__DIR__, "dump_fit_paths.jl")
    dump_cmd = `$(Base.julia_cmd()) --startup-file=no --project=$ROOT $dump_script $(triplet_text(kappa)) $(triplet_text(demand)) $(triplet_text(zeta0)) $(triplet_text(gzeta))`
    run(dump_cmd)

    println("\n", "="^88)
    println(converged ? "CALIBRATION CONVERGED" : "CALIBRATION STOPPED AT MAXIMUM SWEEPS")
    println("  completed sweeps = ", completed_sweeps)
    println("  maximum change   = ", max_change)
    println("  kappa             = ", triplet_text(kappa))
    println("  demand            = ", triplet_text(demand))
    println("  zeta0             = ", triplet_text(zeta0))
    println("  gzeta             = ", triplet_text(gzeta))
    println("  g2zeta (rs, ns)    = ", join(string.(g2zeta), ","))
    if !isempty(knots)
        println("  ln zeta knots rs  = ", triplet_text(knots[1:end÷2]))
        println("  ln zeta knots ns  = ", triplet_text(knots[end÷2+1:end]))
    end
    println("  ln D (", join(PK_NODES, ", "), ") = ", join(string.(pk_defl), ","))
    println("  summary           = ", SUMMARY_TOML)
    println("  history           = ", HISTORY_CSV)
    println("  parameter/moments = ", joinpath(OUTPUT_DIR, "final_calibration.txt"))
    println("  fitted paths      = ", joinpath(OUTPUT_DIR, "fit_paths.csv"))
    println("="^88)

    converged || @warn "Outer calibration did not meet the requested tolerance; inspect the history before treating the parameters as final."
end

main()
