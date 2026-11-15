#==========================================================================================
 An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
 ------------------------------------------------------------------------------------------
 File:    calibrate_service_kappas.jl
 Purpose: Service-kappa block. Solve kappa_rs and kappa_ns so that model RI-S and NRI-S
          labor-productivity growth match the data, re-converging the research-productivity
          block (calibrate_zeta_blocks.jl) at every evaluation.
 ------------------------------------------------------------------------------------------
 Usage:   julia --project=. code/calibrate_service_kappas.jl "<kappa_m>" "<kappa_rs,kappa_ns>" \
              "<eta_m,eta_rs,sigma>" ["<z0_m,z0_rs,z0_ns>" "<g_m,g_rs,g_ns>"]
          Normally called by run_full_calibration.jl.
 ------------------------------------------------------------------------------------------
 Inputs:  results of calibrate_zeta_blocks.jl, run in child processes
 Outputs: TOML file at ENV CALIBRATION_RESULT_PATH (default output/calibration/service_kappa_result.toml)
          <result>_evaluations/evaluation_NNN.toml, .log   one per kappa evaluation
          <result>_evaluations/kappa_evaluations.csv       evaluation history
 ------------------------------------------------------------------------------------------
 Notes:
 Calibrate kappa_rs and kappa_ns around the complete converged zeta loop.

 kappa_m is calibrated inside the goods block of calibrate_zeta_blocks.jl, jointly
 with the goods research productivity, against the goods LP and patent-stock paths.
 For every candidate service-kappa pair this driver launches that script in a fresh
 Julia process, lets its blocks reach their fixed point, and returns the two service
 LP-growth residuals. A numerical 2x2 Jacobian updates (kappa_rs, kappa_ns) until
 both residuals are zero. The kappa_m argument is only a starting value.
==========================================================================================#

import Pkg
import TOML
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "paths.jl"))
const ROOT = CAL_PATHS.root
Pkg.activate(ROOT)

length(ARGS) in (3, 5) || error(
    "usage: julia --project=. code/calibrate_service_kappas.jl \"kappa_m\" " *
    "\"kappa_rs,kappa_ns\" \"eta_m,eta_rs,sigma\" " *
    "[\"zeta0_m,zeta0_rs,zeta0_ns\" \"gzeta_m,gzeta_rs,gzeta_ns\"]"
)

parse_values(text::String, n::Int, label::String) = begin
    values = parse.(Float64, split(text, ','))
    length(values) == n || error("$label must contain exactly $n comma-separated numbers")
    all(isfinite, values) || error("$label must contain finite numbers")
    values
end

value_text(values) = join((@sprintf("%.17g", x) for x in values), ',')
ordered(section) = [Float64(section["m"]), Float64(section["rs"]), Float64(section["ns"])]

const KAPPA_M = Ref(parse(Float64, ARGS[1]))   # latest calibrated kappa_m, used as warm start
const KAPPA_START = parse_values(ARGS[2], 2, "service kappa")
const DEMAND = parse_values(ARGS[3], 3, "demand")
const INITIAL_ZETA = length(ARGS) == 5 ?
    vcat(parse_values(ARGS[4], 3, "zeta0"), parse_values(ARGS[5], 3, "gzeta")) : nothing

const KAPPA_LOWER = parse(Float64, get(ENV, "KAPPA_LOWER", "0.05"))
const KAPPA_UPPER = parse(Float64, get(ENV, "KAPPA_UPPER", "30.0"))   # kappa_ns is 5.94 at the baseline
const KAPPA_TOL = parse(Float64, get(ENV, "KAPPA_TOL", "1e-6"))
const KAPPA_MAXITER = parse(Int, get(ENV, "KAPPA_MAXITER", "10"))
const KAPPA_FD_REL = parse(Float64, get(ENV, "KAPPA_FD_REL", "0.01"))
const RESULT_PATH = get(ENV, "CALIBRATION_RESULT_PATH", joinpath(CAL_PATHS.out_calibration, "service_kappa_result.toml"))
const EVAL_DIR = joinpath(dirname(RESULT_PATH), splitext(basename(RESULT_PATH))[1] * "_evaluations")
const HISTORY_PATH = joinpath(EVAL_DIR, "kappa_evaluations.csv")
# Clear files rather than the directory: a file-sync client can hold a lock on a synced
# folder, and deleting the folder itself then fails with EBUSY.
mkpath(EVAL_DIR)
for f in readdir(EVAL_DIR; join=true)
    isfile(f) && rm(f; force=true)
end

KAPPA_M[] > 0 || error("kappa_m must be positive")
0 < KAPPA_LOWER < KAPPA_UPPER || error("invalid service-kappa bounds")
all((KAPPA_LOWER .<= KAPPA_START) .& (KAPPA_START .<= KAPPA_UPPER)) ||
    error("initial service kappas are outside the bounds")
KAPPA_TOL > 0 || error("KAPPA_TOL must be positive")
KAPPA_MAXITER >= 1 || error("KAPPA_MAXITER must be at least one")
KAPPA_FD_REL > 0 || error("KAPPA_FD_REL must be positive")

const EVAL_COUNT = Ref(0)

open(HISTORY_PATH, "w") do io
    println(io, "evaluation,kappa_m,kappa_rs,kappa_ns,lp_residual_rs,lp_residual_ns,zeta0_m,zeta0_rs,zeta0_ns,gzeta_m,gzeta_rs,gzeta_ns")
end

function evaluate_kappa(kappa_service::Vector{Float64}, zeta_start)
    EVAL_COUNT[] += 1
    id = EVAL_COUNT[]
    kappa = [KAPPA_M[], kappa_service[1], kappa_service[2]]
    result_path = joinpath(EVAL_DIR, @sprintf("evaluation_%03d.toml", id))
    log_path = joinpath(EVAL_DIR, @sprintf("evaluation_%03d.log", id))
    isfile(result_path) && rm(result_path)

    args = [value_text(kappa), value_text(DEMAND)]
    if zeta_start !== nothing
        append!(args, [value_text(zeta_start[1:3]), value_text(zeta_start[4:6])])
    end
    script = joinpath(@__DIR__, "calibrate_zeta_blocks.jl")
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$ROOT $script $args`

    println("  evaluation $id: kappa = ", value_text(kappa))
    open(log_path, "w") do io
        run(pipeline(addenv(cmd, "CALIBRATION_RESULT_PATH" => result_path),
                     stdout=io, stderr=io))
    end
    isfile(result_path) || error("zeta calibration did not write $result_path; inspect $log_path")

    result = TOML.parsefile(result_path)
    kappa[1] = Float64(result["kappa"]["m"])
    KAPPA_M[] = kappa[1]
    zeta0 = ordered(result["zeta0"])
    gzeta = ordered(result["gzeta"])
    g2zeta = haskey(result, "g2zeta") ? [Float64(result["g2zeta"]["rs"]), Float64(result["g2zeta"]["ns"])] : [0.0, 0.0]
    # The next evaluation warm-starts from this curvature (inherited by the child process).
    ENV["GZETA2"] = join(string.(g2zeta), ",")
    # BLOCKB_MODE=knots: the services knot values travel the same way.
    kn = get(result, "zeta_knots", Dict{String,Any}())
    knots = haskey(kn, "ln_rs") ?
        value_text(Float64.(kn["ln_rs"])) * ";" * value_text(Float64.(kn["ln_ns"])) : ""
    isempty(knots) || (ENV["ZETA_KNOTS"] = knots)
    fit = result["fit"]
    residual = [Float64(fit["lp_residual_rs"]), Float64(fit["lp_residual_ns"])]
    all(isfinite, residual) || error("non-finite service LP residual at kappa=$(value_text(kappa))")
    open(HISTORY_PATH, "a") do io
        println(io, join(vcat([id], kappa, residual, zeta0, gzeta), ','))
    end
    @printf("                  residuals = (%+.6e, %+.6e), max = %.3e\n",
            residual[1], residual[2], maximum(abs, residual))
    return (; kappa_service=copy(kappa_service), kappa_m=kappa[1], residual, zeta0, gzeta, g2zeta, knots, fit,
            result_path, log_path)
end

function finite_difference_point(x::Vector{Float64}, j::Int)
    h = KAPPA_FD_REL * max(abs(x[j]), 1.0)
    xp = copy(x)
    xp[j] = min(x[j] + h, KAPPA_UPPER)
    if xp[j] == x[j]
        xp[j] = max(x[j] - h, KAPPA_LOWER)
    end
    xp[j] != x[j] || error("cannot perturb kappa coordinate $j inside bounds")
    xp
end

function save_result(best; converged::Bool, iterations::Int)
    output = Dict(
        "status" => Dict(
            "converged" => converged,
            "iterations" => iterations,
            "evaluations" => EVAL_COUNT[],
            "maximum_lp_residual" => maximum(abs, best.residual),
            "tolerance" => KAPPA_TOL,
        ),
        "kappa" => Dict("m" => best.kappa_m,
                        "rs" => best.kappa_service[1],
                        "ns" => best.kappa_service[2]),
        "zeta0" => Dict("m" => best.zeta0[1], "rs" => best.zeta0[2], "ns" => best.zeta0[3]),
        "gzeta" => Dict("m" => best.gzeta[1], "rs" => best.gzeta[2], "ns" => best.gzeta[3]),
        "g2zeta" => Dict("rs" => best.g2zeta[1], "ns" => best.g2zeta[2]),
        "zeta_knots" => isempty(best.knots) ? Dict{String,Any}() :
            Dict("ln_rs" => parse.(Float64, split(split(best.knots, ";")[1], ",")),
                 "ln_ns" => parse.(Float64, split(split(best.knots, ";")[2], ","))),
        "fit" => best.fit,
    )
    open(RESULT_PATH, "w") do io
        TOML.print(io, output)
    end
end

function main()
    x = clamp.(copy(KAPPA_START), KAPPA_LOWER, KAPPA_UPPER)
    current = evaluate_kappa(x, INITIAL_ZETA)
    best = current
    converged = maximum(abs, current.residual) < KAPPA_TOL
    completed_iterations = 0

    for iteration in 1:KAPPA_MAXITER
        converged && break
        completed_iterations = iteration
        println("\nKappa iteration $iteration")

        p1 = finite_difference_point(x, 1)
        p2 = finite_difference_point(x, 2)
        # Every evaluation of this iteration warm-starts the knots from the current point.
        _warm!() = isempty(current.knots) || (ENV["ZETA_KNOTS"] = current.knots)
        _warm!(); e1 = evaluate_kappa(p1, vcat(current.zeta0, current.gzeta))
        _warm!(); e2 = evaluate_kappa(p2, vcat(current.zeta0, current.gzeta))
        J = hcat((e1.residual .- current.residual) ./ (p1[1] - x[1]),
                 (e2.residual .- current.residual) ./ (p2[2] - x[2]))
        all(isfinite, J) || error("non-finite 2x2 kappa Jacobian")
        s = svdvals(J)
        minimum(s) > 1e-12 || error("service-kappa Jacobian is numerically singular: singular values=$s")
        @printf("  Jacobian = [%+.4e %+.4e; %+.4e %+.4e], condition = %.3e\n",
                J[1,1], J[1,2], J[2,1], J[2,2], maximum(s)/minimum(s))

        step = -(J \ current.residual)
        step_limit = 0.5 .* max.(abs.(x), 0.25)
        step = clamp.(step, .-step_limit, step_limit)

        accepted = nothing
        for scale in (1.0, 0.5, 0.25, 0.125, 0.0625)
            trial_x = clamp.(x .+ scale .* step, KAPPA_LOWER, KAPPA_UPPER)
            maximum(abs.(trial_x .- x)) > 0 || continue
            _warm!(); trial = evaluate_kappa(trial_x, vcat(current.zeta0, current.gzeta))
            if maximum(abs, trial.residual) < maximum(abs, current.residual)
                accepted = trial
                break
            end
        end

        if accepted === nothing
            candidates = (e1, e2)
            candidate_norms = map(e -> maximum(abs, e.residual), candidates)
            idx = argmin(candidate_norms)
            candidate_norms[idx] < maximum(abs, current.residual) ||
                error("2x2 kappa update failed to reduce service LP residuals; inspect $EVAL_DIR")
            accepted = candidates[idx]
        end

        current = accepted
        x = copy(current.kappa_service)
        maximum(abs, current.residual) < maximum(abs, best.residual) && (best = current)
        converged = maximum(abs, current.residual) < KAPPA_TOL
    end

    maximum(abs, current.residual) <= maximum(abs, best.residual) && (best = current)
    converged = maximum(abs, best.residual) < KAPPA_TOL
    isempty(best.knots) || (ENV["ZETA_KNOTS"] = best.knots)
    save_result(best; converged, iterations=completed_iterations)

    println("\n", "="^80)
    println(converged ? "SERVICE KAPPA CALIBRATION CONVERGED" : "SERVICE KAPPA CALIBRATION DID NOT CONVERGE")
    @printf("  kappa = (%.8f, %.8f, %.8f)\n", best.kappa_m, best.kappa_service[1], best.kappa_service[2])
    @printf("  service LP residuals = (%+.8e, %+.8e)\n", best.residual...)
    println("  result = ", RESULT_PATH)
    println("  evaluation logs = ", EVAL_DIR)
    println("  evaluation history = ", HISTORY_PATH)
    println("="^80)

    converged || error("service-kappa calibration failed to reach tolerance $KAPPA_TOL")
end

main()
