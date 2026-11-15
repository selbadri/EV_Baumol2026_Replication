#==========================================================================================
 An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
------------------------------------------------------------------------------------------
 File:    baseline.jl
 Purpose: Load the model at the calibrated parameters. Shared by the calibration figures
          (calibration_plots.jl) and the counterfactual exercises (counterfactuals/).
------------------------------------------------------------------------------------------
 Inputs:  output/calibration/final_calibration.txt     calibrated parameters
          output/calibration/calibration_results.toml  research-productivity knots, the
                                                       capital-price deflator ln D and the
                                                       machine-markup form
          and the data files read by smm_ge_unified.jl (see paths.jl)
 Outputs: none; defines baseline_params, baseline_exog and the solver wrappers
==========================================================================================#

import Pkg
include(joinpath(@__DIR__, "paths.jl"))
const _FINAL_CALIBRATION = joinpath(CAL_PATHS.out_calibration, "final_calibration.txt")
Pkg.activate(CAL_PATHS.root)

"""Read the parameter table in final_calibration.txt and reject missing entries."""
function _read_final_calibration(path::AbstractString)
    isfile(path) || error("Missing final calibration report: $path")
    values = Dict{Symbol,Float64}()
    in_parameters = false
    for line in eachline(path)
        if strip(line) == "PARAMETERS"
            in_parameters = true
            continue
        elseif in_parameters && startswith(strip(line), "DATA AND MODEL MOMENTS")
            break
        elseif in_parameters
            m = match(r"^\s*([A-Za-z][A-Za-z0-9_]*)\s+([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)\s+", line)
            m === nothing || (values[Symbol(m.captures[1])] = parse(Float64, m.captures[2]))
        end
    end
    required = (:kappa_m, :kappa_rs, :kappa_ns, :eta_m, :eta_rs, :eta_ns, :sigma,
                :theta_m, :theta_rs, :theta_ns, :zeta0_m, :zeta0_rs, :zeta0_ns,
                :gzeta_m, :gzeta_rs, :gzeta_ns)
    missing = filter(k -> !haskey(values, k), required)
    isempty(missing) || error("Missing parameters in $path: $(join(string.(missing), ", "))")
    values
end

const _FINAL_VALUES = _read_final_calibration(_FINAL_CALIBRATION)

# Services research productivity: the knots of the calibration (BLOCKB_MODE=knots) when the
# results file has them, otherwise the exponential path. In knots mode zeta0 and gzeta of the
# services in final_calibration.txt are only summaries (first knot, mean growth).
import TOML
const _RESULTS_TOML = joinpath(CAL_PATHS.out_calibration, "calibration_results.toml")
let kn = get(TOML.parsefile(_RESULTS_TOML), "zeta_knots", Dict{String,Any}())
    if haskey(kn, "ln_rs")
        ENV["BLOCKB_MODE"] = "knots"
        ENV["ZETA_KNOTS"] = join(string.(Float64.(kn["ln_rs"])), ",") * ";" * join(string.(Float64.(kn["ln_ns"])), ",")
    else
        ENV["BLOCKB_MODE"] = "legacy"
        ENV["ZETA_KNOTS"] = ""
    end
end
# Capital-price deflator ln D of the calibration (see PK_DEFL in smm_ge_unified.jl), and the
# input knots it was computed at (interior nodes; none in the baseline).
let pk = TOML.parsefile(_RESULTS_TOML)["pk_defl"]
    ENV["PK_DEFL"] = join(string.(Float64.(pk["ln_defl"])), ",")
    ENV["INPUT_KNOTS"] = join(string.(Int.(pk["nodes"][2:end-1])), ",")
end
# Machine-producer markups of the calibration: quadratic trends (baseline) or means (MUX_MODE=mean).
ENV["MUX_MODE"] = get(TOML.parsefile(_RESULTS_TOML), "machine_markups", "quad")

# Force the same specification used to produce final_calibration.txt before loading
# the calibration data and model code.
ENV["ALLOC_MODE"] = "stonegeary"
ENV["PHI_FIXED"] = "1"
# PHI_COMMON from the environment lets the exercises run on a phi robustness calibration.
ENV["PHI_COMMON"] = get(ENV, "PHI_COMMON", "0.50")
ENV["PATENT_SERIES"] = "data_moments"
ENV["N_INIT"] = "one"
ENV["KAPPA_MODE"] = "external"
ENV["KAPPA_EXT"] = join((_FINAL_VALUES[:kappa_m], _FINAL_VALUES[:kappa_rs], _FINAL_VALUES[:kappa_ns]), ",")

include(joinpath(@__DIR__, "smm_ge_unified.jl"))

"""Return a NamedTuple with selected fields replaced."""
with_params(p::NamedTuple; kwargs...) = merge(p, (; kwargs...))
with_exog(e::NamedTuple; kwargs...) = merge(e, (; kwargs...))

const baseline_params = (
    source = _FINAL_CALIBRATION,
    T_data = T_smm, T_ge = T_ge,
    chi = chi_smm, rho = rho_smm, alpha = alpha_fixed,
    # PHI_COMMON, the value the calibration used; phi_m_fixed is a hardcoded 0.5.
    gamma = gamma_fixed, phi = PHI_COMMON,
    kappa_m = _FINAL_VALUES[:kappa_m], kappa_rs = _FINAL_VALUES[:kappa_rs],
    kappa_ns = _FINAL_VALUES[:kappa_ns],
    etam = _FINAL_VALUES[:eta_m], etars = _FINAL_VALUES[:eta_rs],
    etans = _FINAL_VALUES[:eta_ns], sigma = _FINAL_VALUES[:sigma],
    theta_m = _FINAL_VALUES[:theta_m], theta_rs = _FINAL_VALUES[:theta_rs],
    theta_ns = _FINAL_VALUES[:theta_ns],
    zetam_ini = _FINAL_VALUES[:zeta0_m], zetars_ini = _FINAL_VALUES[:zeta0_rs],
    zetans_ini = _FINAL_VALUES[:zeta0_ns],
    gzetam = _FINAL_VALUES[:gzeta_m], gzetars = _FINAL_VALUES[:gzeta_rs],
    gzetans = _FINAL_VALUES[:gzeta_ns],
    dMM = dMM_ge, dMR = dMR_ge, dMN = dMN_ge,
    dRM = dRM_ge, dRR = dRR_ge, dRN = dRN_ge,
    dNM = dNM_ge, dNR = dNR_ge, dNN = dNN_ge,
    lambda_N0 = KPSS_PATENT_SHARES_1950,
    zeta_knots = ZETA_KNOTS_DEFAULT,
)

const baseline_exog = (
    muxm = μxmt_ge, muxrs = μxrst_ge, muxns = μxnst_ge,
    muym = μymt_ge, muyrs = μyrst_ge, muyns = μynst_ge,
    zm = zmt_ge, zrs = zrst_ge, zns = znst_ge,
    hbarm = hbarmt_ge, hbarrs = hbarrst_ge, hbarns = hbarnst_ge,
    Lt = Lt_ge, gl = gl_pop_ge, ls = sci_share_tot_ge,
)

_kappa(p) = (m=p.kappa_m, rs=p.kappa_rs, ns=p.kappa_ns)
_zeta0(p) = (p.zetam_ini, p.zetars_ini, p.zetans_ini)
_gzeta(p) = (p.gzetam, p.gzetars, p.gzetans)
"""Research-productivity paths (Goods, RI-S, NRI-S) over T_ge: knots for the services when the
baseline has them, exponential otherwise."""
zeta_paths(p) = p.zeta_knots === nothing ? nothing : _zeta_knot_paths(_zeta0(p), _gzeta(p), p.zeta_knots)

"""One model solve, with theta either recovered from supplied gross-share targets or imposed as a path."""
function _run_model(p::NamedTuple, e::NamedTuple;
                    mode::Symbol=:calibrate,
                    V_seed=nothing,
                    gross_share_targets=(smt1_target_data, srst1_target_data),
                    tau_p_rs_path=nothing, tau_p_ns_path=nothing,
                    theta_m_fixed=nothing, theta_rs_fixed=nothing,
                    theta_m_path=nothing, theta_rs_path=nothing,
                    N_m_lock=nothing, N_rs_lock=nothing, N_ns_lock=nothing)
    Vkw = V_seed === nothing ? NamedTuple() :
          (; V_init_m=collect(V_seed[1]), V_init_rs=collect(V_seed[2]), V_init_ns=collect(V_seed[3]))
    eq = solve_GE(
        p.T_ge, p.chi, p.rho,
        p.dMM, p.dMR, p.dMN, p.dNM, p.dNR, p.dNN, p.dRM, p.dRR, p.dRN,
        p.alpha, _kappa(p), p.sigma, p.etam, p.etars, p.etans, p.gamma, p.phi,
        p.theta_m, p.theta_rs,
        p.gzetam, p.gzetars, p.gzetans,
        p.zetam_ini, p.zetars_ini, p.zetans_ini,
        e.ls, sci_share_rs_smm, sci_share_ns_smm, prs_pm_smm, pns_pm_smm,
        gross_share_targets[1], gross_share_targets[2],
        1.0, 1.0, 1.0, e.gl,
        e.muxm, e.muxrs, e.muxns, e.muym, e.muyrs, e.muyns,
        e.zm, e.zrs, e.zns, e.hbarm, e.hbarrs, e.hbarns, e.Lt;
        mode=mode, Vkw...,
        v_damp=V_DAMP_SMM, iter_max_v=ITER_MAX_V_SMM, tol_v=TOL_V_SMM,
        T_m=1.0, T_rs_fixed=1.0, T_ns_fixed=1.0,
        alloc=:stonegeary, Lbar_m=0.0, Lbar_rs=0.0, Lbar_ns=0.0,
        lambda_N0=p.lambda_N0,
        tau_p_rs_path=tau_p_rs_path, tau_p_ns_path=tau_p_ns_path,
        theta_m_fixed=theta_m_fixed, theta_rs_fixed=theta_rs_fixed,
        theta_m_path=theta_m_path, theta_rs_path=theta_rs_path,
        N_m_lock=N_m_lock, N_rs_lock=N_rs_lock, N_ns_lock=N_ns_lock,
        φrs=p.phi, φns=p.phi,
        ζ_paths=zeta_paths(p),
    )
    eq["converged"][1] == 1.0 || error("GE solve failed to converge (failure code $(eq["fail_code"][1]))")
    eq
end

# Reproduce the calibrated baseline first. This recovers theta exactly as in the
# calibration and obtains the price-wedge paths held fixed in counterfactuals.
const eq_base = _run_model(baseline_params, baseline_exog)

function _initial_va_shares(eq, e)
    a = baseline_params.alpha
    omega = (1-(1-a)*e.muym[1]*e.muxm[1],
             1-(1-a)*e.muyrs[1]*e.muxrs[1],
             1-(1-a)*e.muyns[1]*e.muxns[1])
    # s_j=p_jY_j/Y is the internally consistent gross-expenditure share. The
    # returned p_j paths are independently rebased to one at t=1 for reporting,
    # so p_j*Y_j must not be reconstructed from those rebased paths.
    nominal = (omega[1]*eq["smt"][1],
               omega[2]*eq["srst"][1],
               omega[3]*eq["snst"][1])
    total = sum(nominal)
    ntuple(j -> nominal[j]/total, 3)
end

const _BASELINE_INITIAL_VA = _initial_va_shares(eq_base, baseline_exog)
const baseline_identified = (
    tau_p_rs = copy(eq_base["τprst"]),
    tau_p_ns = copy(eq_base["τpnst"]),
    initial_va_shares = _BASELINE_INITIAL_VA,
    theta = (eq_base["θm"][1], eq_base["θrs"][1], eq_base["θns"][1]),
)
const baseline_context = (V_seed=(copy(eq_base["Vmt"]), copy(eq_base["Vrst"]), copy(eq_base["Vnst"])),)

"""Convert desired initial VA shares into the gross expenditure-share targets expected by solve_GE."""
function _gross_targets_for_initial_va(e::NamedTuple, va=_BASELINE_INITIAL_VA)
    a = baseline_params.alpha
    omega = (1-(1-a)*e.muym[1]*e.muxm[1],
             1-(1-a)*e.muyrs[1]*e.muxrs[1],
             1-(1-a)*e.muyns[1]*e.muxns[1])
    raw = ntuple(j -> va[j]/omega[j], 3)
    scale = sum(raw)
    (raw[1]/scale, raw[2]/scale)
end

function _theta_for_gross_targets(eq, p::NamedTuple, targets)
    sm, srs = targets
    shares = (sm, srs, 1.0-sm-srs)
    epsilon = 1.0-p.sigma
    if abs(epsilon) <= 1e-12
        return shares
    end
    Y1 = eq["Yt"][1]
    # The reported p_j paths are each rebased to one at t=1. Demand uses the
    # pre-rebasing relative prices, which equal the calibrated t=1 price ratios
    # because counterfactuals retain the baseline production-wedge paths.
    weights = (Y1^((p.etam-1.0)*epsilon),
               prs_pm_smm[1]^epsilon * Y1^((p.etars-1.0)*epsilon),
               pns_pm_smm[1]^epsilon * Y1^((p.etans-1.0)*epsilon))
    raw = ntuple(j -> shares[j]/max(weights[j],1e-300), 3)
    scale = sum(raw)
    ntuple(j -> raw[j]/scale, 3)
end

"""
Re-solve the model while preserving the baseline 1950 value-added shares.

Theta is deliberately omitted from the solve unless an explicit theta path is supplied,
so `solve_GE` recovers theta_m and theta_rs at t=1 for every counterfactual. Explicit paths
are reserved for Exercise 5's intentional demand perturbations.
"""
function run_counterfactual(p::NamedTuple, e::NamedTuple;
                            V_seed=baseline_context.V_seed,
                            label="counterfactual",
                            theta_m_path=nothing, theta_rs_path=nothing,
                            N_m_lock=nothing, N_rs_lock=nothing, N_ns_lock=nothing)
    path_theta = theta_m_path !== nothing || theta_rs_path !== nothing
    (theta_m_path === nothing) == (theta_rs_path === nothing) ||
        error("$label must supply theta_m_path and theta_rs_path together")
    targets = _gross_targets_for_initial_va(e)
    eq = _run_model(p, e; mode=:counterfactual, V_seed=V_seed,
                    gross_share_targets=targets,
                    tau_p_rs_path=baseline_identified.tau_p_rs,
                    tau_p_ns_path=baseline_identified.tau_p_ns,
                    theta_m_path=theta_m_path, theta_rs_path=theta_rs_path,
                    N_m_lock=N_m_lock, N_rs_lock=N_rs_lock, N_ns_lock=N_ns_lock)
    if !path_theta
        got = _initial_va_shares(eq, e)
        gap = maximum(abs.(collect(got) .- collect(_BASELINE_INITIAL_VA)))
        gap <= 1e-8 || error("$label failed to preserve baseline initial VA shares (max gap=$gap)")
    end
    eq
end

# Guard against accidentally loading a stale report or a different model specification.
let reported=(baseline_params.theta_m, baseline_params.theta_rs, baseline_params.theta_ns),
    solved=baseline_identified.theta
    maximum(abs.(collect(reported).-collect(solved))) <= 5e-8 ||
        error("Theta recovered by the baseline does not match final_calibration.txt: reported=$reported solved=$solved")
end

# Largest-remainder rounding for the generated tables: parts are rounded to `digits` decimals so
# that they add exactly to the rounded total (consistency rule of 30 September 2026). The parts
# must add to the total within one unit of the last digit, which covers a total that was itself
# adjusted by this rounding one level up; zero parts stay zero unless the rounding requires otherwise.
function round_to_total(parts::AbstractVector{<:Real}, total::Real; digits::Int=2)
    s = 10.0^digits
    abs(sum(parts) - total) < 1 / s || error("parts sum to $(sum(parts)), total is $total")
    target = round(Int, total * s)
    fl = floor.(Int, parts .* s)
    rem = parts .* s .- fl
    k = target - sum(fl)
    out = copy(fl)
    order = k >= 0 ? sortperm(rem, rev=true) : sortperm(rem)
    for i in 1:abs(k)
        out[order[mod1(i, length(out))]] += sign(k)
    end
    out ./ s
end

# Fixed-decimal string without a negative zero.
fmt_fixed(x::Real, digits::Int) = (r = round(x, digits=digits); abs(r) < 0.5 * 10.0^(-digits) ?
    @sprintf("%.*f", digits, 0.0) : @sprintf("%.*f", digits, r))
