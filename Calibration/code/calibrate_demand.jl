#==========================================================================================
 An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
 ------------------------------------------------------------------------------------------
 File:    calibrate_demand.jl
 Purpose: Demand block. Recover (eta_m, eta_rs, sigma), eta_ns = 1, holding kappa and the
          research-productivity paths fixed.
 ------------------------------------------------------------------------------------------
 Usage:   julia --project=. code/calibrate_demand.jl "<k_m,k_rs,k_ns>" "<z0_m,z0_rs,z0_ns>" \
              "<g_m,g_rs,g_ns>" ["<eta_m,eta_rs,sigma>"]
          Normally called by run_full_calibration.jl.
 ------------------------------------------------------------------------------------------
 Inputs:  data files listed in paths.jl (data_moments.xlsx, Cross_Citations CSVs)
 Outputs: TOML file at ENV CALIBRATION_RESULT_PATH (set by run_full_calibration.jl):
          output/calibration/sweep_XX_demand.toml, with the demand parameters, the fit and
          ln D, the chained VA deflator in final-good units, at 1950, any INPUT_KNOTS and 2023
 ------------------------------------------------------------------------------------------
 Notes:
 Demand block: recover (eta_m, eta_rs, sigma) with eta_ns normalised to one,
 holding the research productivity path and kappa fixed.

 Objective, 148 residuals against three parameters:
   goods value-added share path        74 years
   RI-S value-added share path         74 years
 Model VA shares are omega_j s_j / sum_i omega_i s_i (see _va_shares), since the
 expenditure shares s_j are shares of gross revenue and differ by the markup wedges.

 Prices. The data measure sectoral prices relative to the private VA deflator, a Tornqvist
 index, and the model prices are measured relative to a Fisher chain of the three model
 sector prices with value-added weights (on the data, the two differ by 1e-5 log points
 over 1950-2023). The solver imposes the observed relative service prices through the
 production wedge tau^p (equilibrium_path.jl), so the chained goods price follows from the
 value-added shares. It is therefore reported as a check, not targeted. The calibration
 conditions on observed prices and never tests the model's own price predictions.

 The block also returns ln D(t), the gap between the final-good numeraire and the chained
 deflator, at 1950, any INPUT_KNOTS and 2023. The next sweep expresses capital prices in
 final-good units with it (PK_DEFL in smm_ge_unified.jl).

 sigma is bounded below at 0.005 and settles at this floor in the baseline.

 Run this alternately with calibrate_zeta_blocks.jl until both stop moving. The
 blocks are genuinely coupled in both directions: a 20 per cent change in eta_rs
 moves the RI-S researcher share by up to 0.137, and doubling zeta0 moves the 2023
 RI-S value-added share by 0.042. In practice the loop closes on the second sweep,
 with the demand parameters moving about 0.2 per cent.
==========================================================================================#

# Model specification (read by smm_ge_unified.jl when it is loaded).
ENV["ALLOC_MODE"]="stonegeary"; ENV["PHI_FIXED"]="1"; ENV["PHI_COMMON"]=get(ENV, "PHI_COMMON", "0.50")
length(ARGS) in (3, 4) || error("usage: julia --project=. code/calibrate_demand.jl \"k_m,k_rs,k_ns\" \"zeta0_m,zeta0_rs,zeta0_ns\" \"gzeta_m,gzeta_rs,gzeta_ns\" [\"eta_m,eta_rs,sigma\"]")
ENV["KAPPA_MODE"]="external"; ENV["KAPPA_EXT"] = ARGS[1]
import Pkg
import TOML
include(joinpath(@__DIR__, "paths.jl"))
Pkg.activate(CAL_PATHS.root)
include(joinpath(@__DIR__, "smm_ge_unified.jl"))
using Printf, Statistics
const Z0 = Tuple(parse.(Float64, split(ARGS[2],",")))
const GZ = Tuple(parse.(Float64, split(ARGS[3],",")))
const P0 = length(ARGS) == 4 ? parse.(Float64, split(ARGS[4],",")) :
                              [0.72779181, 1.63321256, 0.20144784]
length(P0) == 3 || error("initial demand must contain eta_m,eta_rs,sigma")
const SM_DAT  = _target_path("VASHARE_GOODS")
const SRS_DAT = _target_path("VASHARE_RI")
const LPM_DAT = log.(pm_m_tgt ./ pm_m_tgt[1])     # goods price, cumulative log change
@printf("data check: sm %d obs %.4f->%.4f | srs %.4f->%.4f | pm ratio %.4f (moment8 %.4f)\n",
        length(SM_DAT), SM_DAT[1], SM_DAT[T_smm], SRS_DAT[1], SRS_DAT[T_smm],
        exp(LPM_DAT[T_smm]), data_moments[8])
sol(p) = try _solved_or_nothing(_solve_GE_p(p, Z0, GZ, (0.0,0.0))) catch; nothing end
# Model goods price relative to a Fisher chain of sector prices with VA weights.
function lpm_chain(eq, va)
    P = (eq["pmt"], eq["prst"], eq["pnst"]); W = (va.m, va.rs, va.ns)
    x = zeros(T_smm)
    for t in 2:T_smm
        rel = [P[j][t] / P[j][t-1] for j in 1:3]
        lasp = sum(W[j][t-1] * rel[j] for j in 1:3)
        paas = 1 / sum(W[j][t] / rel[j] for j in 1:3)
        x[t] = x[t-1] + log(rel[1]) - 0.5 * (log(lasp) + log(paas))
    end
    x
end
# Model goods price in final-good units.
lpm_index(eq) = [log(eq["pmt"][t] / eq["pmt"][1]) for t in 1:T_smm]
function dev(p)
    eq = sol(p); eq === nothing && return nothing
    va = _va_shares(eq)
    r = Float64[]
    for t in 1:T_smm; push!(r, va.m[t]  - SM_DAT[t]);  end
    for t in 1:T_smm; push!(r, va.rs[t] - SRS_DAT[t]); end
    all(isfinite, r) ? r : nothing
end
rms(x) = sqrt(mean(abs2, x))
function main()
    lb = [0.05, 0.05, parse(Float64, get(ENV, "SIGMA_LB", "0.005"))]; ub = [5.0, 5.0, 1.5]
    best = (Inf, nothing)
    for s0 in (copy(P0), [0.5,1.2,0.4], [1.0,2.0,0.15], [0.9,1.4,0.8], [0.6,1.8,0.6])
        x, r, ok = _lm_solve(p -> dev(clamp.(p, lb, ub)), copy(s0), lb, ub;
                             maxiter=300, tol=1e-14, budget=240.0)
        v = (r !== nothing && all(isfinite, r)) ? rms(r) : Inf
        v < best[1] && (best = (v, clamp.(x, lb, ub)))
    end
    v, p = best; p === nothing && (println("failed"); return)
    r = dev(p); eq = sol(p); va = _va_shares(eq)
    pr_chain = rms(lpm_chain(eq, va) .- LPM_DAT)
    # ln D, the chained VA deflator in final-good units, for the next sweep.
    lnD = lpm_index(eq) .- lpm_chain(eq, va)
    lnD_nodes = [lnD[y - 1949] for y in _INPUT_NODES]
    println("  ln D at ", _INPUT_NODES, ": ", join(round.(lnD_nodes, digits=6), ","),
            "   (input PK_DEFL ", join(PK_DEFL, ","), ")")
    @printf("\n  eta_m = %.6f   eta_rs = %.6f   sigma = %.6f   (eta_ns = 1 normalised)\n", p...)
    @printf("  was    %.6f            %.6f           %.6f\n", P0...)
    bnd = [abs(p[i]-lb[i])<1e-6 || abs(p[i]-ub[i])<1e-6 for i in 1:3]
    println("  at a bound: ", any(bnd) ? "YES $bnd" : "no")
    @printf("\n  overall RMSE %.6f  |  goods VA %.6f   RI-S VA %.6f   NRI-S VA %.6f\n",
            v, rms(r[1:T_smm]), rms(r[T_smm+1:2T_smm]), rms(va.ns[1:T_smm] .- (1 .- SM_DAT .- SRS_DAT)))
    @printf("  goods price RMSE (chained VA deflator, check) %.6f\n", pr_chain)
        result_path = get(ENV, "CALIBRATION_RESULT_PATH", "")
        if !isempty(result_path)
                open(result_path, "w") do io
                        TOML.print(io, Dict(
                                "demand" => Dict("eta_m" => p[1], "eta_rs" => p[2], "sigma" => p[3]),
                                "fit" => Dict("overall_rmse" => v,
                                                          "goods_va_rmse" => rms(r[1:T_smm]),
                                                          "ri_services_va_rmse" => rms(r[T_smm+1:2T_smm]),
                                                          "goods_price_rmse" => pr_chain),
                                "at_bound" => any(bnd),
                                "pk_defl" => Dict("nodes" => _INPUT_NODES, "ln_defl" => lnD_nodes),
                        ))
                end
                println("  machine-readable result: ", result_path)
        end
    println("\n  endpoint fit")
    @printf("    goods VA share 1950  %.4f (%.4f)   2023  %.4f (%.4f)\n",
            va.m[1], SM_DAT[1], va.m[T_smm], SM_DAT[T_smm])
    @printf("    RI-S  VA share 1950  %.4f (%.4f)   2023  %.4f (%.4f)\n",
            va.rs[1], SRS_DAT[1], va.rs[T_smm], SRS_DAT[T_smm])
    @printf("    goods price 2023/1950, chained  %.4f (%.4f)\n",
            exp(lpm_chain(eq, va)[T_smm]), data_moments[8])
    println("\n  knock-on at the OLD zeta/kappa (shows whether the outer loop is needed)")
    sci(k) = [eq[k][t]/max(eq["Lst"][t],1e-12) for t in 1:T_smm]
    @printf("    scientist path RMSE %.6f   (was 0.039769 at the old demand parameters)\n",
            sqrt(mean(abs2, vcat(sci("Lsrst").-sci_share_rs_smm, sci("Lsnst").-sci_share_ns_smm))))
    for (nm,Y,Lp,k) in (("goods","Ymt","Lpmt",1),("RI-S","Yrst","Lprst",2),("NRI-S","Ynst","Lpnst",3))
        mo = _ct_lp_ep(eq[Y], eq[Lp])
        @printf("    LP %-6s model %.6f  data %.6f  %+.2f%%\n", nm, mo, _LP_TARGETS[k],
                100*(mo-_LP_TARGETS[k])/_LP_TARGETS[k])
    end
end
main()
