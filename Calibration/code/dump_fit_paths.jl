#==========================================================================================
 An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
 ------------------------------------------------------------------------------------------
 File:    dump_fit_paths.jl
 Purpose: Solve the model once at a given parameter set; write the fitted annual paths and the
          final calibration report.
 ------------------------------------------------------------------------------------------
 Usage:   julia --project=. code/dump_fit_paths.jl "<k_m,k_rs,k_ns>" "<eta_m,eta_rs,sigma>" \
              "<z0_m,z0_rs,z0_ns>" "<g_m,g_rs,g_ns>"
          Normally called by run_full_calibration.jl after the last sweep.
 ------------------------------------------------------------------------------------------
 Inputs:  data files listed in paths.jl (data_moments.xlsx, Cross_Citations CSVs)
 Outputs: output/calibration/fit_paths.csv           model (m_), data (d_) and LP trend (t_)
                                                     paths, 1950-2023
          output/calibration/final_calibration.txt   parameters and data/model moments
 ------------------------------------------------------------------------------------------
 Notes:
 Solve the model once at a given parameter set and write every model and data
 series used in the fit diagnostics to a CSV, one row per year 1950 to 2023.

 Levels are not comparable between model and data for prices, productivity and idea
 stocks, so those columns are cumulative log change since 1950. Shares are levels.
 Columns are prefixed m_ for model, d_ for data and t_ for the LP trends the calibration
 targets.
==========================================================================================#

# Model specification (read by smm_ge_unified.jl when it is loaded).
ENV["ALLOC_MODE"]="stonegeary"; ENV["PHI_FIXED"]="1"; ENV["PHI_COMMON"]=get(ENV, "PHI_COMMON", "0.50")
length(ARGS) == 4 || error("usage: julia --project=. code/dump_fit_paths.jl \"k_m,k_rs,k_ns\" \"eta_m,eta_rs,sigma\" \"zeta0_m,zeta0_rs,zeta0_ns\" \"gzeta_m,gzeta_rs,gzeta_ns\"")
ENV["KAPPA_MODE"]="external"; ENV["KAPPA_EXT"]=ARGS[1]
import Pkg
include(joinpath(@__DIR__, "paths.jl"))
Pkg.activate(CAL_PATHS.root)
include(joinpath(@__DIR__, "smm_ge_unified.jl"))
using Printf
const P  = parse.(Float64, split(ARGS[2], ","))
const Z0 = Tuple(parse.(Float64, split(ARGS[3], ",")))
const GZ = Tuple(parse.(Float64, split(ARGS[4], ",")))
eq = _solve_GE_p(P, Z0, GZ, (0.0,0.0))
_solved(eq) || error(@sprintf("solve_GE did not converge (fail_t %d, last V change %.3e); no fit paths written",
                              Int(eq["fail_t"][1]), eq["diag_diff_V"][1]))
VA = _va_shares(eq)
const NM_DAT = let i = [findfirst(==(y), VARIETY.years) for y in 1950:2023]
    log.(VARIETY.levels.m[i] ./ VARIETY.levels.m[i[1]])
end
yrs = collect(1950:2023)
idx(v) = [log(v[t]/v[1]) for t in 1:T_smm]
lp(Y,L) = eq[Y] ./ eq[L]
SM = _target_path("VASHARE_GOODS"); SRS = _target_path("VASHARE_RI")
dLPm  = [_getval(y,"LP_GOODS") for y in 1950:2023]
dLPrs = [_getval(y,"LP_RI")    for y in 1950:2023]
dLPns = [_getval(y,"LP_NRI")   for y in 1950:2023]
Ls = eq["Lst"]
cols = Dict{String,Vector{Float64}}(
 "year"=>Float64.(yrs),
 "m_va_goods"=>VA.m,           "d_va_goods"=>SM,
 "m_va_ris"=>VA.rs,            "d_va_ris"=>SRS,
 "m_va_nris"=>VA.ns,           "d_va_nris"=>(1 .- SM .- SRS),
 "m_price_goods"=>idx(eq["pmt"]),            "d_price_goods"=>idx(pm_m_tgt),
 "m_lp_goods"=>idx(lp("Ymt","Lpmt")),        "d_lp_goods"=>idx(dLPm),
 "m_lp_ris"=>idx(lp("Yrst","Lprst")),        "d_lp_ris"=>idx(dLPrs),
 "m_lp_nris"=>idx(lp("Ynst","Lpnst")),       "d_lp_nris"=>idx(dLPns),
 "m_sci_ris"=>[eq["Lsrst"][t]/max(Ls[t],1e-12) for t in 1:T_smm],  "d_sci_ris"=>_interp_fill("RDSH_RI"),
 "m_sci_nris"=>[eq["Lsnst"][t]/max(Ls[t],1e-12) for t in 1:T_smm], "d_sci_nris"=>_interp_fill("RDSH_NRI"),
 "m_n_goods"=>idx(eq["Nmt"]), "d_n_goods"=>NM_DAT, "m_n_ris"=>idx(eq["Nrst"]), "m_n_nris"=>idx(eq["Nnst"]),
 "t_lp_goods"=>_LP_TREND.m, "t_lp_ris"=>_LP_TREND.rs, "t_lp_nris"=>_LP_TREND.ns)
ks = ["year","m_va_goods","d_va_goods","m_va_ris","d_va_ris","m_va_nris","d_va_nris",
      "m_price_goods","d_price_goods","m_lp_goods","d_lp_goods","m_lp_ris","d_lp_ris",
      "m_lp_nris","d_lp_nris","m_sci_ris","d_sci_ris","m_sci_nris","d_sci_nris",
      "m_n_goods","d_n_goods","m_n_ris","m_n_nris","t_lp_goods","t_lp_ris","t_lp_nris"]
mkpath(CAL_PATHS.out_calibration)
const OUT_PATH = joinpath(CAL_PATHS.out_calibration, "fit_paths.csv")
open(OUT_PATH,"w") do io
    println(io, join(ks, ","))
    for t in 1:T_smm; println(io, join([@sprintf("%.6f", cols[k][t]) for k in ks], ",")); end
end
println("wrote ", OUT_PATH, "  rows=", T_smm)

# Human-readable final calibration report. The rows below distinguish genuine
# identifying moments/path summaries from normalizations and untargeted prices.
const REPORT_DIR = CAL_PATHS.out_calibration
mkpath(REPORT_DIR)
const REPORT_PATH = joinpath(REPORT_DIR, "final_calibration.txt")
# Block B fits the smoothed annual researcher path, so the table reports that path;
# curv_ends matches the census observations at 1950 and 2023.
const _CURV_ENDS = get(ENV, "BLOCKB_MODE", "knots") in ("curv_ends", "ends")
const _SCI_ROLE = _CURV_ENDS ? "target: census, matched" : "target: annual path"
_sci_tgt(col, t) = _CURV_ENDS ? _getval(1949 + t, col) :
    (col == "RDSH_RI" ? sci_share_rs_smm[t] : sci_share_ns_smm[t])

ct_growth(v) = log(v[T_smm] / v[1]) / (T_smm - 1)
# Sector prices relative to a Fisher chain of the three sector prices with VA weights, the
# counterpart of the private VA deflator that deflates the data (as in calibrate_demand.jl).
const _CHAIN_LN = let P = (eq["pmt"], eq["prst"], eq["pnst"]), W = (VA.m, VA.rs, VA.ns)
    x = zeros(T_smm)
    for t in 2:T_smm
        rel = [P[j][t] / P[j][t-1] for j in 1:3]
        lasp = sum(W[j][t-1] * rel[j] for j in 1:3)
        paas = 1 / sum(W[j][t] / rel[j] for j in 1:3)
        x[t] = x[t-1] + 0.5 * (log(lasp) + log(paas))
    end
    x
end
chained_ratio(v) = (v[T_smm] / v[1]) / exp(_CHAIN_LN[T_smm])
ratio_2023_1950(v) = v[T_smm] / v[1]
model_sci_rs(t) = eq["Lsrst"][t] / max(eq["Lst"][t], 1e-12)
model_sci_ns(t) = eq["Lsnst"][t] / max(eq["Lst"][t], 1e-12)

# LP targets: the slope through 1950 of the cumulative log trend (_ct_lp_ep), compared with the
# same statistic of the model path; endpoint growth of the raw observations is reported beside it.
_lp_rows(nm, k, Y, L, role) =
    [("LP trend slope, $nm, 1950-2023", _LP_TARGETS[k], _ct_lp_ep(eq[Y], eq[L]), role),
     ("LP growth, $nm, 1950-2023 (raw endpoints)", data_moments[k], ct_growth(lp(Y, L)), "untargeted")]
moment_rows = [
        _lp_rows("Goods", 1, "Ymt", "Lpmt", "check: annual trend path targeted")...,
        ("Goods patent stock, ln N(2023)/N(1950)", NM_DAT[end], log(eq["Nmt"][T_smm]/eq["Nmt"][1]), "target: annual path"),
        _lp_rows("RI-S", 2, "Yrst", "Lprst", "target: trend slope")...,
        _lp_rows("NRI-S", 3, "Ynst", "Lpnst", "target: trend slope")...,
        ("Scientist share, RI-S, 1950", _sci_tgt("RDSH_RI", 1), model_sci_rs(1), _SCI_ROLE),
        ("Scientist share, RI-S, 2023", _sci_tgt("RDSH_RI", T_smm), model_sci_rs(T_smm), _SCI_ROLE),
        ("Scientist share, NRI-S, 1950", _sci_tgt("RDSH_NRI", 1), model_sci_ns(1), _SCI_ROLE),
        ("Scientist share, NRI-S, 2023", _sci_tgt("RDSH_NRI", T_smm), model_sci_ns(T_smm), _SCI_ROLE),
        ("Price p_m(2023)/p_m(1950), chained VA deflator", ratio_2023_1950(pm_m_tgt), chained_ratio(eq["pmt"]),
         "check: follows from the VA shares"),
        ("Price p_m(2023)/p_m(1950), final-good units", ratio_2023_1950(pm_m_tgt), ratio_2023_1950(eq["pmt"]),
         "reported: numeraire differs from the data"),
        ("Price p_rs(2023)/p_rs(1950), chained VA deflator", ratio_2023_1950(pm_rs_tgt), chained_ratio(eq["prst"]),
         "relative to p_m imposed (tau^p)"),
        ("Price p_ns(2023)/p_ns(1950), chained VA deflator", ratio_2023_1950(pm_ns_tgt), chained_ratio(eq["pnst"]),
         "relative to p_m imposed (tau^p)"),
        ("Capital-price deflator ln D(2023)", NaN, log(_PK_DEFL_PATH[T_smm]), "equilibrium object"),
        ("VA share, Goods, 1950 (raw)",       va1950_target_data[1], VA.m[1],                  "anchor for theta recovery"),
        ("VA share, RI-S, 1950 (raw)",        va1950_target_data[2], VA.rs[1],                 "anchor for theta recovery"),
        ("VA share, Goods, 2023",             SM[end],        VA.m[T_smm],                     "target: annual path"),
        ("VA share, RI-S, 2023",              SRS[end],       VA.rs[T_smm],                    "target: annual path"),
        ("VA share, NRI-S, 2023",             1 - SM[end] - SRS[end], VA.ns[T_smm],            "implied by adding up"),
        ("Goods LP index, 1950",              1.0,            1.0,                              "normalization, not a moment"),
]

open(REPORT_PATH, "w") do io
        println(io, "FINAL CALIBRATION: PARAMETERS AND MOMENTS")
        println(io, "="^96)
        println(io, "Specification: homogeneous researchers, free mobility, zero research floors")
        println(io, "Sample: 1950-2023")
        println(io)
        println(io, "PARAMETERS")
        println(io, "-"^96)
        @printf(io, "%-28s %18s  %s\n", "parameter", "value", "status")
        @printf(io, "%-28s %18.10f  %s\n", "kappa_m",  kappa_fixed.m, "calibrated")
        @printf(io, "%-28s %18.10f  %s\n", "kappa_rs", kappa_fixed.rs, "calibrated")
        @printf(io, "%-28s %18.10f  %s\n", "kappa_ns", kappa_fixed.ns, "calibrated")
        @printf(io, "%-28s %18.10f  %s\n", "eta_m",    P[1], "calibrated")
        @printf(io, "%-28s %18.10f  %s\n", "eta_rs",   P[2], "calibrated")
        @printf(io, "%-28s %18.10f  %s\n", "eta_ns",   etans_fixed, "fixed normalization")
        @printf(io, "%-28s %18.10f  %s\n", "sigma",    P[3], "calibrated")
        @printf(io, "%-28s %18.10f  %s\n", "theta_m",  eq["θm"][1], "recovered using 1950 VA shares")
        @printf(io, "%-28s %18.10f  %s\n", "theta_rs", eq["θrs"][1], "recovered using 1950 VA shares")
        @printf(io, "%-28s %18.10f  %s\n", "theta_ns", eq["θns"][1], "adding-up residual")
        for (label, value) in (("zeta0_m", Z0[1]), ("zeta0_rs", Z0[2]), ("zeta0_ns", Z0[3]),
                                                   ("gzeta_m", GZ[1]), ("gzeta_rs", GZ[2]), ("gzeta_ns", GZ[3]))
                # Under knots the services' zeta0 and gzeta are summaries of the knot path.
                status = (ZETA_KNOTS_DEFAULT !== nothing && !endswith(label, "_m")) ?
                         "summary of the knots (first knot, mean growth)" : "calibrated"
                @printf(io, "%-28s %18.10f  %s\n", label, value, status)
        end
        if ZETA_KNOTS_DEFAULT !== nothing
                println(io, "Services research productivity at the knots (BLOCKB_MODE=knots; zeta0 and gzeta above")
                println(io, "are the first knot and the mean growth to the last knot):")
                for (k, y) in enumerate(ZETA_KNOT_YEARS)
                        @printf(io, "  %d   zeta_rs %.6f   zeta_ns %.6f\n", y,
                                exp(ZETA_KNOTS_DEFAULT[1][k]), exp(ZETA_KNOTS_DEFAULT[2][k]))
                end
        end

        println(io)
        println(io, "DATA AND MODEL MOMENTS")
        println(io, "-"^112)
        @printf(io, "%-50s %15s %15s %15s  %s\n", "moment", "data", "model", "model-data", "role")
        for (label, data_value, model_value, role) in moment_rows
                if isnan(data_value)
                        @printf(io, "%-50s %15s %15.9f %15s  %s\n", label, "--", model_value, "--", role)
                else
                        @printf(io, "%-50s %15.9f %15.9f %+15.9f  %s\n",
                                        label, data_value, model_value, model_value - data_value, role)
                end
        end

        println(io)
        println(io, "NOTES")
        println(io, "- Goods LP in 1950 equals one only because both data and model LP series are indexed to 1950.")
        println(io, "  It is a normalization/check, not an additional identifying moment.")
        println(io, "- eta_ns is fixed at one; eta_m, rather than eta_ns, is calibrated.")
        println(io, "- theta_m and theta_rs are recovered internally using the two 1950 VA shares as anchors.")
        println(io, "- Model VA shares are omega_j s_j / sum_i omega_i s_i, omega_j = 1-(1-alpha) mu^y_j mu^x_j.")
        println(io, "- mu^y_j = LS_j/(alpha + (1-alpha) LS_j mu^x_j), with LS_j the labor share of sectoral value added.")
        println(io, "- Several blocks fit complete annual paths, so identification is not a one-parameter/one-scalar-moment count.")
        println(io, "- Sector prices are separately indexed to one in 1950 before their endpoint ratios are computed.")
        println(io, "- Data prices are relative to the private VA deflator (Tornqvist). The chained model rows divide")
        println(io, "  by a Fisher chain of model sector prices with VA weights; the final-good rows use the model numeraire.")
        println(io, "- Relative service prices are imposed through the production wedge tau^p, so they never test the model.")
        println(io, "- Model capital prices are data capital prices times D, the chained VA deflator in final-good units.")
        println(io, "- LP targets are the smoothing-spline trends of log LP: Goods fits the annual trend path, and")
        println(io, "  kappa_rs and kappa_ns equate the model and trend slopes of cumulative log LP since 1950.")
end
println("wrote ", REPORT_PATH)
@printf("KPSS endpoint growth factors used by the code: m=%.4f rs=%.4f ns=%.4f\n",
        PAT_GOODS, PAT_RI, PAT_NRI)
@printf("model N growth factors: m=%.4f rs=%.4f ns=%.4f\n",
        eq["Nmt"][T_smm]/eq["Nmt"][1], eq["Nrst"][T_smm]/eq["Nrst"][1], eq["Nnst"][T_smm]/eq["Nnst"][1])
