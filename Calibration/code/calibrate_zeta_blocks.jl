#==========================================================================================
 An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
 ------------------------------------------------------------------------------------------
 File:    calibrate_zeta_blocks.jl
 Purpose: Research-productivity block. Alternate the goods block (zeta0_m, g_m, kappa_m) and
          the services block (zeta0_rs, zeta0_ns, g_rs, g_ns) to a fixed point, given kappa
          and the demand parameters.
 ------------------------------------------------------------------------------------------
 Usage:   julia --project=. code/calibrate_zeta_blocks.jl "<k_m,k_rs,k_ns>" ["<eta_m,eta_rs,sigma>" \
              ["<z0_m,z0_rs,z0_ns>" "<g_m,g_rs,g_ns>"]]
          Normally called by calibrate_service_kappas.jl, once per kappa evaluation.
 ------------------------------------------------------------------------------------------
 Inputs:  data files listed in paths.jl (data_moments.xlsx, Cross_Citations CSVs)
 Outputs: TOML file at ENV CALIBRATION_RESULT_PATH (set by calibrate_service_kappas.jl):
          output/calibration/sweep_XX_kappa_zeta_evaluations/evaluation_NNN.toml
 ------------------------------------------------------------------------------------------
 Notes:
 Block-recursive calibration of the research productivity paths under free mobility
 (alloc=:stonegeary with zero floors) and no research wedges. Goods has
     zeta_m(t) = zeta0_m * exp(g_m * t);
 in the two services ln zeta_j is free at the knots 1950, 1960, ..., 2010, 2023 and linear in
 between (BLOCKB_MODE=knots, the baseline). BLOCKB_MODE=legacy gives the services constant
 growth, zeta0_j * exp(g_j * t), instead (Table J.1).

   Block A  (zeta0_m, g_m, kappa_m)           -> the 74-year goods LP trend path and the
                                                 74-year goods patent-stock path
   Block B  (services zeta at the knots)      -> the 74-year researcher share paths

 alternated to a fixed point (7 to 10 sweeps, parameter change ~1e-7).

 The split is forced by identification, not chosen for convenience. Under free
 mobility the researcher shares are
     L_j / L^s = X_j^{1/(1-phi)} / sum_k X_k^{1/(1-phi)},
 homogeneous of degree zero in X and hence in zeta, so the allocation responds only
 to the cross-sector RATIOS and cannot see a common scaling. The level reaches
 observables solely through the general-equilibrium chain idea stocks -> prices ->
 values, which is exactly what LP growth measures. Verified numerically: holding the
 ratios fixed and varying the level over four orders of magnitude leaves the
 scientist objective unchanged to six digits.

 Fitting the goods LP PATH rather than its 1950-2023 average matters for the same
 reason. A single average growth rate cannot separate zeta0_m from g_m; the two
 trade off exactly along zeta0 * exp(g * 73) = constant.

 kappa_m is calibrated in block A. On the LP path alone it is collinear with the
 goods research productivity, because measured productivity responds to the product
 of e_Y and the change in ln N_m while zeta_m moves that change. The goods
 patent-stock path, ln N_m(t)/N_m(1950) from the same series that sets the 1950
 knowledge shares, pins the idea path directly, which leaves kappa_m to map it into
 productivity. The two paths are equally weighted in log points. The first entry of
 the kappa argument is the starting value for kappa_m; kappa_rs and kappa_ns are
 held fixed.

 The blocks minimise different objectives, so the fixed point is a mutual best
 response and not the minimiser of any single criterion. Inference would need the
 stacked system.

==========================================================================================#

# Model specification (read by smm_ge_unified.jl when it is loaded).
ENV["ALLOC_MODE"]="stonegeary"; ENV["PHI_FIXED"]="1"; ENV["PHI_COMMON"]=get(ENV, "PHI_COMMON", "0.50")
length(ARGS) in (1, 2, 4) || error("usage: julia --project=. code/calibrate_zeta_blocks.jl \"k_m,k_rs,k_ns\" [\"eta_m,eta_rs,sigma\" [\"zeta0_m,zeta0_rs,zeta0_ns\" \"gzeta_m,gzeta_rs,gzeta_ns\"]]")
ENV["KAPPA_MODE"]="external"; ENV["KAPPA_EXT"] = ARGS[1]
import Pkg
import TOML
include(joinpath(@__DIR__, "paths.jl"))
Pkg.activate(CAL_PATHS.root)
include(joinpath(@__DIR__, "smm_ge_unified.jl"))
using Printf, Statistics
const P = length(ARGS)>1 ? parse.(Float64, split(ARGS[2],",")) : [0.72779181, 1.63321256, 0.20144784]
const TH0 = if length(ARGS) == 4
    zstart = parse.(Float64, split(ARGS[3],","))
    gstart = parse.(Float64, split(ARGS[4],","))
    length(zstart) == 3 && length(gstart) == 3 ||
        error("zeta warm start must contain two triples")
    all(>(0.0), zstart) || error("zeta0 warm-start values must be positive")
    vcat(log.(zstart), gstart, kappa_fixed.m, G2ZETA_DEFAULT[2], G2ZETA_DEFAULT[3])
else
    [log(1.0), log(1.4), log(0.73), -0.03, -0.03, -0.04, kappa_fixed.m, G2ZETA_DEFAULT[2], G2ZETA_DEFAULT[3]]
end
# Knot values start from ZETA_KNOTS when given, otherwise from the exponential path at the knots.
const NK = length(ZETA_KNOT_YEARS)
const TH0K = vcat(TH0, ZETA_KNOTS_DEFAULT !== nothing ? vcat(ZETA_KNOTS_DEFAULT...) :
    vcat([TH0[2] + TH0[5] * (y - 1950) for y in ZETA_KNOT_YEARS],
         [TH0[3] + TH0[6] * (y - 1950) for y in ZETA_KNOT_YEARS]))
# theta = [ln z0_m, ln z0_rs, ln z0_ns, gz_m, gz_rs, gz_ns, kappa_m, g2_rs, g2_ns,
#          ln zeta_rs at the knots, ln zeta_ns at the knots]
# BLOCKB_MODE=curv_ends: services research productivity has curvature,
#   ln zeta_j(t) = ln zeta0_j + g_j t + g2_j (t/73)^2, and block B fits the census researcher
#   shares (1950-2020 and 2023) with 1950 and 2023 weighted 1e3, so both are matched. With
#   constant growth the model front-loads the RI-S reallocation and cannot match 1950, 2023 and
#   the path together.
# BLOCKB_MODE=legacy: constant growth, fitted to the 74-year spline path (Table J.1).
# BLOCKB_MODE=ends: constant growth, exactly identified: zeta0_j and g_j match the 1950 and 2023
#   census researcher shares of RI-S and NRI-S (four equations, four unknowns).
# BLOCKB_MODE=knots (default, the baseline): services ln zeta_j free at ZETA_KNOT_YEARS (linear
#   between knots), fitted to the annual researcher path. zeta0_j and g_j of the services are
#   then reported as the first knot and the mean growth between the first and last knots.
const BLOCKB_MODE = get(ENV, "BLOCKB_MODE", "knots")
BLOCKB_MODE in ("curv_ends", "legacy", "ends", "knots") || error("BLOCKB_MODE must be curv_ends, legacy, ends or knots")
const KNOTS = BLOCKB_MODE == "knots"
const IA = [1, 4, 7]        # block A: goods zeta and kappa_m, goods LP and patent paths
const IB = KNOTS ? collect(10:9 + 2NK) : BLOCKB_MODE == "curv_ends" ? [2, 3, 5, 6, 8, 9] : [2, 3, 5, 6]   # block B: services zeta
const LO = [fill(log(1e-12),3); fill(-0.6,3); 0.05; -8.0; -8.0; fill(log(1e-12), 2NK)]
const HI = [fill(log(1e6),3);  fill( 0.6,3); 5.0; 8.0; 8.0; fill(log(1e6), 2NK)]
const CENSUS_YEARS = [1950, 1960, 1970, 1980, 1990, 2000, 2010, 2020, 2023]
const CENSUS_RS = [_getval(y, "RDSH_RI") for y in CENSUS_YEARS]
const CENSUS_NS = [_getval(y, "RDSH_NRI") for y in CENSUS_YEARS]
const CENSUS_W = [1e3; ones(length(CENSUS_YEARS) - 2); 1e3]
# Goods log-LP trend (smoothing spline, smm_ge_unified.jl), indexed to 1950: model and data LP
# levels are not comparable, so what is fitted is cumulative growth since 1950.
const LPM_DAT = _LP_TREND.m
# Goods patent stock, cumulative log growth since 1950 (PATENT_SERIES; default data_moments.xlsx, PATSTOCK_GOODS).
const NM_DAT = let i = [findfirst(==(y), VARIETY.years) for y in 1950:2023]
    any(isnothing, i) && error("variety stocks do not cover every year 1950-2023")
    log.(VARIETY.levels.m[i] ./ VARIETY.levels.m[i[1]])
end
_kap(th) = (m=th[7], rs=kappa_fixed.rs, ns=kappa_fixed.ns)
_zc(th) = BLOCKB_MODE == "curv_ends" ? (0.0, th[8], th[9]) : (0.0, 0.0, 0.0)
_kn(th) = KNOTS ? (th[10:9+NK], th[10+NK:9+2NK]) : nothing
eqof(th) = try _solved_or_nothing(_solve_GE_p(P, (exp(th[1]),exp(th[2]),exp(th[3])), (th[4],th[5],th[6]), (0.0,0.0);
                           kappa=_kap(th), g2zeta=_zc(th), knots=_kn(th))) catch; nothing end
# Services zeta0 and g reported in knots mode: first knot and mean growth to the last knot.
function _sync_knots!(th)
    KNOTS || return th
    for (i0, ig, off) in ((2, 5, 9), (3, 6, 9 + NK))
        th[i0] = th[off + 1]; th[ig] = (th[off + NK] - th[off + 1]) / (ZETA_KNOT_YEARS[end] - 1950)
    end
    th
end
function lpm_dev(eq)
    lp = eq["Ymt"] ./ eq["Lpmt"]
    [log(lp[t]/lp[1]) - LPM_DAT[t] for t in 1:T_smm]
end
nm_dev(eq) = [log(eq["Nmt"][t]/eq["Nmt"][1]) - NM_DAT[t] for t in 1:T_smm]
goods_dev(eq) = vcat(lpm_dev(eq), nm_dev(eq))
function sci_census(eq)
    idx = CENSUS_YEARS .- 1949
    L = [max(eq["Lst"][t], 1e-12) for t in idx]
    ([eq["Lsrst"][t] for t in idx] ./ L .- CENSUS_RS, [eq["Lsnst"][t] for t in idx] ./ L .- CENSUS_NS)
end
function sci_dev(eq)
    if BLOCKB_MODE == "ends"
        rrs, rns = sci_census(eq)
        return [rrs[1], rrs[end], rns[1], rns[end]]
    end
    if BLOCKB_MODE == "curv_ends"
        rrs, rns = sci_census(eq)
        return vcat(CENSUS_W .* rrs, CENSUS_W .* rns)
    end
    r = Float64[]
    for t in 1:T_smm
        L = max(eq["Lst"][t],1e-12)
        push!(r, eq["Lsrst"][t]/L - sci_share_rs_smm[t]); push!(r, eq["Lsnst"][t]/L - sci_share_ns_smm[t])
    end
    r
end
rms(x) = sqrt(mean(abs2, x))
# ZETA_STARTS=warm (default): with a warm start, each block keeps the solution reached from the current
# point and tries the alternative starts only if that solve fails. Choosing the start with the
# lowest block residual can move one block to another region, which the other block then
# follows (seen at phi = 0.75).
const WARM_ONLY = get(ENV, "ZETA_STARTS", "warm") == "warm" && length(ARGS) == 4   # multi: previous rule
# Solve the free coordinates `idx` of theta against `devfun`, holding the rest fixed.
function solve_block(th, idx, devfun, starts; budget=100.0)
    lb = LO[idx]; ub = HI[idx]; best = (Inf, nothing)
    for (n, s0) in enumerate(starts)
        WARM_ONLY && n > 1 && isfinite(best[1]) && break
        F = x -> begin
            t = copy(th); t[idx] = clamp.(x, lb, ub)
            eq = eqof(t); eq === nothing ? nothing : devfun(eq)
        end
        x, r, ok = _lm_solve(F, copy(s0), lb, ub; maxiter=250, tol=1e-14, budget=budget)
        v = (r !== nothing && all(isfinite, r)) ? rms(r) : Inf
        v < best[1] && (best = (v, clamp.(x, lb, ub)))
    end
    best
end
function main()
    th = copy(TH0K)
    @printf("kappa = %s\n\n", ARGS[1])
    @printf("%5s  %13s  %13s   %s\n", "sweep", "goodsLP RMSE", "sci RMSE", "max param change")
    local vA = Inf
    local vB = Inf
    for k in 1:10
        prev = copy(th)
        vA, xA = solve_block(th, IA, goods_dev,
                    (th[IA], [th[1]+log(0.5), th[4]+0.01, th[7]],
                     [th[1], th[4], 0.5*th[7]+1.0]); budget=120.0)
        xA === nothing && (println("  block A failed"); return)
        th[IA] = xA
        vB, xB = solve_block(th, IB, sci_dev,
                    KNOTS ? (th[IB],) :
                    (th[IB], BLOCKB_MODE != "curv_ends" ? [th[1]+log(1.41), th[1]+log(0.73), th[4], th[4]] :
                              [th[1]+log(1.41), th[1]+log(0.73), th[4], th[4], 0.0, 0.0]);
                    budget=KNOTS ? 1800.0 : BLOCKB_MODE == "curv_ends" ? 300.0 : 120.0)
        xB === nothing && (println("  block B failed"); return)
        th[IB] = xB
        _sync_knots!(th)
        d = maximum(abs.(th .- prev))
        @printf("%5d  %13.7f  %13.7f   %.3e\n", k, vA, vB, d)
        d < 1e-7 && break
    end
    eq = eqof(th); eq === nothing && return
    z0 = (exp(th[1]), exp(th[2]), exp(th[3])); g = (th[4], th[5], th[6]); km = th[7]
    lp_model = (
        m  = _ct_lp_ep(eq["Ymt"],  eq["Lpmt"]),
        rs = _ct_lp_ep(eq["Yrst"], eq["Lprst"]),
        ns = _ct_lp_ep(eq["Ynst"], eq["Lpnst"]),
    )
    @printf("\n  zeta0 = (%.6g, %.6g, %.6g)    ratios  1 : %.3f : %.3f\n", z0..., z0[2]/z0[1], z0[3]/z0[1])
    @printf("  gzeta = (%+.6f, %+.6f, %+.6f)\n", g...)
    @printf("  kappa_m = %.6f\n", km)
    @printf("  goods LP path RMSE = %.6f   goods patent path RMSE = %.6f   (log points)\n",
            rms(lpm_dev(eq)), rms(nm_dev(eq)))
    census_rmse = rms(vcat(sci_census(eq)...))
    @printf("  scientist RMSE at census points = %.6f   curvature g2_rs = %+.4f  g2_ns = %+.4f\n",
            census_rmse, _zc(th)[2], _zc(th)[3])
        result_path = get(ENV, "CALIBRATION_RESULT_PATH", "")
        if !isempty(result_path)
            open(result_path, "w") do io
                TOML.print(io, Dict(
                    "kappa" => Dict("m" => km, "rs" => kappa_fixed.rs, "ns" => kappa_fixed.ns),
                    "zeta0" => Dict("m" => z0[1], "rs" => z0[2], "ns" => z0[3]),
                    "gzeta" => Dict("m" => g[1], "rs" => g[2], "ns" => g[3]),
                    "g2zeta" => Dict("rs" => _zc(th)[2], "ns" => _zc(th)[3]),
                    "zeta_knots" => KNOTS ? Dict("years" => ZETA_KNOT_YEARS, "ln_rs" => th[10:9+NK],
                                                 "ln_ns" => th[10+NK:9+2NK]) : Dict{String,Any}(),
                    "fit" => Dict("goods_lp_rmse" => rms(lpm_dev(eq)),
                                  "goods_patent_rmse" => rms(nm_dev(eq)),
                                  "scientist_rmse" => census_rmse,
                                  "lp_model_m" => lp_model.m,
                                  "lp_model_rs" => lp_model.rs,
                                  "lp_model_ns" => lp_model.ns,
                                  "lp_data_m" => _LP_TARGETS[1],
                                  "lp_data_rs" => _LP_TARGETS[2],
                                  "lp_data_ns" => _LP_TARGETS[3],
                                  "lp_residual_rs" => lp_model.rs - _LP_TARGETS[2],
                                  "lp_residual_ns" => lp_model.ns - _LP_TARGETS[3]),
                ))
            end
            println("  machine-readable result: ", result_path)
        end
    println("\n  LP trend slope 1950-2023   (goods path targeted here; services set by the kappa loop)")
    for (nm,Y,Lp,k) in (("goods","Ymt","Lpmt",1),("RI-S","Yrst","Lprst",2),("NRI-S","Ynst","Lpnst",3))
        mo = _ct_lp_ep(eq[Y], eq[Lp])
        @printf("    %-6s model %.6f   data %.6f   %+7.2f%%   %s\n", nm, mo, _LP_TARGETS[k],
                100*(mo-_LP_TARGETS[k])/_LP_TARGETS[k], k==1 ? "path targeted" : "kappa target")
    end
    println("\n  scientist shares, model (census data)")
    rrs, rns = sci_census(eq)
    for (k, y) in enumerate(CENSUS_YEARS)
        @printf("    %d   RI-S %.4f (%.4f)    NRI-S %.4f (%.4f)\n", y,
                rrs[k] + CENSUS_RS[k], CENSUS_RS[k], rns[k] + CENSUS_NS[k], CENSUS_NS[k])
    end
end
main()
