#==========================================================================================
 An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
 ------------------------------------------------------------------------------------------
 File:    smm_ge_unified.jl
 Purpose: Model set-up shared by every calibration driver and by baseline.jl. Loads and
          transforms the data, fixes the structural constants, defines the targets, the
          GE wrapper _solve_GE_p and the numerical helpers. Loads equilibrium_path.jl.
          Not run on its own: the drivers set the specification through ENV variables
          (ALLOC_MODE, BLOCKB_MODE, PHI_COMMON, KAPPA_MODE, ...) and then include it.
 ------------------------------------------------------------------------------------------
 Inputs:  ../Empirical/output/dataset/data_moments.xlsx   (sheet DataMoments; or ENV DATA_MOMENTS)
          ../Cross_Citations/output/citations/spillover_backward1.csv
          ../Cross_Citations/output/robustness/stocks_summary.csv
          ../Cross_Citations/output/patents/class_stocks.csv  (only for PATENT_SERIES other than
                                                            the default data_moments)
 Outputs: none
==========================================================================================#

using DataFrames
using Printf
using Dates
using XLSX
using SmoothingSplines
using Random
using LinearAlgebra
using DelimitedFiles

include(joinpath(@__DIR__, "paths.jl"))

# Load the GE model.
include(joinpath(@__DIR__, "equilibrium_path.jl"))

#------------------------------------------------------------------------------
# FIXED PARAMETERS
#------------------------------------------------------------------------------

const T_smm       = 74          # 1950-2023 (number of periods in data window)
const T_ge        = 100         # total simulation horizon
const chi_smm     = 2.000
const rho_smm     = 0.040
const etans_fixed = 1.000
const gzetam_fix  = 0.0   # gzeta_m fixed at 0 (numeraire); gzeta_rs, gzeta_ns are estimated

# Fixed structural parameters (not part of SMM vector).
const gamma_fixed = parse(Float64, get(ENV, "GAMMA_COMMON", "0.5"))   # robustness: 0.25, 0.75, 0.90
# ALPHA: labor elasticity in variety production (also used to build μy from data). 0.70 follows
# Gollin's labor share of all income; 0.7344 removes land and inventories from the capital side
# (0.7/(0.7 + 0.3(1-s_j)), value-added weighted), since machines are reproducible assets only.
const alpha_fixed = parse(Float64, get(ENV, "ALPHA", "0.7344205648207067"))
const phi_m_fixed = 0.500   # scientist-labor elasticity in innovation (sector-neutral)
const phi_fixed   = phi_m_fixed   # sector-neutral φ (Roy efficiency index assumes common φ)
const theta_fr_fixed = 2.000   # Fréchet research-talent shape θ_F (fixed, not calibrated)
const T_m_fixed      = 1.000   # Fréchet goods-sector scale (normalization)

const W_BAD_PENALTY = 1.0   # penalty weight per non-finite simulated moment
# V_DAMP and ITER_MAX_V override the weight on the new V iterate and the iteration cap, for
# convergence diagnostics (phi = 0.75 stalls at the defaults).
const V_DAMP_SMM     = parse(Float64, get(ENV, "V_DAMP", "0.90"))
const ITER_MAX_V_SMM = parse(Int, get(ENV, "ITER_MAX_V", "45"))
const TOL_V_SMM      = 1e-3

#------------------------------------------------------------------------------
# LOAD DATA MOMENTS FROM EXCEL
#------------------------------------------------------------------------------

# Read from ../Empirical/output/dataset (see paths.jl). ENV DATA_MOMENTS points a run at
# another data file (e.g. an alternative capital-price series).
const DATA_MOMENTS_PATH = CAL_PATHS.data_moments
const df_data  = DataFrame(XLSX.readtable(DATA_MOMENTS_PATH, "DataMoments"))
const df_years = Int.(df_data[!, "Year"])
const smm_yrs  = collect(1950:2023)
const SPLINE_LAMBDA = 1e4
# SMOOTH_MODE=legacy: one spline at SPLINE_LAMBDA for inputs and targets alike.
# SMOOTH_MODE=role: exogenous inputs keep a spline, at SPLINE_LAMBDA_INPUT (6.25 removes
# business-cycle frequencies only, so endpoints stay close to the data and the post-2023
# tail, which holds 2023 values, is not biased); targets are the raw annual data, since
# the model has no cycle and path fits already average it out. Series observed only in
# census years are interpolated linearly, which keeps every observation and the adding-up
# of shares, and smoothing nine points with little sampling error has no purpose.
const SMOOTH_MODE = get(ENV, "SMOOTH_MODE", "legacy")
SMOOTH_MODE in ("legacy", "role") || error("SMOOTH_MODE must be legacy or role, got $SMOOTH_MODE")
const SPLINE_LAMBDA_INPUT = SMOOTH_MODE == "role" ?
    parse(Float64, get(ENV, "SPLINE_LAMBDA_INPUT", "6.25")) : SPLINE_LAMBDA

function _getval(year::Int, colname::String)
    idx = findfirst(==(year), df_years)
    return idx !== nothing ? Float64(coalesce(df_data[idx, colname], NaN)) : NaN
end

function _spline_fill(colname::String; lam::Float64=SPLINE_LAMBDA)
    col = Float64.(coalesce.(df_data[!, colname], NaN))
    obs = .!isnan.(col)
    count(obs) >= 2 || error("Need >=2 observations for spline fit: $colname")
    spl = fit(SmoothingSpline, Float64.(df_years[obs]), col[obs], lam)
    return predict(spl, Float64.(smm_yrs))
end

function _interp_fill(colname::String)
    col = Float64.(coalesce.(df_data[!, colname], NaN))
    obs = .!isnan.(col)
    yo, vo = df_years[obs], col[obs]
    (first(yo) <= first(smm_yrs) && last(yo) >= last(smm_yrs)) ||
        error("Interpolation would extrapolate $colname beyond $(first(yo))-$(last(yo))")
    map(smm_yrs) do y
        k = searchsortedlast(yo, y)
        yo[k] == y && return vo[k]
        w = (y - yo[k]) / (yo[k+1] - yo[k])
        (1 - w) * vo[k] + w * vo[k+1]
    end
end

# Series observed only in census years (researcher shares).
_census_path(colname::String) = SMOOTH_MODE == "role" ? _interp_fill(colname) : _spline_fill(colname)

# Data path a model path is fitted to.
function _target_path(colname::String)
    SMOOTH_MODE == "legacy" && return _spline_fill(colname)
    v = [_getval(y, colname) for y in smm_yrs]
    any(isnan, v) ? _interp_fill(colname) : v
end

#------------------------------------------------------------------------------
# INPUT TRENDS
#------------------------------------------------------------------------------
# INPUT_TRENDS selects how the exogenous inputs are smoothed; targets are unaffected.
#   spline  the paper's spline (penalty 1e4); capital prices on logs.
#   fp      fractional polynomial FP2 for every input: powers in {0, 0.5, 1, 2, 3} (negative
#           powers explode near 1950), x = year - 1949, a repeated power adds x^p log x; logs
#           for prices, indexes and population, levels for shares; no episodes excluded.
#   case    a form per input (the paper, Section 4.1). Default.
# Every fit uses 1950-2023 with equal weights and returns values on smm_yrs.
const INPUT_TRENDS = get(ENV, "INPUT_TRENDS", "case")
INPUT_TRENDS in ("spline", "fp", "case") || error("INPUT_TRENDS must be spline, fp or case, got $INPUT_TRENDS")
# INPUT_KNOTS: interior years (comma separated) at which labor quality and machine productivity
# may change growth rate, for INPUT_TRENDS=case. Growth is constant between the observations at
# 1950, the knots and 2023; with no knots (the baseline) it is constant over 1950-2023. After
# 2023 both grow at their last-segment rate. Appendix H.2 uses INPUT_KNOTS=1980.
const INPUT_KNOTS = let s = strip(get(ENV, "INPUT_KNOTS", ""))
    isempty(s) ? Int[] : sort(parse.(Int, strip.(split(s, ","))))
end
all(1950 .< INPUT_KNOTS .< 2023) || error("INPUT_KNOTS must lie strictly inside 1950-2023")
isempty(INPUT_KNOTS) || INPUT_TRENDS == "case" || error("INPUT_KNOTS requires INPUT_TRENDS=case")
const _INPUT_NODES = vcat(1950, INPUT_KNOTS, 2023)

_col(name::String) = Float64.(coalesce.(df_data[!, name], NaN))
function _obs(col, years; windows=Tuple{Int,Int}[])
    keep = [(1950 <= y <= 2023) && isfinite(v) && !any(a <= y <= b for (a, b) in windows)
            for (y, v) in zip(years, col)]
    return Float64.(years[keep]), col[keep]
end
_at(col, years, y) = col[findfirst(==(y), years)]
const _TT = (Float64.(smm_yrs) .- 1950.0) ./ 73.0

const _FP_POWERS = (0.0, 0.5, 1.0, 2.0, 3.0)
_fpterm(u, p) = p == 0 ? log.(u) : u .^ p
function _fp_design(u, p1, p2)
    t2 = _fpterm(u, p2)
    p2 == p1 && (t2 = t2 .* log.(u))
    return hcat(ones(length(u)), _fpterm(u, p1), t2)
end
function _fp2_trend(col, years)
    x, v = _obs(col, years)
    u = x .- 1949.0
    best = (Inf, 0.0, 0.0, Float64[])
    for (i, p1) in enumerate(_FP_POWERS), p2 in _FP_POWERS[i:end]
        X = _fp_design(u, p1, p2); β = X \ v
        r = sum(abs2, v .- X * β)
        r < best[1] && (best = (r, p1, p2, β))
    end
    return _fp_design(Float64.(smm_yrs) .- 1949.0, best[2], best[3]) * best[4]
end
_fp2_log(col, years) = exp.(_fp2_trend(log.(col), years))

# Constant growth between consecutive observations at _INPUT_NODES (1950, any knots, 2023).
function _endpoint_trend(col, years)
    lv = [log(_at(col, years, y)) for y in _INPUT_NODES]
    out = similar(_TT)
    for (i, y) in enumerate(smm_yrs)
        k = min(searchsortedlast(_INPUT_NODES, y), length(_INPUT_NODES) - 1)
        a, b = _INPUT_NODES[k], _INPUT_NODES[k+1]
        out[i] = lv[k] + (lv[k+1] - lv[k]) * (y - a) / (b - a)
    end
    return exp.(out)
end

# Quadratic through the 1950 and 2023 observations; only the curvature is estimated.
function _endq_trend(col, years; logs::Bool=false)
    f = logs ? log : identity
    x, v = _obs(col, years)
    y0, y1 = f(_at(col, years, 1950)), f(_at(col, years, 2023))
    t = (x .- 1950.0) ./ 73.0
    z = t .* (t .- 1.0)
    c = sum(z .* (f.(v) .- y0 .- (y1 - y0) .* t)) / sum(z .* z)
    out = y0 .+ (y1 - y0) .* _TT .+ c .* _TT .* (_TT .- 1.0)
    return logs ? exp.(out) : out
end


function _spline_fill(col::Vector{Float64}, years::Vector{Int};
                      λ::Float64=SPLINE_LAMBDA_INPUT,
                      target::Vector{Int}=smm_yrs,
                      min_obs::Int=2)::Vector{Float64}
    obs = .!isnan.(col)
    count(obs) >= min_obs || error("Need >=$(min_obs) non-NaN observations for spline fit")
    spl = fit(SmoothingSpline, Float64.(years[obs]), col[obs], λ)
    return predict(spl, Float64.(target))
end

# Annualized continuous-time growth from raw endpoints: dlog(x)/dt.
function _raw_ct_growth(colname::String, y0::Int, y1::Int)
    v0 = _getval(y0, colname)
    v1 = _getval(y1, colname)
    return log(v1 / v0) / (y1 - y0)
end

# Build a level series normalized at 1 from a sequence of log growth rates.
_norm_level_from_loggrowth(g::Vector{Float64}) = cumprod([1.0; exp.(g)])

function _normalize_triplet(a::Vector{Float64}, b::Vector{Float64}, c::Vector{Float64}; eps::Float64=1e-10)
    sa = clamp.(a, 0.0, 1.0)
    sb = clamp.(b, 0.0, 1.0)
    sc = clamp.(c, 0.0, 1.0)
    s = max.(sa .+ sb .+ sc, eps)
    return sa ./ s, sb ./ s, sc ./ s
end

#------------------------------------------------------------------------------
# EXOGENOUS INPUT PATHS FROM DATA (fed to solve_GE)
#------------------------------------------------------------------------------
# Naming convention:
#   *_smm : in-sample path (length T_smm, 1950-2023)
#   *_ge  : full model path (length T_ge), with a tail rule for t>T_smm

# KPSS market-value patent stocks. The CSV entries are natural logs. Endpoint
# growth factors identify κ, while 1950 level shares seed only N^Λ in Λ.
function _load_kpss_patents(path::String)
    isfile(path) || error("Missing KPSS patent-stock file: $path")
    lines = readlines(path)
    isempty(lines) && error("KPSS patent-stock file is empty: $path")
    header = split(chomp(lines[1]), ','; keepempty=true)
    expected = ["variant", "label", "sector", "1950", "1960", "1976",
                "1990", "2000", "2010", "2019", "2023"]
    header == expected || error("Unexpected stocks_summary.csv header: $(join(header, ','))")

    sector_map = Dict("Goods" => :m, "RI_Services" => :rs, "NRI_Services" => :ns)
    logs1950 = Dict{Symbol,Float64}()
    logs2023 = Dict{Symbol,Float64}()
    # Interior years, columns 5..10. Previously parsed and discarded. They are the
    # only patent information not already consumed by the kappa identity, which
    # uses the 1950 and 2023 endpoints only.
    logs_int = Dict{Symbol,Vector{Float64}}()
    for line in lines[2:end]
        startswith(line, "value_weighted,KPSS market value,") || continue
        fields = split(chomp(line), ','; keepempty=true)
        length(fields) == length(header) || error("Malformed KPSS row: $line")
        sector = get(sector_map, fields[3], nothing)
        sector === nothing && error("Unexpected KPSS sector: $(fields[3])")
        haskey(logs1950, sector) && error("Duplicate KPSS row for sector $(fields[3])")
        log1950 = parse(Float64, fields[4])
        log2023 = parse(Float64, fields[11])
        all(isfinite, (log1950, log2023)) || error("Non-finite KPSS log stock for $(fields[3])")
        logs1950[sector] = log1950
        logs2023[sector] = log2023
        li = [parse(Float64, fields[k]) for k in 5:10]
        all(isfinite, li) || error("Non-finite interior KPSS log stock for $(fields[3])")
        logs_int[sector] = li
    end
    all(haskey(logs1950, j) for j in (:m, :rs, :ns)) ||
        error("stocks_summary.csv must contain one KPSS market-value row for each sector")

    log0 = (m=logs1950[:m], rs=logs1950[:rs], ns=logs1950[:ns])
    log1 = (m=logs2023[:m], rs=logs2023[:rs], ns=logs2023[:ns])
    growth = (m=exp(log1.m-log0.m), rs=exp(log1.rs-log0.rs), ns=exp(log1.ns-log0.ns))
    growth_ratio = (m=exp(log1.m)/exp(log0.m),
                    rs=exp(log1.rs)/exp(log0.rs),
                    ns=exp(log1.ns)/exp(log0.ns))
    all(isapprox(growth[j], growth_ratio[j]; atol=1e-12, rtol=1e-12) for j in (:m, :rs, :ns)) ||
        error("KPSS growth factors fail exp(log2023)/exp(log1950) check")
    maxlog = max(log0.m, log0.rs, log0.ns)
    levels = (m=exp(log0.m-maxlog), rs=exp(log0.rs-maxlog), ns=exp(log0.ns-maxlog))
    total = levels.m + levels.rs + levels.ns
    shares = (m=levels.m/total, rs=levels.rs/total, ns=levels.ns/total)

    all(isfinite, values(growth)) && all(>(0.0), values(growth)) ||
        error("KPSS patent growth factors must be finite and positive")
    all(isfinite, values(shares)) && all(>(0.0), values(shares)) ||
        error("KPSS 1950 patent shares must be finite and positive")
    isapprox(sum(values(shares)), 1.0; atol=1e-12, rtol=1e-12) ||
        error("KPSS 1950 patent shares do not sum to one")
    # Interior growth factors N_j(t)/N_j(1950), same units as `growth`.
    iyears = (1960, 1976, 1990, 2000, 2010, 2019)
    igrowth = (m  = exp.(logs_int[:m]  .- log0.m),
               rs = exp.(logs_int[:rs] .- log0.rs),
               ns = exp.(logs_int[:ns] .- log0.ns))
    for j in (:m, :rs, :ns)
        gj = getfield(igrowth, j)
        length(gj) == length(iyears) || error("interior KPSS row for $j has wrong length")
        all(isfinite, gj) && all(>(0.0), gj) ||
            error("interior KPSS growth factors for $j must be finite and positive")
        issorted(gj) || @warn "interior KPSS stock for $j is not monotone" gj
    end
    return (; log1950=log0, log2023=log1, growth, shares,
              interior_years=iyears, interior_growth=igrowth)
end

# ---------------------------------------------------------------------------
# Annual sectoral variety stocks from `class_stocks.csv`.
#
# Supersedes the eight-year `stocks_summary.csv` panel. Two things changed that
# matter. The sectoral split now begins in 1977 rather than 1990, so there are 47
# annual observations of ln N_j - ln N_Goods instead of three or four intervals.
# And the headline series excludes business-method patents, which were carrying
# about 35 percent of measured NRI-S idea growth (ln N(2023)-ln N(1950) falls from
# 2.727 to 2.288 for NRI-S when they are removed, against 0.089 for RI-S).
#
# Rows before 1977 hold the AGGREGATE path only: ln N_j - ln N_Goods is constant
# there by construction. Use those years for the aggregate and for the 1950 initial
# condition, never as a relative-growth target.
# ---------------------------------------------------------------------------
function _load_variety_stocks(path::String; aggregation::String="ex_cpc_hist")
    isfile(path) || error("Missing variety-stock file: $path")
    lines = readlines(path)
    header = split(chomp(lines[1]), ','; keepempty=true)
    header == ["year","sector","N","lnN","aggregation"] ||
        error("Unexpected class_stocks.csv header: $(join(header, ','))")
    sector_map = Dict("Goods" => :m, "RI_Services" => :rs, "NRI_Services" => :ns)
    raw = Dict{Symbol,Dict{Int,Float64}}(:m=>Dict(), :rs=>Dict(), :ns=>Dict())
    for line in lines[2:end]
        f = split(chomp(line), ','; keepempty=true)
        length(f) == 5 || error("Malformed row: $line")
        f[5] == aggregation || continue
        j = get(sector_map, f[2], nothing); j === nothing && continue
        y = parse(Int, f[1]); v = parse(Float64, f[3])
        isfinite(v) && v > 0.0 || error("Non-positive stock at $y, $(f[2])")
        raw[j][y] = v
    end
    years = sort(collect(keys(raw[:m])))
    isempty(years) && error("No rows with aggregation=$aggregation")
    for j in (:rs, :ns)
        sort(collect(keys(raw[j]))) == years || error("Sector year grids disagree")
    end
    y0 = years[1]
    lev  = (m=[raw[:m][y] for y in years], rs=[raw[:rs][y] for y in years],
            ns=[raw[:ns][y] for y in years])
    grow = (m=lev.m ./ lev.m[1], rs=lev.rs ./ lev.rs[1], ns=lev.ns ./ lev.ns[1])
    rel  = (rs=log.(lev.rs ./ lev.m), ns=log.(lev.ns ./ lev.m))
    # First year the relative series actually moves; before it the split is imputed.
    split_year = y0
    for k in 2:length(years)
        if abs(rel.rs[k]-rel.rs[k-1]) > 1e-9 || abs(rel.ns[k]-rel.ns[k-1]) > 1e-9
            split_year = years[k]; break
        end
    end
    tot = lev.m[1] + lev.rs[1] + lev.ns[1]
    shares = (m=lev.m[1]/tot, rs=lev.rs[1]/tot, ns=lev.ns[1]/tot)
    isapprox(shares.m+shares.rs+shares.ns, 1.0; atol=1e-12) ||
        error("1950 variety shares do not sum to one")
    return (; years, levels=lev, growth=grow, rel, shares, split_year, aggregation)
end

# Variety stocks from data_moments.xlsx (PATSTOCK_*), the default: the same source as every
# other calibration input. Same fields as _load_variety_stocks, on 1950-2023.
function _load_variety_moments()
    yrs = df_years
    keep = [y in 1950:2023 for y in yrs]
    col(c) = Float64.(coalesce.(df_data[keep, c], NaN))
    lev = (m=col("PATSTOCK_GOODS"), rs=col("PATSTOCK_RI"), ns=col("PATSTOCK_NRI"))
    all(v -> all(isfinite, v) && all(>(0.0), v), values(lev)) ||
        error("PATSTOCK_* in data_moments.xlsx must be positive for 1950-2023")
    years = yrs[keep]
    grow = (m=lev.m ./ lev.m[1], rs=lev.rs ./ lev.rs[1], ns=lev.ns ./ lev.ns[1])
    rel  = (rs=log.(lev.rs ./ lev.m), ns=log.(lev.ns ./ lev.m))
    split_year = years[1]
    for k in 2:length(years)
        if abs(rel.rs[k]-rel.rs[k-1]) > 1e-9 || abs(rel.ns[k]-rel.ns[k-1]) > 1e-9
            split_year = years[k]; break
        end
    end
    tot = lev.m[1] + lev.rs[1] + lev.ns[1]
    shares = (m=lev.m[1]/tot, rs=lev.rs[1]/tot, ns=lev.ns[1]/tot)
    return (; years, levels=lev, growth=grow, rel, shares, split_year, aggregation="data_moments")
end

# Which patent series drives the ENTIRE patent block: the variety stocks, the
# kappa endpoints, the spillover weights and the inner-loop conditions. One switch
# so they can never come from different series.
#   data_moments                            DEFAULT. PATSTOCK_* in data_moments.xlsx, which is
#                                           the ex_cpc_hist series exported by the replication
#                                           package (step 10b)
#   ex_cpc_hist                             ex_cpc accumulated without obsolescence
#                                           (blueprints do not depreciate in the model), with
#                                           the pre-1976 sector composition measured on the
#                                           USPC route (Cross_Citations, code/05_patent_stocks.py
#                                           and code/04_hpdf_sector_flows.py)
#   ex_cpc / ex_bm / harmonised / baseline  annual counts from class_stocks.csv, 15 percent
#                                           obsolescence, 1976 shares held before 1976
#   kpss_vw                                 value-weighted, old eight-year panel
const PATENT_SERIES = get(ENV, "PATENT_SERIES", "data_moments")
const N_INIT = get(ENV, "N_INIT", "one")   # "one" or "shares"
# :stonegeary solves the allocation; :data imposes the observed researcher path and
# reports the wage dispersion it implies, which is the required distortion.
const ALLOC_MODE = Symbol(get(ENV, "ALLOC_MODE", "stonegeary"))

# The count aggregations all live in class_stocks.csv; kpss_vw falls back to the
# old eight-year panel and is handled separately below.
const VARIETY = PATENT_SERIES == "data_moments" ? _load_variety_moments() :
    _load_variety_stocks(CAL_PATHS.class_stocks;
                         aggregation = PATENT_SERIES == "kpss_vw" ? "ex_cpc" : PATENT_SERIES)
const VARIETY_SPLIT_YEAR = VARIETY.split_year   # relative targets are valid from here

const KPSS_PATENTS = _load_kpss_patents(CAL_PATHS.stocks_summary)
const KPSS_PATENT_GROWTH = KPSS_PATENTS.growth
# lambda_N0 now comes from the new annual file. The service sectors swap rank in
# 1950: NRI-S is the larger of the two here (0.01726) where the old panel made it
# the smaller (0.00867). This scales every Lambda_j, so it moves zeta and the wedge.
const KPSS_PATENT_SHARES_1950 = PATENT_SERIES == "kpss_vw" ? KPSS_PATENTS.shares : VARIETY.shares
const _NINI = N_INIT == "shares" ?
    (KPSS_PATENT_SHARES_1950.m, KPSS_PATENT_SHARES_1950.rs, KPSS_PATENT_SHARES_1950.ns) :
    (1.0, 1.0, 1.0)

# Raw columns from data_moments.xlsx (NaN where missing).
const profit_share_GOODS = Float64.(coalesce.(df_data[!, "profit_share_GOODS"], NaN))
const profit_share_RI    = Float64.(coalesce.(df_data[!, "profit_share_RI"],    NaN))
const profit_share_NRI   = Float64.(coalesce.(df_data[!, "profit_share_NRI"],   NaN))
const CAPPRICE_GOODS     = Float64.(coalesce.(df_data[!, "CAPPRICE_GOODS"],     NaN))
const CAPPRICE_RI        = Float64.(coalesce.(df_data[!, "CAPPRICE_RI"],        NaN))
const CAPPRICE_NRI       = Float64.(coalesce.(df_data[!, "CAPPRICE_NRI"],       NaN))
const LABQLTY_GOODS      = Float64.(coalesce.(df_data[!, "LABQLTY_GOODS"],      NaN))
const LABQLTY_RI         = Float64.(coalesce.(df_data[!, "LABQLTY_RI"],         NaN))
const LABQLTY_NRI        = Float64.(coalesce.(df_data[!, "LABQLTY_NRI"],        NaN))
# LS_NET="s_m,s_rs,s_ns": labor share of value added net of land and inventory income, with s_j the
# land+inventory share of capital income (BLS 1987-2023 means, held over the whole sample):
# LS' = LS / (1 - s_j (1 - LS)). The model's machines are reproducible assets only.
# LS_NET=none switches the adjustment off (alpha 0.70 setting).
const _LS_NET = let v = get(ENV, "LS_NET", "0.140,0.048,0.207")
    v in ("", "none") ? (0.0, 0.0, 0.0) : Tuple(parse.(Float64, split(v, ",")))
end
_ls_net(v, sj) = v ./ (1 .- sj .* (1 .- v))
const LABVASHARE_GOODS   = _ls_net(Float64.(coalesce.(df_data[!, "LABVASHARE_GOODS"],   NaN)), _LS_NET[1])
const LABVASHARE_RI      = _ls_net(Float64.(coalesce.(df_data[!, "LABVASHARE_RI"],      NaN)), _LS_NET[2])
const LABVASHARE_NRI     = _ls_net(Float64.(coalesce.(df_data[!, "LABVASHARE_NRI"],     NaN)), _LS_NET[3])
# Endpoint growth now from the annual ex_cpc file, so kappa and lambda_N0 come from
# the SAME source. Excluding business methods changes NRI-S dramatically: the
# 1950-2023 growth factor falls from 52.46 on the old value-weighted panel to about
# 9.94 here. kappa_ns rises correspondingly, which is the intended consequence: each
# NRI-S idea does more for output, so fewer are needed and zeta_ns falls.
# PATENT_SERIES switches the whole patent block at once: the kappa endpoints, the
# spillover weights, and the inner-loop conditions. They must never come from
# different series.
#
#   ex_cpc    annual counts excluding business methods, sectoral split from 1977
#   kpss_vw   KPSS value-weighted, eight years only, split from 1990, business
#             methods included
#
# The choice matters a great deal for the allocation, because N_j is the denominator
# of profit per variety. Dilution relative to goods at 2023 is 8.83 (rs) and 1.20
# (ns) under ex_cpc, against 13.50 and 6.98 under kpss_vw: NRI-S is squeezed 5.8
# times harder under the value-weighted series, RI-S only 1.5 times harder.
const _VIDX_2023         = findfirst(==(2023), VARIETY.years)
const _VIDX_1990         = findfirst(==(1990), VARIETY.years)
const PAT_GOODS          = PATENT_SERIES == "kpss_vw" ? KPSS_PATENTS.growth.m  : VARIETY.growth.m[_VIDX_2023]
const PAT_RI             = PATENT_SERIES == "kpss_vw" ? KPSS_PATENTS.growth.rs : VARIETY.growth.rs[_VIDX_2023]
const PAT_NRI            = PATENT_SERIES == "kpss_vw" ? KPSS_PATENTS.growth.ns : VARIETY.growth.ns[_VIDX_2023]
const PAT_G90 = PATENT_SERIES == "kpss_vw" ?
    (m=KPSS_PATENTS.interior_growth.m[3], rs=KPSS_PATENTS.interior_growth.rs[3],
     ns=KPSS_PATENTS.interior_growth.ns[3]) :
    (m=VARIETY.growth.m[_VIDX_1990], rs=VARIETY.growth.rs[_VIDX_1990],
     ns=VARIETY.growth.ns[_VIDX_1990])

# Interior patent stocks, N_j(t)/N_j(1950) at 1960, 1976, 1990, 2000, 2010, 2019.
# UNUSED by the kappa identity, which consumes only the 1950 and 2023 endpoints, so
# these are the one genuinely independent block of patent information available.
const PAT_INTERIOR_YEARS = KPSS_PATENTS.interior_years
const PAT_INTERIOR       = KPSS_PATENTS.interior_growth
const pop                = Float64.(coalesce.(df_data[!, "pop"],                NaN))

# Relative value-added prices (in-sample targets, length T_smm).
# solve_GE recovers the price wedges to match these for t≤T_smm and holds the
# wedges fixed at their 2023 values thereafter; no _ge extension is fed.
const tail_ge          = T_ge - T_smm

# Knowledge-spillover elasticities δ_ji, constant parameters. Each is the mean over citing years
# 1978-2020 of the probability that a patent in sector j cites a patent in sector i (business
# methods excluded; 1977 is a truncated citation window and citing years after 2020 are
# incomplete), with each citing sector's row renormalized to sum to one.
const spill_raw, spill_header = readdlm(CAL_PATHS.spillovers, ',', header=true)
const spill_cols = String.(vec(spill_header))
_spill_col(name::String) = Float64.(spill_raw[:, findfirst(==(name), spill_cols)])
const spill_years = Int.(round.(_spill_col("start")))
const SPILL_START = 1978
const SPILL_END   = 2020
function _spill_row(r::String)
    keep = (spill_years .>= SPILL_START) .& (spill_years .<= SPILL_END)
    m = Tuple(sum(_spill_col("delta_$(r)$(j)")[keep]) / count(keep) for j in ("m", "r", "n"))
    return m ./ sum(m)
end
const dMM, dMR, dMN = _spill_row("m")
const dRM, dRR, dRN = _spill_row("r")
const dNM, dNR, dNN = _spill_row("n")
@assert all(0.0 .< (dMM, dMR, dMN, dRM, dRR, dRN, dNM, dNR, dNN) .< 1.0)
# The solver takes per-period arrays; the elasticities are the same in every period.
const dMM_ge = fill(dMM, T_ge); const dMR_ge = fill(dMR, T_ge); const dMN_ge = fill(dMN, T_ge)
const dRM_ge = fill(dRM, T_ge); const dRR_ge = fill(dRR, T_ge); const dRN_ge = fill(dRN, T_ge)
const dNM_ge = fill(dNM, T_ge); const dNR_ge = fill(dNR, T_ge); const dNN_ge = fill(dNN, T_ge)

# Service/goods relative prices are inputs: the production wedges reproduce them exactly.
const pm_m_smm         = _spline_fill("VAPRICE_GOODS"; lam=SPLINE_LAMBDA_INPUT)
const pm_rs_smm        = _spline_fill("VAPRICE_RI";    lam=SPLINE_LAMBDA_INPUT)
const pm_ns_smm        = _spline_fill("VAPRICE_NRI";   lam=SPLINE_LAMBDA_INPUT)
# The two relative service prices follow the same splines as the price targets.
const prs_pm_smm       = pm_rs_smm ./ pm_m_smm
const pns_pm_smm       = pm_ns_smm ./ pm_m_smm
# The goods price relative to the aggregate is a target, and the ratios are reported.
const pm_m_tgt         = _target_path("VAPRICE_GOODS")
const pm_rs_tgt        = _target_path("VAPRICE_RI")
const pm_ns_tgt        = _target_path("VAPRICE_NRI")

# Machine-producer markups μ^x_j (Appendix F). Baseline: each sector's gross markup 1/μ^x
# follows the OLS quadratic trend in time of the measured series 1/(1-ψ) over 1950-2023 and
# stays at its 2023 value afterwards. MUX_MODE=mean holds it at the 1950-2023 mean (Table J.1).
const MUX_MODE = get(ENV, "MUX_MODE", "quad")
MUX_MODE in ("mean", "quad") || error("MUX_MODE must be mean or quad")
function _mux_const(ps)
    x, v = _obs(ps, df_years)
    return length(v) / sum(1 ./ (1 .- v))
end
const μxm  = _mux_const(profit_share_GOODS)
const μxrs = _mux_const(profit_share_RI)
const μxns = _mux_const(profit_share_NRI)
function _mux_quad(ps)
    x, v = _obs(ps, df_years)
    keep = (x .>= smm_yrs[1]) .& (x .<= smm_yrs[end])
    t = (x[keep] .- smm_yrs[1]) ./ (smm_yrs[end] - smm_yrs[1])
    b = hcat(ones(length(t)), t, t .^ 2) \ (1.0 ./ (1.0 .- v[keep]))
    τ = collect(0:T_smm-1) ./ (T_smm - 1)
    gross = b[1] .+ b[2] .* τ .+ b[3] .* τ .^ 2
    all(>(1.0), gross) || error("quadratic machine markup trend falls below one")
    return 1.0 ./ gross
end
# The solver takes per-period arrays.
const μxmt_smm  = MUX_MODE == "quad" ? _mux_quad(profit_share_GOODS) : fill(μxm,  T_smm)
const μxrst_smm = MUX_MODE == "quad" ? _mux_quad(profit_share_RI)    : fill(μxrs, T_smm)
const μxnst_smm = MUX_MODE == "quad" ? _mux_quad(profit_share_NRI)   : fill(μxns, T_smm)

# Final-producer markups μ^y_j, constant parameters. Sectoral value added consolidates the variety
# producers with the machine producers that supply them, VA_j = p_jY_j - I_j^x = ω_j p_jY_j with
# ω_j = 1-(1-α)μ^y_j μ^x_j. Labor receives αμ^y_j p_jY_j, so its VA share is LS_j = αμ^y_j/ω_j,
# which inverts to μ^y_j = LS_j / (α + (1-α) LS_j μ^x_j). Each sector's gross markup 1/μ^y is the
# 1950-2023 mean of the values implied by the measured labor shares, so model labor shares are
# constant.
_μy_from_labva(ls, μx, α) = ls / (α + (1.0 - α) * ls * μx)
# μx is the year's machine markup on the 1950-2023 path.
function _muy_const(ls_raw, μx_path)
    x, v = _obs(ls_raw, df_years)
    keep = (x .>= smm_yrs[1]) .& (x .<= smm_yrs[end])
    μx_obs = μx_path[Int.(x[keep]) .- smm_yrs[1] .+ 1]
    gross = 1.0 ./ _μy_from_labva.(v[keep], μx_obs, alpha_fixed)
    return length(gross) / sum(gross)
end
const μym  = _muy_const(LABVASHARE_GOODS, μxmt_smm)
const μyrs = _muy_const(LABVASHARE_RI,    μxrst_smm)
const μyns = _muy_const(LABVASHARE_NRI,   μxnst_smm)
all(x -> 0.0 < x < 1.0, (μym, μyrs, μyns)) || error("μ^y implied by the labor VA shares lies outside (0,1)")

const μymt_smm = fill(μym,  T_smm); const μyrst_smm = fill(μyrs, T_smm); const μynst_smm = fill(μyns, T_smm)
const μxmt_ge  = vcat(μxmt_smm,  fill(μxmt_smm[end],  T_ge - T_smm))
const μxrst_ge = vcat(μxrst_smm, fill(μxrst_smm[end], T_ge - T_smm))
const μxnst_ge = vcat(μxnst_smm, fill(μxnst_smm[end], T_ge - T_smm))
const μymt_ge  = fill(μym,  T_ge);  const μyrst_ge  = fill(μyrs, T_ge);  const μynst_ge  = fill(μyns, T_ge)

# Gross-revenue-to-value-added factors ω_j.
_ω_va(μy, μx, α) = 1.0 - (1.0 - α) * μy * μx
const ωmt_ge                = _ω_va.(μymt_ge,  μxmt_ge,  alpha_fixed)
const ωrst_ge               = _ω_va.(μyrst_ge, μxrst_ge, alpha_fixed)
const ωnst_ge               = _ω_va.(μynst_ge, μxnst_ge, alpha_fixed)

# z_j paths from sectoral p_k and μx.
# Capital prices are smoothed in logs (the default; CAPPRICE_LOG_SPLINE=0 restores levels). A
# spline in levels undershoots at the end of a series that falls by an order of magnitude:
# the quality-adjusted RI-S price (CJK) falls 30-fold, and its level spline ends at 0.001
# against 0.033 in the data, doubling g_z,rs.
const CAPPRICE_LOG_SPLINE = get(ENV, "CAPPRICE_LOG_SPLINE", "1") == "1"
function _pk_fill(v)
    INPUT_TRENDS == "fp" && return _fp2_log(v, df_years)
    INPUT_TRENDS == "case" && return _endpoint_trend(v, df_years)
    return CAPPRICE_LOG_SPLINE ? exp.(_spline_fill(log.(v), df_years)) : _spline_fill(v, df_years)
end
const pₖmt_s                = _pk_fill(CAPPRICE_GOODS)
const pₖrst_s               = _pk_fill(CAPPRICE_RI)
const pₖnst_s               = _pk_fill(CAPPRICE_NRI)
const pₖmt                  = pₖmt_s  ./ pₖmt_s[1]
const pₖrst                 = pₖrst_s ./ pₖrst_s[1]
const pₖnst                 = pₖnst_s ./ pₖnst_s[1]
# The data deflate capital prices by the private VA deflator (Tornqvist); machines in the
# model cost final goods, whose price drifts against that deflator. D(t) is the model's chained VA
# deflator in final-good units (1950 = 1), so the model capital price is the data price times D.
# PK_DEFL holds ln D at _INPUT_NODES, linear in between. It is an equilibrium object:
# run_full_calibration.jl updates it every sweep, and baseline.jl reads the calibrated value.
# The default is the fixed point of the published calibration (nodes 1950 and 2023); with
# INPUT_KNOTS set and PK_DEFL unset it starts from zero.
const PK_DEFL_BASELINE = "0,0.06380282626824207"
const PK_DEFL = let s = strip(get(ENV, "PK_DEFL", isempty(INPUT_KNOTS) ? PK_DEFL_BASELINE : ""))
    isempty(s) ? zeros(length(_INPUT_NODES)) : parse.(Float64, strip.(split(s, ",")))
end
length(PK_DEFL) == length(_INPUT_NODES) || error("PK_DEFL needs one value per input node $(_INPUT_NODES)")
const _PK_DEFL_PATH = let d = zeros(T_smm)
    for k in 1:length(_INPUT_NODES)-1
        a, b = _INPUT_NODES[k], _INPUT_NODES[k+1]
        for y in a:b
            d[y-1949] = PK_DEFL[k] + (PK_DEFL[k+1] - PK_DEFL[k]) * (y - a) / (b - a)
        end
    end
    exp.(d)
end
const zmt_s                 = 1.0 ./ (pₖmt  .* _PK_DEFL_PATH .* μxmt_smm)
const zrst_s                = 1.0 ./ (pₖrst .* _PK_DEFL_PATH .* μxrst_smm)
const znst_s                = 1.0 ./ (pₖnst .* _PK_DEFL_PATH .* μxnst_smm)

const _dt_smm               = T_smm - 1

# Piecewise constant log growth through the values at _INPUT_NODES, as annual rates on the
# smm grid (length T_smm - 1). With no interior knots it is one constant rate.
function _pw_growth(x::Vector{Float64})
    g = zeros(_dt_smm)
    for k in 1:length(_INPUT_NODES)-1
        a, b = _INPUT_NODES[k], _INPUT_NODES[k+1]
        ia, ib = a - 1949, b - 1949
        g[ia:ib-1] .= log(x[ib] / x[ia]) / (b - a)
    end
    return g
end
const _gzm_path             = _pw_growth(zmt_s)
const _gzrs_path            = _pw_growth(zrst_s)
const _gzns_path            = _pw_growth(znst_s)
# Growth after 2023: the last segment's rate (the 1950-2023 rate without knots).
const gzm_smm               = _gzm_path[end]
const gzrs_smm              = _gzrs_path[end]
const gzns_smm              = _gzns_path[end]

const zmt_smm               = _norm_level_from_loggrowth(_gzm_path)
const zrst_smm              = _norm_level_from_loggrowth(_gzrs_path)
const znst_smm              = _norm_level_from_loggrowth(_gzns_path)
# Tail rule: extend z_j with constant continuous-time growth g_zj.
const zmt_ge                = vcat(zmt_smm, zmt_smm[end] .* cumprod(exp.(fill(gzm_smm, tail_ge))))
const zrst_ge               = vcat(zrst_smm, zrst_smm[end] .* cumprod(exp.(fill(gzrs_smm, tail_ge))))
const znst_ge               = vcat(znst_smm, znst_smm[end] .* cumprod(exp.(fill(gzns_smm, tail_ge))))

# hbar_j paths.
_h_path(v) = INPUT_TRENDS == "spline" ? _spline_fill(v, df_years) :
             INPUT_TRENDS == "fp" ? _fp2_log(v, df_years) : _endpoint_trend(v, df_years)
const hbarmt_s              = _h_path(LABQLTY_GOODS)
const hbarrst_s             = _h_path(LABQLTY_RI)
const hbarnst_s             = _h_path(LABQLTY_NRI)
const _ghm_path             = _pw_growth(hbarmt_s)
const _ghrs_path            = _pw_growth(hbarrst_s)
const _ghns_path            = _pw_growth(hbarnst_s)
const ghbarm_smm            = _ghm_path[end]
const ghbarrs_smm           = _ghrs_path[end]
const ghbarns_smm           = _ghns_path[end]

const hbarmt_smm            = _norm_level_from_loggrowth(_ghm_path)
const hbarrst_smm           = _norm_level_from_loggrowth(_ghrs_path)
const hbarnst_smm           = _norm_level_from_loggrowth(_ghns_path)
# Tail rule: extend hbar_j with constant continuous-time growth g_hbar_j.
const hbarmt_ge             = vcat(hbarmt_smm, hbarmt_smm[end] .* cumprod(exp.(fill(ghbarm_smm, tail_ge))))
const hbarrst_ge            = vcat(hbarrst_smm, hbarrst_smm[end] .* cumprod(exp.(fill(ghbarrs_smm, tail_ge))))
const hbarnst_ge            = vcat(hbarnst_smm, hbarnst_smm[end] .* cumprod(exp.(fill(ghbarns_smm, tail_ge))))

# Sector-specific love-of-variety parameters from the 1950-2023 LP identity:
#   LP_j(t) = p_j(t)^β N_j(t)^eY_j(t) A_j(t),  β=(1-α)/α,
#   eY_j(t) = (κ_j + (1-μx_j(t))/μx_j(t)) β.
# N-growth factors use the KPSS market-value log-stock endpoints.
const _KAPPA_START_YEAR = 1950
const _KAPPA_END_YEAR   = 2023

function _year_index(years::Vector{Int}, year::Int, source::String)::Int
    idx = findfirst(==(year), years)
    idx === nothing && error("Missing $year from $source")
    return idx
end

const _kappa_smm_t0  = _year_index(smm_yrs,  _KAPPA_START_YEAR, "the SMM year grid")
const _kappa_smm_t1  = _year_index(smm_yrs,  _KAPPA_END_YEAR,   "the SMM year grid")

function _require_positive_finite(x::Real, label::String)::Float64
    value = Float64(x)
    isfinite(value) && value > 0.0 || error("$label must be finite and strictly positive; got $value")
    return value
end

function _backout_kappa(sector::String, lp_col::String,
                        N_growth_factor::Float64, price_path::Vector{Float64},
                        μx_path::Vector{Float64}, μy_path::Vector{Float64},
                        z_path::Vector{Float64}, hbar_path::Vector{Float64})
    lp0 = _require_positive_finite(_getval(_KAPPA_START_YEAR, lp_col), "$lp_col[$(_KAPPA_START_YEAR)]")
    lp1 = _require_positive_finite(_getval(_KAPPA_END_YEAR, lp_col),   "$lp_col[$(_KAPPA_END_YEAR)]")
    N0 = 1.0
    N1 = _require_positive_finite(N_growth_factor, "$sector N-growth factor")
    p0 = _require_positive_finite(price_path[_kappa_smm_t0], "$sector price[$(_KAPPA_START_YEAR)]")
    p1 = _require_positive_finite(price_path[_kappa_smm_t1], "$sector price[$(_KAPPA_END_YEAR)]")
    μx0 = _require_positive_finite(μx_path[_kappa_smm_t0], "$sector μx[$(_KAPPA_START_YEAR)]")
    μx1 = _require_positive_finite(μx_path[_kappa_smm_t1], "$sector μx[$(_KAPPA_END_YEAR)]")
    A0 = _require_positive_finite(A_(μy_path[_kappa_smm_t0], μx0, z_path[_kappa_smm_t0],
                                    alpha_fixed, hbar_path[_kappa_smm_t0]),
                                  "$sector A[$(_KAPPA_START_YEAR)]")
    A1 = _require_positive_finite(A_(μy_path[_kappa_smm_t1], μx1, z_path[_kappa_smm_t1],
                                    alpha_fixed, hbar_path[_kappa_smm_t1]),
                                  "$sector A[$(_KAPPA_END_YEAR)]")

    β = (1.0 - alpha_fixed) / alpha_fixed
    dln_lp = log(lp1 / lp0)
    dln_p  = log(p1 / p0)
    dln_N  = log(N1 / N0)
    dln_A  = log(A1 / A0)
    abs(dln_N) > 1e-10 || error("$sector N index has no $(_KAPPA_START_YEAR)-$(_KAPPA_END_YEAR) log change; κ is unidentified")

    markup_component_end = (1.0 - μx1) / μx1
    κ = (dln_lp - β*dln_p - dln_A) / (β*dln_N) - markup_component_end
    # Sign is checked only where the backed-out value is used (kappa_fixed below): with
    # KAPPA_MODE=external it is a diagnostic, and a nonpositive value then means the
    # measured non-idea terms already exceed LP growth over the sample.
    isfinite(κ) || error("Backed-out κ_$sector must be finite; got $κ")

    reconstructed_dln_lp = β*dln_p + β*(κ + markup_component_end)*dln_N + dln_A
    isapprox(reconstructed_dln_lp, dln_lp; atol=1e-12, rtol=1e-12) ||
        error("κ_$sector failed the $(_KAPPA_START_YEAR)-$(_KAPPA_END_YEAR) LP accounting check")
    return (; κ, N0, N1, dln_lp, dln_p, dln_A, reconstructed_dln_lp)
end

const _kappa_m_fit = _backout_kappa("m", "LP_GOODS", PAT_GOODS, pm_m_smm,
                                    μxmt_smm, μymt_smm, zmt_smm, hbarmt_smm)
const _kappa_rs_fit = _backout_kappa("rs", "LP_RI", PAT_RI, pm_rs_smm,
                                     μxrst_smm, μyrst_smm, zrst_smm, hbarrst_smm)
const _kappa_ns_fit = _backout_kappa("ns", "LP_NRI", PAT_NRI, pm_ns_smm,
                                     μxnst_smm, μynst_smm, znst_smm, hbarnst_smm)
# kappa_j: BACKED OUT from the LP identity, or FIXED EXTERNALLY.
#
# Backing out is what forced kappa_rs to 0.362. The identity has to reconcile 73x
# patent growth in R&D-intensive services with modest measured productivity growth,
# and the only way is a low output elasticity of the idea stock. That in turn makes
# e_Y_rs = (kappa_rs + markup)(1-alpha)/alpha only 0.32, so profit per variety falls
# as N^-0.68 and V_rs collapses by a factor of 3.3 over the sample. The dilution is
# therefore not a modelling choice; it is the identity's answer to the two data
# series, and it is what stops researchers moving into RI-S in the model.
#
# KAPPA_MODE=external breaks that. kappa is set from outside, the LP-growth moments
# move into the objective, and the identity becomes a testable restriction rather
# than a definition. KAPPA_EXT sets the values; the default makes services share the
# goods elasticity, which is close to the kappa_rs = 1.94 that removes RI-S dilution
# entirely.
# DEFAULTS ARE NOW THE CHOSEN SPECIFICATION, so a JuliaHub run (which carries no
# environment) uses it. Set KAPPA_MODE=backout / ZFIT_TARGET=patents to recover
# the earlier configurations.
const KAPPA_MODE = get(ENV, "KAPPA_MODE", "external")
# kappa_rs = 0.362 is the value the LP identity backs out of the patent data. At the
# common 1.81 the inner system cannot hold LP growth and the scientist reallocation
# at once: RI-S productivity overshoots by 36 percent. At 0.362 all eight conditions
# hold to 5.6e-13. This value needs its own justification FROM THE PATENT DATA; taken
# from the moments being fitted it would be circular.    # "backout" or "external"
const _kappa_backed = (m=_kappa_m_fit.κ, rs=_kappa_rs_fit.κ, ns=_kappa_ns_fit.κ)
const kappa_fixed = if KAPPA_MODE == "external"
    v = split(get(ENV, "KAPPA_EXT", "1.0,1.0,1.0"), ',')
    if length(v) == 3
        (m=parse(Float64, v[1]), rs=parse(Float64, v[2]), ns=parse(Float64, v[3]))
    else
        (m=_kappa_backed.m, rs=_kappa_backed.m, ns=_kappa_backed.m)
    end
else
    _kappa_backed
end
all(>(0.0), values(kappa_fixed)) || error("κ used by the model must be positive; got $kappa_fixed")

println("Patent-LP calibration of fixed κ ($(_KAPPA_START_YEAR)-$(_KAPPA_END_YEAR)):")
@printf("  KPSS 1950 shares: m=%.6f, rs=%.6f, ns=%.6f\n",
    KPSS_PATENT_SHARES_1950.m, KPSS_PATENT_SHARES_1950.rs, KPSS_PATENT_SHARES_1950.ns)
for (label, fit) in (("m", _kappa_m_fit), ("rs", _kappa_rs_fit), ("ns", _kappa_ns_fit))
    @printf("  %-2s: N index %.6f -> %.6f | Δln LP=% .6f | Δln p=% .6f | Δln A=% .6f | κ=% .6f\n",
            label, fit.N0, fit.N1, fit.dln_lp, fit.dln_p, fit.dln_A, fit.κ)
end
println()

# Population growth and normalized Lt.
# Tail rule: hold last in-sample population log growth constant for t>T_smm.
const pop_smooth_smm        = INPUT_TRENDS == "spline" ? _spline_fill(pop, df_years) :
                              INPUT_TRENDS == "fp" ? _fp2_log(pop, df_years) :
                              _endq_trend(pop, df_years; logs=true)
const gl_smm                = diff(log.(pop_smooth_smm))
const gl_tail_value         = gl_smm[end]
const gl_pop_ge             = vcat(gl_smm, fill(gl_tail_value, tail_ge))
const Lt_ge                 = _norm_level_from_loggrowth(gl_pop_ge)

# Economy-wide scientist share path ls_t.
const sci_share_tot_smm     = clamp.(INPUT_TRENDS == "spline" ? _census_path("SCI_SH") :
                                     INPUT_TRENDS == "fp" ? _fp2_trend(_col("SCI_SH"), df_years) :
                                     _endq_trend(_col("SCI_SH"), df_years), 0.001, 0.999)
const sci_share_tot_ge      = vcat(sci_share_tot_smm, fill(sci_share_tot_smm[end], tail_ge))

# RS/NS scientist share targets (in-sample, length T_smm).
# Passed directly to solve_GE; the scientist wedges are recovered to match these
# for t≤T_smm and held fixed at their 2023 values thereafter (no _ge extension).
const _sci_share_m0         = clamp.(_census_path("RDSH_GOODS"), 1e-4, 0.98)
const _sci_share_rs0        = clamp.(_census_path("RDSH_RI"),    1e-4, 0.98)
const _sci_share_scale      = min.(1.0, 0.98 ./ max.(_sci_share_m0 .+ _sci_share_rs0, 1e-6))
const sci_share_m_smm       = clamp.(_sci_share_m0 .* _sci_share_scale, 1e-4, 0.98)
const sci_share_rs_smm      = clamp.(_sci_share_rs0 .* _sci_share_scale, 1e-4, 0.98)
const sci_share_ns_smm      = clamp.(1.0 .- sci_share_m_smm .- sci_share_rs_smm, 1e-4, 0.98)

# Sanity checks for normalized paths used by solve_GE.
@assert isapprox(Lt_ge[1], 1.0; atol=1e-12)
@assert isapprox(zmt_ge[1], 1.0; atol=1e-12) && isapprox(zrst_ge[1], 1.0; atol=1e-12) && isapprox(znst_ge[1], 1.0; atol=1e-12)
@assert isapprox(hbarmt_ge[1], 1.0; atol=1e-12) && isapprox(hbarrst_ge[1], 1.0; atol=1e-12) && isapprox(hbarnst_ge[1], 1.0; atol=1e-12)

#------------------------------------------------------------------------------
# DATA MOMENTS (targeted + untargeted reporting)
#------------------------------------------------------------------------------

# Aggregate LP growth (U1) from aggregate LP raw endpoints.
const data_agg_lp_growth = _raw_ct_growth("LP_AGGREGATE", 1950, 2023)

# Ratios/shares computed from the target paths (smoothed only under SMOOTH_MODE=legacy).
const data_pm_ratio          = pm_m_tgt[end] / pm_m_tgt[1]
const data_prs_ratio         = (pm_rs_tgt[end] / pm_m_tgt[end]) / (pm_rs_tgt[1] / pm_m_tgt[1])  # U2
const data_pns_ratio         = (pm_ns_tgt[end] / pm_m_tgt[end]) / (pm_ns_tgt[1] / pm_m_tgt[1])  # U3
const data_sci_share_rs_1950 = sci_share_rs_smm[1]                          # U4
const data_sci_share_rs_2023 = sci_share_rs_smm[end]                        # U5
const data_sci_share_ns_1950 = sci_share_ns_smm[1]                          # U6
const data_sci_share_ns_2023 = sci_share_ns_smm[end]                        # U7

# Targeted data moments — computed from RAW data where intended.
const data_moments = [
    _raw_ct_growth("LP_GOODS", 1950, 2023),   # 1:  Goods LP growth (cont. avg.)
    _raw_ct_growth("LP_RI",    1950, 2023),   # 2:  RI-S LP growth
    _raw_ct_growth("LP_NRI",   1950, 2023),   # 3:  NRI-S LP growth
    _getval(1950, "VASHARE_GOODS"),            # 4:  Goods VA share, 1950  [reported only]
    _getval(1950, "VASHARE_RI"),               # 5:  RI-S VA share, 1950   [reported only]
    _getval(2023, "VASHARE_GOODS"),            # 6:  Goods VA share, 2023
    _getval(2023, "VASHARE_RI"),               # 7:  RI-S VA share, 2023
    data_pm_ratio,                              # 8:  pm[2023]/pm[1950]
    _raw_ct_growth("LP_GOODS", 1950, 1955),   # 9:  Goods LP growth, 1950-1955
    _raw_ct_growth("LP_RI",    1950, 1955),   # 10: RI-S LP growth, 1950-1955
    _raw_ct_growth("LP_NRI",   1950, 1955),   # 11: NRI-S LP growth, 1950-1955
    data_sci_share_rs_2023,                     # 12: RI-S scientist share, 2023
    data_sci_share_ns_2023,                     # 13: NRI-S scientist share, 2023
    data_sci_share_rs_2023 - data_sci_share_rs_1950,  # 14: RI-S scientist share change (2023-1950)
    data_sci_share_rs_1950,                     # 15: RI-S scientist share, 1950  [now TARGETED]
    data_sci_share_ns_1950,                     # 16: NRI-S scientist share, 1950 [now TARGETED]
]
# The data anchors are 1950 VALUE-ADDED shares, whereas θ recovery inside solve_GE
# pins the 1950 EXPENDITURE shares s_j. Invert s^VA_j = ω_j s_j / Σ_i ω_i s_i at 1950.
function _gross_from_va(va_m, va_rs, ωm, ωrs, ωns)
    raw = (va_m/ωm, va_rs/ωrs, (1.0 - va_m - va_rs)/ωns)
    tot = sum(raw)
    return (raw[1]/tot, raw[2]/tot)
end
const va1950_target_data = (data_moments[4], data_moments[5])
const smt1_target_data, srst1_target_data =
    _gross_from_va(va1950_target_data..., ωmt_ge[1], ωrst_ge[1], ωnst_ge[1])

# Model value-added shares s^VA_j = ω_j s_j / Σ_i ω_i s_i, first T periods.
function _va_shares(eq; T::Int=T_smm)
    vm  = ωmt_ge[1:T]  .* eq["smt"][1:T]
    vrs = ωrst_ge[1:T] .* eq["srst"][1:T]
    vns = ωnst_ge[1:T] .* eq["snst"][1:T]
    tot = vm .+ vrs .+ vns
    return (m=vm ./ tot, rs=vrs ./ tot, ns=vns ./ tot)
end

const moment_labels = [
    "Goods LP growth (cont. avg.)",
    "RI-S LP growth (cont. avg.)",
    "NRI-S LP growth (cont. avg.)",
    "Goods VA share, 1950",
    "RI-S VA share, 1950",
    "Goods VA share, 2023",
    "RI-S VA share, 2023",
    "pm[2023]/pm[1950] (VAPRICE_GOODS, SMOOTH_MODE target path)",
    "Goods LP growth (cont. avg. 1950-1955)",
    "RI-S LP growth (cont. avg. 1950-1955)",
    "NRI-S LP growth (cont. avg. 1950-1955)",
    "RI-S scientist share, 2023",
    "NRI-S scientist share, 2023",
    "RI-S scientist share change (2023-1950)",
    "RI-S scientist share, 1950",
    "NRI-S scientist share, 1950",
]

# Over-identified objective. The three LP-growth moments (1,2,3) are matched
# internally by the concentrated zeta solve. The outer objective targets the 2023 VA
# shares (6,7), the goods price change (8), the 2023 scientist shares (12,13) and,
# now that T_rs and T_ns are fixed rather than recovered, the 1950 scientist shares
# (15,16). Seven moments against five outer parameters (etam, etars, sigma,
# gzeta_rs, gzeta_ns), so two over-identifying restrictions.
#
# The assignment of parameters to moments is unchanged: zeta_j(0) match sectoral LP
# growth, gzeta_rs and gzeta_ns shape the 2023 scientist shares, and the Fréchet
# scales carry the 1950 research allocation. The difference is that T_rs and T_ns
# are normalized to one instead of inverting the 1950 shares, so those shares are
# targeted.
# Just identified. The three LP-growth moments (1,2,3) AND the two 1950 scientist
# shares (15,16) are matched inside `_concentrate_zetas`, which solves five
# residuals in five unknowns: the three zeta_j(0) levels and the two scientist
# mobility wedges tau^s_rs, tau^s_ns. That leaves five outer moments against the
# five outer parameters. The 5x5 inner Jacobian is full rank with condition
# number about 6.8e2, unlike the closed-form T recovery which collapsed to
# rank 1: a wedge shifts the allocation while leaving zeta_j, and hence the
# idea flow, untouched, so the LP moments separate the two.
const target_moment_idx = [6, 7, 8, 12, 13]
const target_moment_set = Set(target_moment_idx)

println("Data moments:")
for (i, (label, val)) in enumerate(zip(moment_labels, data_moments))
    @printf("  [%d] %-44s = %.6f\n", i, label, val)
end
println("Untargeted moments (data):")
@printf("  [U1] %-44s = %.6f\n", "Aggregate LP growth (LP_AGGREGATE, cont. avg.)", data_agg_lp_growth)
@printf("  [U2] %-44s = %.6f\n", "(p_rs/p_m)[2023]/(p_rs/p_m)[1950]", data_prs_ratio)
@printf("  [U3] %-44s = %.6f\n", "(p_ns/p_m)[2023]/(p_ns/p_m)[1950]", data_pns_ratio)
@printf("  [U4] %-44s = %.6f\n", "RI-S scientist share, 1950", data_sci_share_rs_1950)
@printf("  [U5] %-44s = %.6f\n", "RI-S scientist share, 2023", data_sci_share_rs_2023)
@printf("  [U6] %-44s = %.6f\n", "NRI-S scientist share, 1950", data_sci_share_ns_1950)
@printf("  [U7] %-44s = %.6f\n", "NRI-S scientist share, 2023", data_sci_share_ns_2023)
println()

#------------------------------------------------------------------------------
# PARAMETER CONFIGURATION
#------------------------------------------------------------------------------

# Stone-Geary specification, 3 September 2026 (revision: zeta0 normalised).
#
#   OUTER, drawn        etam, etars, sigma      demand side, moments 6, 7, 8
#                       phi_rs, phi_ns          returns to research scale in the two
#                                               service sectors; phi_m normalised to
#                                               0.5. These minimise the distance
#                                               between model and data researcher
#                                               allocations over the whole path.
#
#   INNER, concentrated gzeta_m, gzeta_rs, gzeta_ns  -> the three 2023 patent
#                                                       endpoints
#                       Lbar_rs, Lbar_ns             -> the two 1950 scientist shares
#
#   NORMALISED          zeta0_j = 1 in every sector.
#
# Why the normalisation. With zeta0 free AND floors present there were eight objects
# against six conditions and the Jacobian nullspace had dimension 2: the floors and
# the zeta0 levels were both aimed at the 1950 cross-section, so zeta0 simply
# re-solved to undo whatever the floors did and the equilibrium never moved. Fixing
# zeta0 leaves the floors as the only instrument on the 1950 shares, which is what
# makes them bite. The patent levels are then carried by gzeta alone, one condition
# per sector rather than two.
#
#   PRE-COMPUTED        kappa_j, from the LP accounting identity at load time.
#
# Scientist wedges are OFF (tau^s_j = 1).
# phi_m restored to the drawn set. Fixing it at 0.5 was a restriction rather than a
# normalisation: phi_j are structural elasticities and nothing normalises them.
# PHI_FIXED=1 holds phi common at PHI_COMMON and drops it from the search, leaving
# etam, etars and sigma against moments 6, 7 and 8: just identified. This is safe
# only because every innovation parameter has exactly zero leverage on those three,
# so the outer and inner blocks do not compete for the same moments.
const PHI_FIXED  = get(ENV, "PHI_FIXED", "1") != "0"
const PHI_COMMON = parse(Float64, get(ENV, "PHI_COMMON", "0.50"))

const param_lower = PHI_FIXED ? Float64[0.100, 1.001, 0.001] :
                                Float64[0.100, 1.001, 0.001, 0.010, 0.010, 0.010]
const param_upper = PHI_FIXED ? Float64[0.999, 3.500, 0.999] :
                                Float64[0.999, 3.500, 0.999, 0.900, 0.900, 0.900]
const param_names = PHI_FIXED ? ["etam","etars","sigma"] :
                                ["etam","etars","sigma","phi_m","phi_rs","phi_ns"]

# (phi_m, phi_rs, phi_ns) for a drawn outer vector, whichever mode is active.
@inline _phis(p) = PHI_FIXED ? (PHI_COMMON, PHI_COMMON, PHI_COMMON) :
                               (p[4], p[5], p[6])
const ZETA0_NORM = 1.0     # zeta0_j in every sector

# Subsistence researchers at their 1950 observed headcounts. Goods is the reference
# sector and carries no floor, so the mobile pool is what the model allocates.
const SUBSISTENCE = get(ENV, "SUBSISTENCE", "0") != "0"
const _LS_1950    = sci_share_tot_ge[1] * Lt_ge[1]
const LBAR_RS = SUBSISTENCE ? sci_share_rs_smm[1] * _LS_1950 : 0.0
const LBAR_NS = SUBSISTENCE ? sci_share_ns_smm[1] * _LS_1950 : 0.0
const theta_m_init = 0.42
const theta_rs_init = 0.08

const search_range = collect(zip(param_lower, param_upper))

println("Free parameters ($(length(param_names))): $(join(param_names, ", "))")
println("Targeted moments ($(length(target_moment_idx))): $(target_moment_idx)\n")

function _soft_failure_loss(params::Vector{Float64}; base::Float64=5e3)::Float64
    # Keep failures informative (graded by location in parameter space) rather than flat 1e10.
    span = param_upper .- param_lower
    span = map(x -> x <= 0 ? 1.0 : x, span)
    center = (param_upper .+ param_lower) ./ 2
    z = (params .- center) ./ span
    return base + dot(z, z)
end

#------------------------------------------------------------------------------
# NESTED CONCENTRATION OF ZETA
#------------------------------------------------------------------------------
# The three zeta levels are concentrated out of the SMM search: for a given
# outer point (etam, etars, sigma, gzeta_rs, gzeta_ns) an inner Broyden solve
# finds zeta0_{m,rs,ns} so that the three continuous-time LP-growth rates equal
# their data targets (moments 1,2,3). This preserves the exact LP fit while the
# trend parameters gzeta_rs, gzeta_ns shape the scientist-share dynamics
# (gzeta_m=0).
const _ZETA_LB    = 1.0e-6
const _ZETA_UB    = 100.0
# Labor-productivity targets: the smoothing-spline trend of ln LP in each sector (all
# observations, the same penalty as the other targets), cumulative since 1950. The scalar
# statistic is the slope through 1950 of a cumulative log path, sum(t x_t)/sum(t^2); equal
# statistics for model and trend leave no linear drift between them.
_lp_trend(c) = (s = _spline_fill(log.(Float64.(coalesce.(df_data[!, c], NaN))), df_years;
                                 λ=SPLINE_LAMBDA, target=collect(1950:2023)); s .- s[1])
_slope0(x) = (t = collect(0.0:length(x)-1); sum(t .* x) / sum(t .^ 2))
const _LP_TREND   = (m=_lp_trend("LP_GOODS"), rs=_lp_trend("LP_RI"), ns=_lp_trend("LP_NRI"))
const _LP_TARGETS = [_slope0(_LP_TREND.m), _slope0(_LP_TREND.rs), _slope0(_LP_TREND.ns)]
# With T_j fixed the objective is no longer flat in zeta_rs and zeta_ns, so the
# warm start affects speed rather than the solution.
const _ZETA_WARM  = Ref{Vector{Float64}}([0.25, 0.37, 0.24, 1.0, 1.0])  # warm start: 3 zeta + 2 wedges

# The same statistic for a model path.
_ct_lp_ep(Y, Lp) = _slope0([log((Y[t] / Lp[t]) / (Y[1] / Lp[1])) for t in 1:T_smm])

# One GE solve at explicit zeta levels. φ is sector-neutral (phi_m_fixed). The
# science-sector research productivities trend at gzeta_rs, gzeta_ns (gzeta_m=0 as
# numeraire): ζ_j(t)=ζ0_j·exp(gzeta_j·(t-1)). The Fréchet research block uses
# θ_F=theta_fr_fixed and the scales T_j normalized to one. The
# closed-form recovery of T_rs, T_ns from the 1950 shares is switched off, so those
# two shares become targeted moments (15, 16) rather than identities.
# One GE solve at drawn parameters p, with an explicit zeta level and trend.
# kappa is the pre-computed kappa_fixed. alloc=:stonegeary, no wedges, and T
# and theta_F are inert in this block.

# An unconverged V iteration returns the last iterate, which is not an equilibrium; the
# calibration blocks and the fit export treat it as a failed solve.
_solved(eq) = eq !== nothing && haskey(eq, "converged") && eq["converged"][1] == 1.0
_solved_or_nothing(eq) = _solved(eq) ? eq : nothing

function _solve_GE_p(p::AbstractVector{<:Real}, z0::NTuple{3,Float64},
                     g::NTuple{3,Float64}, Lb::NTuple{2,Float64};
                     phis::Union{Nothing,NTuple{3,Float64}}=nothing,
                     lbar_common::Union{Nothing,Float64}=nothing,
                     kappa::NamedTuple=kappa_fixed,
                     g2zeta::NTuple{3,Float64}=G2ZETA_DEFAULT,
                     knots::Union{Nothing,NTuple{2,Vector{Float64}}}=ZETA_KNOTS_DEFAULT,
                     V_init::Union{Nothing,NTuple{3,Vector{Float64}}}=nothing,
                     v_damp::Float64=V_DAMP_SMM, iter_max_v::Int=ITER_MAX_V_SMM)
    return solve_GE(
        T_ge, chi_smm, rho_smm,
        dMM_ge, dMR_ge, dMN_ge, dNM_ge, dNR_ge, dNN_ge,
        dRM_ge, dRR_ge, dRN_ge,
        alpha_fixed, kappa,
        p[3], p[1], p[2], etans_fixed, gamma_fixed,
        (phis === nothing ? _phis(p) : phis)[1], theta_m_init, theta_rs_init,
        g[1], g[2], g[3],
        z0[1], z0[2], z0[3],
        sci_share_tot_ge,
        sci_share_rs_smm, sci_share_ns_smm,
        prs_pm_smm, pns_pm_smm,
        smt1_target_data, srst1_target_data,
        # N_INIT=shares puts the 1950 patent cross-section into N everywhere, not
        # just the spillover aggregator. N_INIT=one keeps the normalised index.
        _NINI[1], _NINI[2], _NINI[3], gl_pop_ge,
        μxmt_ge,  μxrst_ge,  μxnst_ge,
        μymt_ge,  μyrst_ge,  μynst_ge,
        zmt_ge,   zrst_ge,   znst_ge,
        hbarmt_ge, hbarrst_ge, hbarnst_ge,
        Lt_ge;
        v_damp=v_damp, iter_max_v=iter_max_v, tol_v=TOL_V_SMM,
        V_init_m=(V_init === nothing ? nothing : V_init[1]),
        V_init_rs=(V_init === nothing ? nothing : V_init[2]),
        V_init_ns=(V_init === nothing ? nothing : V_init[3]),
        θfr=theta_fr_fixed, T_m=1.0, T_rs_fixed=1.0, T_ns_fixed=1.0,
        φrs=(phis === nothing ? _phis(p) : phis)[2],
        φns=(phis === nothing ? _phis(p) : phis)[3],
        tau_s_rs=1.0, tau_s_ns=1.0,
        alloc=ALLOC_MODE,
        Lbar_m=(lbar_common === nothing ? 0.0 : lbar_common),
        Lbar_rs=(lbar_common === nothing ? Lb[1] : lbar_common),
        Lbar_ns=(lbar_common === nothing ? Lb[2] : lbar_common),
        lambda_N0=KPSS_PATENT_SHARES_1950,
        ζ_paths=(knots !== nothing ? _zeta_knot_paths(z0, g, knots) :
                 all(iszero, g2zeta) ? nothing : _zeta_paths(z0, g, g2zeta)),
    )
end

# Curvature used when a caller does not pass g2zeta. The calibration driver sets GZETA2
# ("g2_rs,g2_ns") after each services block, so the demand block, the service-kappa loop and the
# fit export solve the model at the current estimate.
const G2ZETA_DEFAULT = let v = parse.(Float64, split(get(ENV, "GZETA2", "0,0"), ","))
    length(v) == 2 || error("GZETA2 must be g2_rs,g2_ns")
    (0.0, v[1], v[2])
end

# Non-parametric services research productivity (BLOCKB_MODE=knots): ln ζ_j at the knot years,
# linear in between, and after the last knot ln ζ_j keeps the slope of the last segment (the
# same tail rule as the exponential). Goods keeps ζ0_m exp(g_m τ). ZETA_KNOTS carries the values
# between processes as "ln ζ_rs at the knots;ln ζ_ns at the knots".
const ZETA_KNOT_YEARS = [1950, 1960, 1970, 1980, 1990, 2000, 2010, 2023]
# Knots of the published calibration, used when ZETA_KNOTS is not set. BLOCKB_MODE=legacy gives
# the services constant growth instead (Table J.1).
const ZETA_KNOTS_BASELINE = "-0.7669740071351017,-1.2246486957382428,-1.7066290909628488,-1.933301989645065,-1.757400354499845,-1.422632552992822,-1.0950590952549566,-0.6978972200356192;-1.867534748813675,-2.0532632295014563,-2.3559760266359833,-2.697943305777403,-2.8519096341934174,-2.882093593761576,-2.893396867768713,-2.9139183600694882"
const ZETA_KNOTS_DEFAULT = let s = get(ENV, "ZETA_KNOTS", get(ENV, "BLOCKB_MODE", "knots") == "knots" ? ZETA_KNOTS_BASELINE : "")
    if isempty(s)
        nothing
    else
        parts = split(s, ";")
        v = Tuple(parse.(Float64, split(p, ",")) for p in parts)
        length(v) == 2 && all(length(x) == length(ZETA_KNOT_YEARS) for x in v) ||
            error("ZETA_KNOTS must hold two lists of $(length(ZETA_KNOT_YEARS)) values")
        (v[1], v[2])
    end
end
function _knot_path(v)
    K = ZETA_KNOT_YEARS; n = length(K)
    [begin
        y = 1949 + t
        k = clamp(searchsortedlast(K, y), 1, n - 1)
        exp(v[k] + (v[k+1] - v[k]) * (y - K[k]) / (K[k+1] - K[k]))
    end for t in 1:T_ge]
end
_zeta_knot_paths(z0, g, kn) = ([z0[1] * exp(g[1] * (t - 1)) for t in 1:T_ge], _knot_path(kn[1]), _knot_path(kn[2]))

# Research productivity with curvature: ln ζ_j(t) = ln ζ0_j + g_j τ + g2_j (τ/73)^2, τ = t-1, in
# sample; after 2023 ln ζ_j continues linearly at its 2023 growth rate, so the curvature does
# not compound over the tail. g2_j = 0 gives the constant-growth path.
function _zeta_paths(z0, g, c)
    τc = T_smm - 1
    path(z, gg, cc) = [begin
            τ = t - 1; s = min(τ, τc)
            z * exp(gg * τ + cc * (s / τc)^2 + (τ > τc ? 2cc * (τ - τc) / τc : 0.0))
        end for t in 1:T_ge]
    return (path(z0[1], g[1], c[1]), path(z0[2], g[2], c[2]), path(z0[3], g[3], c[3]))
end

#------------------------------------------------------------------------------
# Bounded Levenberg-Marquardt, used by calibrate_demand.jl and calibrate_zeta_blocks.jl.
# The damping falls back toward a gradient step when the Gauss-Newton step is bad, and
# the accept test never moves uphill, so the result is a deterministic function of the start.
#------------------------------------------------------------------------------
function _lm_solve(F, u0::Vector{Float64}, lb::Vector{Float64}, ub::Vector{Float64};
                   maxiter::Int=60, tol::Float64=1e-9, fd::Float64=1e-6,
                   budget::Float64=Inf, t0::Float64=time())
    u = clamp.(copy(u0), lb, ub)
    r = F(u); (r === nothing || !all(isfinite, r)) && return (u, fill(Inf, length(u0)), false)
    ssr = sum(abs2, r); n = length(u); lam = 1e-3
    for _ in 1:maxiter
        (time() - t0 > budget) && break
        sqrt(ssr) <= tol && break
        J = Matrix{Float64}(undef, length(r), n)
        okJ = true
        for j in 1:n
            h = fd * max(abs(u[j]), 1.0)
            up = copy(u); up[j] = clamp(up[j] + h, lb[j], ub[j])
            hj = up[j] - u[j]; hj == 0.0 && (hj = h)
            rp = F(up)
            if rp === nothing || !all(isfinite, rp); okJ = false; break; end
            @inbounds for i in eachindex(r); J[i, j] = (rp[i] - r[i]) / hj; end
        end
        okJ || break
        JtJ = J' * J; Jtr = J' * r
        stepped = false
        for _ in 1:14
            step = try
                (JtJ + lam * Diagonal(diag(JtJ) .+ 1e-14)) \ (-Jtr)
            catch
                nothing
            end
            if step === nothing; lam *= 10; continue; end
            ut = clamp.(u .+ step, lb, ub)
            rt = F(ut)
            if rt !== nothing && all(isfinite, rt)
                st = sum(abs2, rt)
                if st < ssr
                    u, r, ssr = ut, rt, st
                    lam = max(lam / 10, 1e-14); stepped = true; break
                end
            end
            lam *= 10
        end
        stepped || break
    end
    return (u, r, all(isfinite, r) && maximum(abs, r) <= 1e-6)
end
