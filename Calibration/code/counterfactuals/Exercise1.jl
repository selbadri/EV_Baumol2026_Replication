#==========================================================================================
 An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
 ------------------------------------------------------------------------------------------
 File:    Exercise1.jl
 Purpose: Exercise 1. Exact decomposition of sectoral labor productivity growth into five
          channels along the baseline path (Section 5.1, Appendix H.1).
 ------------------------------------------------------------------------------------------
 Usage:   julia --project=. code/counterfactuals/Exercise1.jl
 ------------------------------------------------------------------------------------------
 Inputs:  output/calibration/final_calibration.txt, calibration_results.toml (via ../baseline.jl)
 Outputs: output/counterfactuals/Table_3.tex        1950-2023
          output/counterfactuals/Table_H1.tex       1950-1980 and 1980-2023
          output/counterfactuals/Table_H2.tex       all three periods, only with INPUT_KNOTS set
                                                    (the Appendix H.2 calibration, robustness.jl)
          output/counterfactuals/Exercise1_contributions.csv   unrounded contributions (read by
                                                    robustness.jl for Tables J.1 and H.2)
          output/counterfactuals/Exercise1_channels.pdf   cumulative channel contributions (not in the paper)
 ------------------------------------------------------------------------------------------
 Notes:
 Exercise 1: Exact decomposition of sectoral labor productivity growth
 g_{A_j}(t) along the baseline equilibrium path.

 Anchor:
     A_j(t) = Y_j/L^p_j
            = N_j^{e_Y,j} · [μ^y_j μ^x_j (1-α) p_j z_j]^{(1-α)/α} · h̄_j

 The final markups are constant and the machine markups follow quadratic trends, so
 annual log differences give an EXACT decomposition into five channels (Equation 18):

   (1) Endogenous innovation         e_Y,j · Δln N_j
   (2) Embodied technical progress   (1-α)/α · Δln z_j
   (3) Relative price                (1-α)/α · Δln p_j
   (4) Labor efficiency              Δln h̄_j
   (5) Machine market power          (1-α)/α · Δln μ^x_j + Δe_Y,j · ln N_j(t-1)

 The fifth channel is zero when the machine markups are held at their means
 (MUX_MODE=mean, Table J.1); the tables then leave it out.
==========================================================================================#

# Load the baseline (exposes eq_base, baseline_params, baseline_exog).
include(joinpath(@__DIR__, "..", "baseline.jl"))

using Printf
using Statistics
import Plots

#------------------------------------------------------------------------------
# Configuration
#------------------------------------------------------------------------------
const YEAR0 = 1950                       # base year of the solved sample (t index 1)
const ND    = baseline_params.T_data     # 74 data years, 1950–2023
const ALPHA = baseline_params.alpha

# Reporting periods (start year, end year).
const PERIODS = [(1950, 2023), (1950, 1980), (1980, 2023)]

# Model-simulated LP-growth benchmark for validation check 2. Computed
# from the baseline solver output exactly as the calibration builds its SIMULATED
# sectoral LP moment (smm_ge_unified.jl `ann_lp`): the endpoint CAGR of Y_j/L^p_j.
# This keeps Exercise 1 purely model-based — no empirical data enters.
ann_lp_model(Y, Lp) = ((Y[ND] / Lp[ND]) / (Y[1] / Lp[1]))^(1.0 / (ND - 1)) - 1.0

const OUTDIR = CAL_PATHS.out_counterfact

yr_to_idx(y::Int) = y - YEAR0 + 1        # 1950 -> 1

#------------------------------------------------------------------------------
# 1. Assemble per-sector input series (sliced to the data window 1:ND)
#------------------------------------------------------------------------------
struct SectorInputs
    name :: String
    N    :: Vector{Float64}   # normalized variety stock N_j(1950)=1 (solver)
    p    :: Vector{Float64}   # sectoral price p_j           (solver)
    Y    :: Vector{Float64}   # gross output Y_j             (solver)
    Lp   :: Vector{Float64}   # production labor L^p_j       (solver)
    eY   :: Vector{Float64}   # love-for-variety exponent    (solver)
    z    :: Vector{Float64}   # machine productivity z_j     (exogenous)
    hbar :: Vector{Float64}   # labor efficiency h̄_j         (exogenous)
    mux  :: Vector{Float64}   # inverse machine markup μ^x_j (exogenous)
end

_slice(v) = collect(v[1:ND])

const SECTORS = [
    SectorInputs("Goods",
        _slice(eq_base["Nmt"]),  _slice(eq_base["pmt"]),  _slice(eq_base["Ymt"]),
        _slice(eq_base["Lpmt"]), _slice(eq_base["e_Y_mt"]),
        _slice(baseline_exog.zm),   _slice(baseline_exog.hbarm), _slice(baseline_exog.muxm)),
    SectorInputs("RI-Serv.",
        _slice(eq_base["Nrst"]), _slice(eq_base["prst"]), _slice(eq_base["Yrst"]),
        _slice(eq_base["Lprst"]), _slice(eq_base["e_Y_rst"]),
        _slice(baseline_exog.zrs),   _slice(baseline_exog.hbarrs), _slice(baseline_exog.muxrs)),
    SectorInputs("NRI-Serv.",
        _slice(eq_base["Nnst"]), _slice(eq_base["pnst"]), _slice(eq_base["Ynst"]),
        _slice(eq_base["Lpnst"]), _slice(eq_base["e_Y_nst"]),
        _slice(baseline_exog.zns),   _slice(baseline_exog.hbarns), _slice(baseline_exog.muxns)),
]

#------------------------------------------------------------------------------
# 2. Core decomposition
#    Returns per-year contribution vectors (length ND). Index t holds the
#    contribution to the change from year t-1 to year t; entry 1 is 0.
#------------------------------------------------------------------------------
struct Decomp
    row1  :: Vector{Float64}
    row2  :: Vector{Float64}
    row3  :: Vector{Float64}
    row4  :: Vector{Float64}
    row5  :: Vector{Float64}
    total :: Vector{Float64}          # sum of rows (1)-(5)
    gA    :: Vector{Float64}          # direct Δln(Y/L^p) from solver
    resid :: Vector{Float64}          # total - gA (validation 1)
end

function decompose(s::SectorInputs)
    n = ND
    c = (1.0 - ALPHA) / ALPHA
    r1 = zeros(n); r2 = zeros(n); r3 = zeros(n); r4 = zeros(n); r5 = zeros(n)
    tot = zeros(n); gA = zeros(n); res = zeros(n)

    lnN = log.(s.N)
    for t in 2:n
        # (1) Endogenous innovation at the current elasticity e_Y,j(t)
        r1[t]  = s.eY[t] * (lnN[t] - lnN[t-1])

        # (2) Embodied technical progress
        r2[t]  = c * (log(s.z[t]) - log(s.z[t-1]))

        # (3) Relative price effect
        r3[t]  = c * (log(s.p[t]) - log(s.p[t-1]))

        # (4) Labor efficiency
        r4[t]  = log(s.hbar[t]) - log(s.hbar[t-1])

        # (5) Machine market power. e_Y,j = (κ_j + (1-μ^x_j)/μ^x_j)(1-α)/α moves with μ^x_j, so
        # the drift of the elasticity on the inherited stock belongs to this channel too.
        r5[t]  = c * (log(s.mux[t]) - log(s.mux[t-1])) + (s.eY[t] - s.eY[t-1]) * lnN[t-1]

        tot[t] = r1[t] + r2[t] + r3[t] + r4[t] + r5[t]
        gA[t]  = log(s.Y[t] / s.Lp[t]) - log(s.Y[t-1] / s.Lp[t-1])
        res[t] = tot[t] - gA[t]
    end
    return Decomp(r1, r2, r3, r4, r5, tot, gA, res)
end

const DECOMPS = [decompose(s) for s in SECTORS]
# The machine-market-power row is reported only when the markups vary.
const MARKUPS_VARY = any(any(!iszero, d.row5) for d in DECOMPS)

#------------------------------------------------------------------------------
# 3. Validation
#------------------------------------------------------------------------------
function run_validation()
    println("\n" * "="^72)
    println("  Exercise 1 — validation log")
    println("="^72)

    # Exercise 1 and all LP accounting use normalized N, never the Lambda-only
    # share-seeded stocks N^Λ.
    all(isapprox(s.N[1], 1.0; atol=1e-12, rtol=1e-12) for s in SECTORS) ||
        error("Exercise 1 requires normalized N_j(1950)=1")

    # Check 1: rows reproduce the direct log-difference of Y/L^p.
    println("\n[Check 1] max |Σ rows − Δln(Y/L^p)| over 1950–2023 (should be ~1e-12):")
    for (s, d) in zip(SECTORS, DECOMPS)
        mx = maximum(abs.(d.resid[2:end]))
        flag = mx < 1e-10 ? "OK" : "FAIL"
        @printf("  %-28s  max residual = %.3e   [%s]\n", s.name, mx, flag)
        mx < 1e-10 || error("Exercise 1 identity failed for $(s.name): residual=$mx")
    end

    # Check 2: the decomposition's g_A series reproduces the model-simulated LP
    # growth used in the calibration (model-vs-model consistency, §1.8). Both are
    # the endpoint CAGR of Y_j/L^p_j from the same Run D solve, so the gap is 0
    # to numerical precision.
    println("\n[Check 2] full-sample sectoral LP growth (annualized, model-simulated):")
    @printf("  %-28s  %10s  %10s  %10s\n", "sector", "decomp gA", "solver LP", "gap(pp)")
    for (s, d) in zip(SECTORS, DECOMPS)
        decomp_cagr = exp(mean(d.gA[2:end])) - 1.0   # CAGR implied by decomposed gA
        solver_cagr = ann_lp_model(s.Y, s.Lp)         # calibration's simulated moment
        @printf("  %-28s  %9.4f%%  %9.4f%%  %+8.2e\n",
                s.name, 100decomp_cagr, 100solver_cagr, 100 * (decomp_cagr - solver_cagr))
        abs(decomp_cagr - solver_cagr) < 1e-12 ||
            error("Exercise 1 CAGR check failed for $(s.name)")
    end
    println()
end

#------------------------------------------------------------------------------
# 4. Period aggregation + console table
#    Annualized percentage-point contribution over [y0, y1]:
#        100 · Σ_{t=i0+1}^{i1} vec[t] / (i1 - i0)
#------------------------------------------------------------------------------
function agg(vec::Vector{Float64}, y0::Int, y1::Int)
    i0 = yr_to_idx(y0); i1 = yr_to_idx(y1)
    return 100.0 * sum(vec[(i0 + 1):i1]) / (i1 - i0)
end

# (label, indent-level, field-symbol) row template shared by console & LaTeX.
const ROWS = [
    ("(1) Endogenous innovation",       0, :row1),
    ("(2) Embodied technical progress", 0, :row2),
    ("(3) Relative price effect",       0, :row3),
    ("(4) Labor efficiency",            0, :row4),
    (MARKUPS_VARY ? [("(5) Machine market power", 0, :row5)] : [])...,
    ("Total",                           0, :total),
]

function print_console_table()
    println("="^72)
    println("  Exercise 1 — decomposition of sectoral labor productivity growth")
    println("  (annualized contributions, percentage points)")
    println("="^72)
    hdrs = [@sprintf("%d–%d", y0, y1) for (y0, y1) in PERIODS]
    for (s, d) in zip(SECTORS, DECOMPS)
        println("\n  ── $(s.name) " * "─"^max(0, 52 - length(s.name)))
        @printf("  %-34s %10s %10s %10s\n", "", hdrs[1], hdrs[2], hdrs[3])
        for (label, indent, sym) in ROWS
            vec = getfield(d, sym)
            vals = [agg(vec, y0, y1) for (y0, y1) in PERIODS]
            pad = "  " * "  "^indent
            name = pad * label
            @printf("  %-34s %10.2f %10.2f %10.2f\n", name, vals[1], vals[2], vals[3])
        end
    end
    println()
end

#------------------------------------------------------------------------------
# 5. LaTeX booktabs table
#------------------------------------------------------------------------------
"""Unrounded contributions by period and channel (read by robustness.jl for Table J.1)."""
function write_contributions_csv(path::String)
    open(path, "w") do io
        println(io, "start,end,channel,goods,ri_services,nri_services")
        for (y0, y1) in PERIODS, (name, sym) in (("endogenous_innovation", :row1),
                ("embodied_technical_progress", :row2), ("relative_price", :row3),
                ("labor_efficiency", :row4), ("machine_market_power", :row5), ("total", :total))
            vals = [agg(getfield(d, sym), y0, y1) for d in DECOMPS]
            println(io, "$y0,$y1,$name," * join((@sprintf("%.10f", v) for v in vals), ","))
        end
    end
    println("  Wrote contributions -> $path")
end

const MECH_BLOCKS = [("I.\\ Endogenous innovation", :row1),
                     ("II.\\ Embodied technical progress", :row2),
                     ("III.\\ Relative price effect", :row3),
                     ("IV.\\ Labor efficiency", :row4),
                     (MARKUPS_VARY ? [("V.\\ Machine market power", :row5)] : [])...]

"""Printed entries for one period: channels rounded so that they add to the printed total,
sector by sector. Returns a matrix (channel or total) x (sector)."""
function printed_panel(y0, y1)
    tot = [round(agg(d.total, y0, y1), digits=2) for d in DECOMPS]
    out = zeros(length(MECH_BLOCKS) + 1, 3)
    for j in 1:3
        out[1:end-1, j] = round_to_total([agg(getfield(DECOMPS[j], b[2]), y0, y1) for b in MECH_BLOCKS], tot[j])
        out[end, j] = tot[j]
    end
    out
end

"""Write one booktabs table for the given periods: Table 3 has one period, Table H.1 two
periods side by side, and Table H.2 (stacked=true) one panel per period, stacked."""
function write_latex_table(path::String, periods, caption, label, note; stacked::Bool=false)
    br = "\\\\"
    fmt2(x) = fmt_fixed(x, 2)
    np = stacked ? 1 : length(periods)
    panels = [printed_panel(y0, y1) for (y0, y1) in periods]
    io = IOBuffer()
    println(io, "\\begin{table}[htbp]")
    println(io, "  \\centering")
    println(io, "  \\caption{$caption}")
    println(io, "  \\label{$label}")
    println(io, "  \\begin{threeparttable}")
    println(io, "  \\centering")
    println(io, "  \\begin{tabular}{l *{$(3np)}{S[table-format=-1.2]}}")
    println(io, "    \\toprule")
    if np == 1
        println(io, "     & \\multicolumn{3}{c}{Sector} " * br)
        println(io, "    \\cmidrule(lr){2-4}")
    else
        heads = [@sprintf("\\multicolumn{3}{c}{Panel %s. %d--%d}", 'A' + k - 1, y0, y1)
                 for (k, (y0, y1)) in enumerate(periods)]
        println(io, "     & " * join(heads, " & ") * " " * br)
        println(io, "    " * join([@sprintf("\\cmidrule(lr){%d-%d}", 3k - 1, 3k + 1) for k in 1:np], " "))
    end
    println(io, "     & " * join(repeat(["{Goods}", "{RI-Serv.}", "{NRI-Serv.}"], np), " & ") * " " * br)
    println(io, "    \\midrule")
    rows(ps) = begin
        for (i, (lab, _)) in enumerate(MECH_BLOCKS)
            println(io, "    $lab & " * join([fmt2(p[i, j]) for p in ps for j in 1:3], " & ") * " " * br)
        end
        println(io, "    \\midrule")
        println(io, "    Total & " * join([fmt2(p[end, j]) for p in ps for j in 1:3], " & ") * " " * br)
    end
    if stacked
        for (k, (y0, y1)) in enumerate(periods)
            k > 1 && println(io, "    \\addlinespace")
            println(io, @sprintf("    \\multicolumn{4}{l}{\\textit{Panel %s. %d--%d}} ", 'A' + k - 1, y0, y1) * br)
            rows([panels[k]])
        end
    else
        rows(panels)
    end
    println(io, "    \\bottomrule")
    println(io, "  \\end{tabular}")
    println(io, "  \\begin{tablenotes}[flushleft]")
    println(io, "    \\footnotesize")
    println(io, "    \\RaggedRight")
    println(io, "    \\item[] \\textit{Notes:} $note")
    println(io, "  \\end{tablenotes}")
    println(io, "  \\end{threeparttable}")
    println(io, "\\end{table}")
    open(f -> write(f, String(take!(io))), path, "w")
    println("  Wrote LaTeX table -> $path")
end

const NOTE_3 = "Entries report annualized contributions, in percentage points, to the growth of sectoral labor productivity \$\\mathcal{A}_j = Y_j/L_j^p\$ over 1950--2023. The channels follow Equation (18): endogenous innovation is \$e_j^Y(t) g_{N_j}(t)\$, embodied technical progress reflects machine productivity, the relative price effect reflects sectoral output prices" *
    (MARKUPS_VARY ? ", and machine market power collects the change in the machine-producer markup and the change it induces in \$e_j^Y(t)\$." : ".") *
    " Entries are rounded so that they add to the printed totals."
const NOTE_H2 = "Entries report annualized contributions, in percentage points, to sectoral labor productivity growth in a recalibration in which labor efficiency and machine productivity grow at constant rates between their " *
    join(string.(vcat(1950, INPUT_KNOTS)), ", ") * ", and 2023 observations. The decomposition follows Table 3. Entries are rounded so that they add to the printed totals."
const NOTE_H1 = "Entries report annualized contributions, in percentage points, to sectoral labor productivity growth over 1950--1980 (Panel A) and 1980--2023 (Panel B), where sectoral labor productivity is \$\\mathcal{A}_j = Y_j/L_j^p\$. The decomposition follows Table 3. Entries are rounded so that they add to the printed totals."

#------------------------------------------------------------------------------
# 6. Companion figure
#    Cumulative contribution of each channel to ln A_j(t),
#    normalized to zero in the initial year, one panel per sector.
#------------------------------------------------------------------------------
function write_figure(path::String)
    years = collect(YEAR0:(YEAR0 + ND - 1))          # 1950 .. 2023
    channels = [(:row1, "Love-for-variety", :steelblue),
                (:row2, "Embodied tech.",   :seagreen),
                (:row3, "Relative price",    :indianred),
                (:row4, "Labor efficiency",  :goldenrod),
                (:row5, "Machine market power", :purple)]

    panels = Plots.Plot[]
    for (sidx, (s, d)) in enumerate(zip(SECTORS, DECOMPS))
        plt = Plots.plot(; title = s.name, titlefontsize = 9,
                         legend = (sidx == 1 ? :topleft : false),
                         legendfontsize = 6,
                         xlabel = "year",
                         ylabel = (sidx == 1 ? "cumulative contribution to ln \$\\mathcal{A}_j\$" : ""))
        for (sym, lab, col) in channels
            cum = cumsum(getfield(d, sym))           # per-year contributions; entry 1 = 0
            Plots.plot!(plt, years, cum; label = lab, color = col, lw = 1.6)
        end
        Plots.hline!(plt, [0.0]; color = :black, lw = 0.7, label = false)
        push!(panels, plt)
    end
    fig = Plots.plot(panels...; layout = (1, 3), size = (1300, 380),
                     left_margin = 6Plots.mm, bottom_margin = 6Plots.mm)
    Plots.savefig(fig, path)
    println("  Wrote figure -> $path")
end

#------------------------------------------------------------------------------
# 7. Driver
#------------------------------------------------------------------------------
function main()
    isdir(OUTDIR) || mkpath(OUTDIR)
    run_validation()
    print_console_table()
    write_latex_table(joinpath(OUTDIR, "Table_3.tex"), [(1950, 2023)],
        "Decomposition of Sectoral Labor Productivity Growth, 1950--2023", "tab:ex1_decomp", NOTE_3)
    write_latex_table(joinpath(OUTDIR, "Table_H1.tex"), [(1950, 1980), (1980, 2023)],
        "Decomposition of Sectoral Labor Productivity Growth by Subperiod", "tab:ex1_decomp_sub", NOTE_H1)
    isempty(INPUT_KNOTS) || write_latex_table(joinpath(OUTDIR, "Table_H2.tex"), PERIODS,
        "Decomposition of Sectoral Labor Productivity Growth, Inputs with a Break in " * join(INPUT_KNOTS, ", "),
        "tab:ex1_decomp_break", NOTE_H2; stacked=true)
    write_contributions_csv(joinpath(OUTDIR, "Exercise1_contributions.csv"))
    write_figure(joinpath(OUTDIR, "Exercise1_channels.pdf"))
    println("\nExercise 1 complete.")
end

main()
