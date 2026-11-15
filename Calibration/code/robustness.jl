#==========================================================================================
 An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
 ------------------------------------------------------------------------------------------
 File:    robustness.jl
 Purpose: Alternative calibrations of Appendix J (Table J.1) and Appendix H.2 (Table H.2).
          Each case is a full recalibration from the baseline values with all other external
          inputs unchanged. Table J.1 columns:
            phi = 0.25, phi = 0.75, gamma = 0.25, gamma = 0.75, gamma = 0.90,
            services research productivity with constant growth (BLOCKB_MODE=legacy),
            machine-producer markups at their means (MUX_MODE=mean).
          Further cases reported in the text of Appendix J and in Appendix H.2:
            mean markups with phi = 0.25 and 0.75; inputs with a break in 1980 (INPUT_KNOTS).
          For each case, runs run_full_calibration.jl and Exercises 1, 3 and 5 in child
          processes, then assembles the tables with the baseline column.
 ------------------------------------------------------------------------------------------
 Usage:   julia --project=. code/robustness.jl                 recalibrate and write the tables
          julia --project=. code/robustness.jl --table-only    rewrite the tables from stored runs
          julia --project=. code/robustness.jl --only phi_075  run only these cases (comma-separated)
          julia --project=. code/robustness.jl --redo          rerun cases that already finished
          A case that already finished (all its outputs present) is reused, so an interrupted
          run picks up where it stopped.
 ------------------------------------------------------------------------------------------
 Inputs:  output/calibration/ and output/counterfactuals/ (baseline)
          data files listed in paths.jl
 Outputs: output/robustness/<case>/calibration/      calibration of each case
          output/robustness/<case>/counterfactuals/  Exercises 1, 3 and 5 of each case
          output/robustness/<case>/*.log
          output/robustness/Table_J1.tex
          output/robustness/Table_H2.tex             copied from inputs_1980 (Exercise1.jl)
          output/robustness/robustness_summary.csv   every case, unrounded
==========================================================================================#

import Pkg
include(joinpath(@__DIR__, "paths.jl"))
Pkg.activate(CAL_PATHS.root)
include(joinpath(@__DIR__, "read_results.jl"))
using Printf

# At phi = 0.75 the idea-value iteration needs a smaller weight on the new iterate and more
# iterations, in the calibration and in the counterfactuals.
const PHI075 = ["PHI_COMMON" => "0.75", "V_DAMP" => "0.3", "ITER_MAX_V" => "150"]
# (folder, column head or "" when not a Table J.1 column, environment)
const CASES = [
    ("phi_025",       "\$\\phi = 0.25\$",           ["PHI_COMMON" => "0.25"]),
    ("phi_075",       "\$\\phi = 0.75\$",           PHI075),
    ("gamma_025",     "\$\\gamma = 0.25\$",         ["GAMMA_COMMON" => "0.25"]),
    ("gamma_075",     "\$\\gamma = 0.75\$",         ["GAMMA_COMMON" => "0.75"]),
    ("gamma_090",     "\$\\gamma = 0.90\$",         ["GAMMA_COMMON" => "0.90"]),
    ("services_constant", "{\\shortstack{Services \$\\zeta_j\$\\\\constant growth}}", ["BLOCKB_MODE" => "legacy"]),
    ("mean_markups",  "{\\shortstack{Mean machine\\\\markups}}", ["MUX_MODE" => "mean"]),
    ("mean_phi_025",  "",                           ["MUX_MODE" => "mean", "PHI_COMMON" => "0.25"]),
    ("mean_phi_075",  "",                           vcat(["MUX_MODE" => "mean"], PHI075)),
    ("inputs_1980",   "",                           ["INPUT_KNOTS" => "1980"]),
]
const TABLE_ONLY = "--table-only" in ARGS
const REDO = "--redo" in ARGS
# --only name1,name2  recalibrates just these cases; the tables use the stored runs of the rest
const ONLY = let i = findfirst(==("--only"), ARGS)
    i === nothing ? [c[1] for c in CASES] : split(ARGS[i + 1], ",")
end
const STEPS = (("calibration", "run_full_calibration.jl", "calibration/final_calibration.txt"),
               ("exercise1", joinpath("counterfactuals", "Exercise1.jl"), "counterfactuals/Exercise1_contributions.csv"),
               ("exercise3", joinpath("counterfactuals", "Exercise3.jl"), "counterfactuals/Exercise3_results.csv"),
               ("exercise5", joinpath("counterfactuals", "Exercise5.jl"), "counterfactuals/Exercise5_feedback.csv"))

function run_case(name, env)
    dir = joinpath(CAL_PATHS.out_robustness, name)
    cal = joinpath(dir, "calibration"); cf = joinpath(dir, "counterfactuals")
    mkpath(cal); mkpath(cf)
    env = vcat(env, ["CALIBRATION_OUTPUT_DIR" => cal, "COUNTERFACTUAL_OUTPUT_DIR" => cf,
                     "GKSwstype" => "100"])
    for (step, script, _) in STEPS
        log = joinpath(dir, "$step.log")
        println("  $name: $step -> $log")
        cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(CAL_PATHS.root) $(joinpath(@__DIR__, script))`
        open(log, "w") do io
            run(pipeline(addenv(cmd, env...), stdout=io, stderr=io))
        end
    end
end

"""Unrounded results of one calibration folder and its counterfactual outputs."""
function results(cal_dir, cf_dir)
    r = read_final_calibration(cal_dir)
    toml = read_results_toml(cal_dir)
    contrib = Dict{String,Vector{Float64}}()
    for line in Iterators.drop(eachline(joinpath(cf_dir, "Exercise1_contributions.csv")), 1)
        f = split(line, ",")
        (f[1], f[2]) == ("1950", "2023") && (contrib[f[3]] = parse.(Float64, f[4:6]))
    end
    ex3 = Dict(split(l, ",")[1] => parse(Float64, split(l, ",")[2])
               for l in Iterators.drop(eachline(joinpath(cf_dir, "Exercise3_results.csv")), 1))
    fb = parse.(Float64, split(collect(eachline(joinpath(cf_dir, "Exercise5_feedback.csv")))[2], ","))
    lp = r.moments["LP growth, Goods, 1950-2023 (raw endpoints)"]
    pat = r.moments["Goods patent stock, ln N(2023)/N(1950)"]
    (; p = r.params, innov = contrib["endogenous_innovation"],
       lp_model = 100lp[2], lp_data = 100lp[1], pat_model = pat[2], pat_data = pat[1],
       rmse = Float64(toml["zeta_fit"]["scientist_rmse"]),
       gain = ex3["E+S"] - ex3["Baseline"], fixed = ex3["Arithmetic benchmark"] - ex3["Baseline"],
       feedback = 100fb[3] / fb[1])
end

if !TABLE_ONLY
    for (name, _, env) in CASES
        name in ONLY || continue
        done = [joinpath(CAL_PATHS.out_robustness, name, s[3]) for s in STEPS]
        if !REDO && all(isfile, done)
            println("  $name: already done, reusing (pass --redo to rerun)")
            continue
        end
        run_case(name, env)
    end
end

const BASE = results(CAL_PATHS.out_calibration, CAL_PATHS.out_counterfact)
const RES = Dict(n => results(joinpath(CAL_PATHS.out_robustness, n, "calibration"),
                              joinpath(CAL_PATHS.out_robustness, n, "counterfactuals"))
                 for (n, _, _) in CASES
                 if isfile(joinpath(CAL_PATHS.out_robustness, n, STEPS[end][3])))
# Goods LP growth along the trend the calibration targets, from the baseline fit paths.
const LP_TREND_GOODS = let lines = collect(eachline(joinpath(CAL_PATHS.out_calibration, "fit_paths.csv")))
    k = findfirst(==("t_lp_goods"), split(lines[1], ","))
    100 * parse(Float64, split(lines[end], ",")[k]) / (length(lines) - 2)
end

# Summary of every case, unrounded.
mkpath(CAL_PATHS.out_robustness)
open(joinpath(CAL_PATHS.out_robustness, "robustness_summary.csv"), "w") do io
    println(io, "case,kappa_m,kappa_rs,kappa_ns,eta_m,eta_rs,sigma,innov_m,innov_rs,innov_ns,",
                "goods_lp_growth,goods_patent_stock,scientist_rmse,gain,fixed_shares_gain,feedback_pct")
    for (n, c) in vcat([("baseline", BASE)], [(n, RES[n]) for (n, _, _) in CASES if haskey(RES, n)])
        vals = (c.p["kappa_m"], c.p["kappa_rs"], c.p["kappa_ns"], c.p["eta_m"], c.p["eta_rs"], c.p["sigma"],
                c.innov..., c.lp_model, c.pat_model, c.rmse, c.gain, c.fixed, c.feedback)
        println(io, n, ",", join((@sprintf("%.8f", v) for v in vals), ","))
    end
end

# Table J.1
const JCASES = [(n, h) for (n, h, _) in CASES if !isempty(h)]
missing_cases = [n for (n, _) in JCASES if !haskey(RES, n)]
isempty(missing_cases) || error("Table J.1 needs the cases $(join(missing_cases, ", "))")
const COLS = vcat([BASE], [RES[n] for (n, _) in JCASES])
const HEADS = vcat(["Baseline"], [h for (_, h) in JCASES])

io = IOBuffer()
ncol = length(COLS)
row(label, f, fmt) = println(io, "    $label & " *
    join((Printf.format(Printf.Format(fmt), f(c)) for c in COLS), " & ") * " \\\\")
panel(t) = println(io, "    \\multicolumn{$(ncol + 1)}{l}{\\textit{$t}} \\\\")
println(io, "% Auto-generated by robustness.jl — do not edit by hand.")
println(io, "\\begin{table}[htbp]")
println(io, "  \\centering")
println(io, "  \\caption{Alternative Calibrations}")
println(io, "  \\label{tab:robustness}")
println(io, "  \\begin{threeparttable}")
println(io, "  \\resizebox{\\textwidth}{!}{%")
println(io, "  \\begin{tabular}{l *{$ncol}{S[table-format=1.3]}}")
println(io, "    \\toprule")
println(io, "     & " * join((startswith(h, "{") ? h : "{$h}" for h in HEADS), " & ") * " \\\\")
println(io, "    \\midrule")
panel("Panel A. Internally calibrated parameters")
for (sym, key) in (("\$\\kappa_m\$", "kappa_m"), ("\$\\kappa_{rs}\$", "kappa_rs"), ("\$\\kappa_{ns}\$", "kappa_ns"),
                   ("\$\\eta_m\$", "eta_m"), ("\$\\eta_{rs}\$", "eta_rs"), ("\$\\sigma\$", "sigma"))
    row(sym, c -> c.p[key], "%.3f")
end
println(io, "    \\addlinespace")
panel("Panel B. Endogenous innovation, contribution to labor productivity growth (pp a year)")
for (j, s) in enumerate(("Goods", "RI-Services", "NRI-Services"))
    row(s, c -> c.innov[j], "%.2f")
end
println(io, "    \\addlinespace")
panel("Panel C. Fit")
row(@sprintf("Goods LP growth, \\%% (trend %.2f, data %.2f)", LP_TREND_GOODS, BASE.lp_data), c -> c.lp_model, "%.2f")
row(@sprintf("Goods patent stock (data %.3f)", BASE.pat_data), c -> c.pat_model, "%.3f")
row("Researcher shares, RMSE", c -> c.rmse, "%.3f")
println(io, "    \\addlinespace")
panel("Panel D. Structural change and aggregate labor productivity growth (pp a year)")
row("Full-shutdown gain", c -> c.gain, "%.3f")
row("Fixed-shares gain", c -> c.fixed, "%.3f")
row("Productivity feedback (\\% of gain)", c -> c.feedback, "%.0f")
println(io, "    \\bottomrule")
println(io, "  \\end{tabular}}")
println(io, "  \\begin{tablenotes}[flushleft]")
println(io, "    \\footnotesize")
println(io, "    \\RaggedRight")
println(io, "    \\item[] \\textit{Notes:} Each column is a full recalibration. LP denotes labor productivity. The baseline sets")
println(io, "  \$\\phi = \\gamma = 0.5\$, research productivity in RI-Services and NRI-Services is free at the knots 1950, 1960,")
println(io, "  \\ldots, 2010, and 2023 and log-linear in between, Goods research productivity grows at a constant rate, and")
println(io, "  machine-producer markups follow quadratic trends. The seventh column gives research productivity in RI-Services")
println(io, "  and NRI-Services constant growth rates. The last column holds each machine-producer markup at its 1950--2023")
println(io, "  mean. In every column the model matches the slopes of the trends of labor productivity in RI-Services and")
println(io, "  NRI-Services exactly, and Goods labor productivity growth is the model growth between 1950 and 2023. The")
println(io, "  Goods patent stock row reports \$\\ln N_m(2023)/N_m(1950)\$. The researcher RMSE is computed over the sectoral")
println(io, "  shares of RI-Services and NRI-Services at the census years 1950--2020 and 2023. Panel D reports the difference")
println(io, "  between the full shutdown of Section 5.3 and the baseline, the corresponding difference for the fixed-shares")
println(io, "  benchmark, and the share of the full-shutdown gain due to the productivity feedback of Section 5.5. With")
println(io, "  \$\\phi = 0.75\$ the counterfactual equilibria are solved with stronger damping of the value iteration.")
println(io, "  \\end{tablenotes}")
println(io, "  \\end{threeparttable}")
println(io, "\\end{table}")

const OUT = joinpath(CAL_PATHS.out_robustness, "Table_J1.tex")
open(f -> write(f, String(take!(io))), OUT, "w")
println("Wrote ", OUT)

# Table H.2: written by Exercise1.jl in the inputs_1980 case.
const H2 = joinpath(CAL_PATHS.out_robustness, "inputs_1980", "counterfactuals", "Table_H2.tex")
if isfile(H2)
    cp(H2, joinpath(CAL_PATHS.out_robustness, "Table_H2.tex"); force=true)
    println("Wrote ", joinpath(CAL_PATHS.out_robustness, "Table_H2.tex"))
end
