#==========================================================================================
 An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
 ------------------------------------------------------------------------------------------
 File:    run_all.jl
 Purpose: Run the whole Calibration step in order:
            1. motivation_plots.jl       data figures (Figures 1-5, B.1)
            2. run_full_calibration.jl   internal calibration
               calibration_table.jl      Table 2
            3. calibration_plots.jl      calibration figures (Figures 6-7)
            4. counterfactuals/          quantitative exercises (Tables 3-8, H.1, H.3, H.4)
            5. robustness.jl             alternative calibrations (Tables J.1, H.2)
          Each step runs in its own Julia process, because the model defines
          configuration-dependent constants when it is loaded.
 ------------------------------------------------------------------------------------------
 Usage:   julia --project=. run_all.jl                     everything
          julia --project=. run_all.jl --skip-calibration  keep output/calibration as is and
                                                           redo figures and counterfactuals
          julia --project=. run_all.jl --skip-robustness   leave out step 5 (ten recalibrations)
          The two flags can be combined.
 ------------------------------------------------------------------------------------------
 Inputs:  ../Empirical/output/dataset/data_moments.xlsx and the other files in code/paths.jl
 Outputs: output/figures/, output/calibration/, output/counterfactuals/, output/robustness/
          output/logs/<step>.log   one log per step
==========================================================================================#

import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

include(joinpath(@__DIR__, "code", "paths.jl"))
using Dates

const SKIP_CALIBRATION = "--skip-calibration" in ARGS
const SKIP_ROBUSTNESS  = "--skip-robustness" in ARGS

const STEPS = [
    ("motivation_plots",      "motivation_plots.jl"),
    ("run_full_calibration",  "run_full_calibration.jl"),
    ("calibration_table",     "calibration_table.jl"),
    ("calibration_plots",     "calibration_plots.jl"),
    ("Exercise1",             joinpath("counterfactuals", "Exercise1.jl")),
    ("Exercise2",             joinpath("counterfactuals", "Exercise2.jl")),
    ("Exercise3",             joinpath("counterfactuals", "Exercise3.jl")),
    ("Exercise4",             joinpath("counterfactuals", "Exercise4.jl")),
    ("Exercise4b",            joinpath("counterfactuals", "Exercise4b.jl")),
    ("Exercise5",             joinpath("counterfactuals", "Exercise5.jl")),
    ("robustness",            "robustness.jl"),
]

for p in (CAL_PATHS.out_calibration, CAL_PATHS.out_figures,
          CAL_PATHS.out_counterfact, CAL_PATHS.out_logs)
    mkpath(p)
end
for f in (CAL_PATHS.data_moments, CAL_PATHS.nsf_benchmarks, CAL_PATHS.spillovers,
          CAL_PATHS.stocks_summary)
    isfile(f) || error("Missing input $f. Run the Empirical and Cross_Citations steps first.")
end

const T0 = now()
for (name, script) in STEPS
    if SKIP_CALIBRATION && name == "run_full_calibration"
        println("[skip] $name (using output/calibration as stored)")
        continue
    end
    if SKIP_ROBUSTNESS && name == "robustness"
        println("[skip] $name")
        continue
    end
    log_path = joinpath(CAL_PATHS.out_logs, "$name.log")
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(CAL_PATHS.root) $(joinpath(CAL_PATHS.code, script))`
    t = now()
    println("[$(Dates.format(t, "HH:MM:SS"))] $name -> $log_path")
    open(log_path, "w") do io
        run(pipeline(addenv(cmd, "GKSwstype" => "100"), stdout=io, stderr=io))
    end
    println("           done in $(round(now() - t, Second))")
end
println("All steps finished in $(round(now() - T0, Minute)).")
