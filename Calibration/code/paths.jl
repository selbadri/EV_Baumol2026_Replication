#==========================================================================================
 An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
 ------------------------------------------------------------------------------------------
 File:    paths.jl
 Purpose: Single place for every file path used by the calibration code. All paths are
          relative to this folder, so nothing needs editing on a new machine.
 ------------------------------------------------------------------------------------------
 Inputs:  produced by the sibling folders of the replication package
   ../Empirical/output/dataset/data_moments.xlsx        data moments (sheet DataMoments)
   ../Empirical/input/manual/NSF/nsf_rd_benchmarks.csv  NSF R&D benchmarks (Figure B1)
   ../Cross_Citations/output/citations/spillover_backward1.csv   citation-based spillover weights
   ../Cross_Citations/output/robustness/stocks_summary.csv       KPSS patent-stock summary
   ../Cross_Citations/output/patents/class_stocks.csv            annual patent stocks by aggregation
                                                          (only for non-default PATENT_SERIES)
 Outputs:
   output/calibration/       calibrated parameters, fit report, fitted paths
   output/figures/           motivation figures and calibration figures
   output/counterfactuals/   tables and figures of the quantitative exercises
   output/robustness/        alternative calibrations (Tables J.1, H.2)
   output/logs/              logs written by run_all.jl
==========================================================================================#

if !isdefined(@__MODULE__, :CAL_PATHS)
    const CAL_PATHS = let
        root      = normpath(joinpath(@__DIR__, ".."))
        package   = normpath(joinpath(root, ".."))
        empirical = joinpath(package, "Empirical")
        crosscit  = joinpath(package, "Cross_Citations")
        (
            root            = root,
            code            = joinpath(root, "code"),
            # Inputs
            data_moments    = get(ENV, "DATA_MOMENTS",
                                  joinpath(empirical, "output", "dataset", "data_moments.xlsx")),
            nsf_benchmarks  = joinpath(empirical, "input", "manual", "NSF", "nsf_rd_benchmarks.csv"),
            spillovers      = joinpath(crosscit, "output", "citations", "spillover_backward1.csv"),
            stocks_summary  = joinpath(crosscit, "output", "robustness", "stocks_summary.csv"),
            class_stocks    = joinpath(crosscit, "output", "patents", "class_stocks.csv"),
            # Outputs
            # The robustness runs (robustness.jl) redirect these two with ENV variables.
            out_calibration = get(ENV, "CALIBRATION_OUTPUT_DIR",
                                  joinpath(root, "output", "calibration")),
            out_figures     = joinpath(root, "output", "figures"),
            out_counterfact = get(ENV, "COUNTERFACTUAL_OUTPUT_DIR",
                                  joinpath(root, "output", "counterfactuals")),
            out_robustness  = joinpath(root, "output", "robustness"),
            out_logs        = joinpath(root, "output", "logs"),
        )
    end
end
