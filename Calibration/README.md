# An Anatomy of the Baumol Cost Disease: Calibration and Quantitative Exercises

**Selim Elbadri & Juan I. Vizcaino (2026)**

## Description

This folder holds the model side of *An Anatomy of the Baumol Cost Disease*. A set of Julia scripts:

1. draws the data figures of Section 2 and Appendix B;
2. calibrates the three-sector endogenous-growth model with structural change to US data for 1950–2023;
3. draws the calibration figures;
4. runs the quantitative exercises of Section 5 and Appendix H;
5. runs the alternative calibrations of Appendices H.2 and J.

The three sectors are Goods, research-intensive services (RI-Services) and non-research-intensive services (NRI-Services), defined as in `../Empirical`.

The data come from the two sibling folders of the replication package, and nothing is copied here:

| Input | Built by |
|-------|----------|
| `../Empirical/output/dataset/data_moments.xlsx` | `Empirical/master.do` |
| `../Empirical/input/manual/NSF/nsf_rd_benchmarks.csv` | Hand-collected NSF benchmarks (see the README there) |
| `../Cross_Citations/output/citations/spillover_backward1.csv` | Cross-sector citation probabilities |
| `../Cross_Citations/output/robustness/stocks_summary.csv` | KPSS patent-stock summary |
| `../Cross_Citations/output/patents/class_stocks.csv` | Annual patent stocks. Read only when `PATENT_SERIES` is set to something other than the default |

All paths are set in one file, `code/paths.jl`, relative to this folder.

## Attribution & License

If you use this code, please cite the paper:

> Elbadri, S. & Vizcaino, J. I. (2026). *An Anatomy of the Baumol Cost Disease*. Working paper.

The code in this folder is released under the MIT License (see `LICENSE`).

## Repository Structure

```
Calibration/
├── run_all.jl                 Runs every step below, in order
├── Project.toml, Manifest.toml   Julia environment (Julia 1.12)
├── README.md
├── LICENSE
├── code/
│   ├── paths.jl               All input and output paths
│   ├── motivation_plots.jl    Step 1: data figures
│   ├── run_full_calibration.jl        Step 2: calibration driver (outer loop)
│   │   ├── calibrate_service_kappas.jl    kappa_rs, kappa_ns block
│   │   ├── calibrate_zeta_blocks.jl       research-productivity block (zeta, kappa_m)
│   │   ├── calibrate_demand.jl            demand block (eta_m, eta_rs, sigma) and ln D
│   │   └── dump_fit_paths.jl              final report and fitted paths
│   ├── calibration_table.jl   Table 2 from the stored calibration
│   ├── calibration_plots.jl   Step 3: calibration figures
│   ├── baseline.jl            Loads the model at the calibrated parameters (steps 3–5)
│   ├── counterfactuals/       Step 4: Exercise1.jl ... Exercise5.jl
│   ├── robustness.jl          Step 5: alternative calibrations (Tables J.1, H.2)
│   ├── read_results.jl        Reads a stored calibration (used by the two table scripts)
│   ├── smm_ge_unified.jl      Model set-up: data, inputs, targets, constants, GE wrapper
│   ├── equilibrium_path.jl    GE solver: 100-period transition path
│   └── v_seed.jl              Warm start for the patent values
└── output/
    ├── figures/               Figures 1–7, B.1
    ├── calibration/           Parameters, fit report, fitted paths, Table 2
    ├── counterfactuals/       Tables 3–8, H.1, H.3, H.4
    ├── robustness/            One folder per alternative calibration, Tables J.1 and H.2, robustness_summary.csv
    └── logs/                  One log per step (not tracked)
```

## How the calibration works

**Inputs.** The exogenous inputs are built in `smm_ge_unified.jl` (Section 4.1 of the paper):

- Machine-producer markups `1/mu^x_j` follow quadratic time trends fitted to the measured profit shares of produced nonresidential capital, 1950–2023 (`MUX_MODE=quad`; `MUX_MODE=mean` holds each at its 1950–2023 mean). The final-producer markups `1/mu^y_j` are constant at the means implied by the labor shares.
- Labor efficiency `hbar_j` and machine productivity `z_j = 1/(p^K_j D mu^x_j)` grow at constant rates between their 1950 and 2023 observations (`INPUT_KNOTS=1980` adds a break in 1980, Appendix H.2).
- `D(t)` is the model's chained value-added deflator in final-good units. The data deflate capital prices by the private VA deflator (a Törnqvist index), while machines in the model cost final goods; `D` converts one into the other. It is an equilibrium object, updated at every sweep of the calibration.

**Parameters.** The model has 26 internally calibrated parameters (Table 2):

| Block | Parameters | Main targets |
|-------|------------|--------------|
| Goods innovation | `kappa_m`, `zeta_m(1950)`, `gzeta_m` | Goods labor-productivity trend path and Goods patent-stock path, 1950–2023 |
| Services innovation | `kappa_rs`, `kappa_ns` | Slopes of the RI- and NRI-Services labor-productivity trends |
| | `zeta_rs(t_k)`, `zeta_ns(t_k)` at the knots 1950, 1960, …, 2010, 2023 | RI- and NRI-Services R&D-worker share paths |
| Demand | `eta_m`, `eta_rs`, `sigma` | Goods and RI-Services value-added share paths |

The consumption weights `theta_j` are recovered from the 1950 value-added shares and `eta_ns = 1` is a normalization. Goods research productivity grows at a constant rate; in the services `ln zeta_j` is linear between the knots (`BLOCKB_MODE=knots`; `BLOCKB_MODE=legacy` gives constant growth, Appendix J). The labor-productivity targets are the smoothing-spline trends of log labor productivity; for the services the statistic is the slope through the origin of cumulative log growth on years since 1950. Relative service prices are imposed through production wedges, and model prices are measured against a Fisher chain of the three sector prices with value-added weights, the counterpart of the data deflator (on the data the two differ by 1e-5 log points over 1950–2023), so the Goods price is a check rather than a target. Externally set parameters (Table 1) are constants in `smm_ge_unified.jl`: `rho = 0.04`, `chi = 2`, `gamma = 0.5`, `phi = 0.5` and `alpha = 0.7344`.

`run_full_calibration.jl` solves the parameters as a nested fixed point:

1. **Innovation block** (`calibrate_service_kappas.jl`). For a given demand triple, it searches over `(kappa_rs, kappa_ns)` with a damped Newton step until the model and trend slopes of RI- and NRI-Services labor productivity are equal. At each candidate, `calibrate_zeta_blocks.jl` re-fits research productivity. It alternates two sub-blocks:
   - the Goods sub-block (`zeta_m(1950), gzeta_m, kappa_m`), fitted to the Goods labor-productivity trend and patent-stock paths;
   - the services sub-block (the knots), fitted to the R&D-worker share paths.
2. **Demand block** (`calibrate_demand.jl`). Given the innovation parameters, it fits `(eta_m, eta_rs, sigma)` by bounded Levenberg–Marquardt to the Goods and RI-Services value-added share paths, and returns the deflator `ln D` for the next sweep.
3. **Repeat** steps 1–2 until the largest change across all parameters and `ln D` is below `1e-5`.
4. **Export.** `dump_fit_paths.jl` writes:
   - `final_calibration.txt`, the parameters with data and model moments side by side;
   - `fit_paths.csv`, the annual model and data paths;
   - `calibration_results.toml` (written by the driver) holds the knots, `ln D` and the convergence status.

Each block runs in a fresh Julia process, because `smm_ge_unified.jl` fixes the model specification through environment variables when it is loaded. The drivers set these variables themselves, so no setup is needed.

`baseline.jl` reads `output/calibration/final_calibration.txt` and `calibration_results.toml` and solves the model once at those parameters. The calibration figures and every counterfactual exercise start from this solved baseline.

## Outputs and the paper

| Paper item | File | Script |
|------------|------|--------|
| Figures 1, 2, 3, 4, 5 | `output/figures/Figure_1.pdf` … `Figure_5.pdf` | `motivation_plots.jl` |
| Figure B.1 | `output/figures/Figure_B1.pdf` | `motivation_plots.jl` |
| Figure 6 (externally calibrated inputs) | `output/figures/Figure_6.pdf` | `calibration_plots.jl` |
| Figure 7 (internal calibration fit) | `output/figures/Figure_7.pdf` | `calibration_plots.jl` |
| Table 2 (internally calibrated parameters) | `output/calibration/Table_2.tex` | `calibration_table.jl` |
| Table 3, Table H.1 | `output/counterfactuals/Table_3.tex`, `Table_H1.tex` | `Exercise1.jl` |
| Table 4 | `output/counterfactuals/Table_4.tex` | `Exercise2.jl` |
| Tables 5, 6, H.3, H.4 | `output/counterfactuals/Table_5.tex`, `Table_6.tex`, `Table_H3.tex`, `Table_H4.tex` | `Exercise3.jl` |
| Table 7 | `output/counterfactuals/Table_7.tex` | `Exercise4.jl` |
| Table 8 | `output/counterfactuals/Table_8.tex` | `Exercise5.jl` |
| Table H.2 | `output/robustness/Table_H2.tex` | `robustness.jl` (case `inputs_1980`) |
| Table J.1 | `output/robustness/Table_J1.tex` | `robustness.jl` |
| Appendix F, footnotes 41 and 42 (markup checks: real estate excluded; BEA data in every year) | `output/robustness/robustness_summary.csv`, rows `markups_noRE` and `markups_BEA` | `robustness.jl` (cases `markups_noRE`, `markups_BEA`) |

Table 1 lists externally set constants (see above). Other files in `output/counterfactuals/` (`Exercise4_shapley.tex`, `Exercise4b_*`, `Exercise3_subperiods.tex`, `Exercise2_stock_appendix.tex`, `Exercise1_channels.pdf`, `Figure_H1.pdf` and the CSVs) are supporting results that do not appear in the paper. `output/robustness/robustness_summary.csv` also holds the mean-markup runs with `phi = 0.25` and `0.75` quoted in Appendix J.

## Instructions

**Requirements:** Julia 1.12 (the versions of all packages are pinned in `Manifest.toml`). Run `../Empirical/master.do` first. Make sure `../Cross_Citations/output/` holds the three CSVs listed above (they are tracked, or run `../Cross_Citations/run_all.py`).

1. Open a terminal in this folder and install the packages once:
   ```
   julia --project=. -e 'using Pkg; Pkg.instantiate()'
   ```
2. Run everything:
   ```
   julia --project=. run_all.jl
   ```
   Step 5 builds Tables J.1 and H.2 and the markup checks of Appendix F from twelve alternative calibrations. It reuses the runs stored in `output/robustness/`; to recompute them, which takes about 4 hours, run `julia --project=. code/robustness.jl --redo`. To skip parts:
   ```
   julia --project=. run_all.jl --skip-calibration    # keep output/calibration as stored
   julia --project=. run_all.jl --skip-robustness     # leave out Tables J.1 and H.2
   ```

Each script can also be run on its own, e.g. `julia --project=. code/counterfactuals/Exercise3.jl`. Usage and arguments are in each file's header. Logs are written to `output/logs/`.

## Contact

Selim Elbadri: selim.elbadri@nottingham.ac.uk  
Juan I. Vizcaino: juan.vizcaino@nottingham.ac.uk
