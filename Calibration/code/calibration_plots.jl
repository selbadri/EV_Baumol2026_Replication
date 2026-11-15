#==========================================================================================
 An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
 ------------------------------------------------------------------------------------------
 File:    calibration_plots.jl
 Purpose: Draw the calibration figures: the externally calibrated inputs (Figure 6) and
          the fit of the internal calibration (Figure 7). Solves the baseline once.
 ------------------------------------------------------------------------------------------
 Usage:   julia --project=. code/calibration_plots.jl
 ------------------------------------------------------------------------------------------
 Inputs:  ../Empirical/output/dataset/data_moments.xlsx
          output/calibration/final_calibration.txt, calibration_results.toml (via baseline.jl)
 Outputs: output/figures/Figure_6.pdf, .png   Externally calibrated inputs
          output/figures/Figure_7.pdf, .png   Internal calibration model fit
==========================================================================================#

import Pkg
include(joinpath(@__DIR__, "paths.jl"))
Pkg.activate(CAL_PATHS.root)

using DataFrames
using Plots
using SmoothingSplines
using XLSX

const YEARS = collect(1950:2023)
const SPLINE_LAMBDA = 1.0e4
const OUTPUT_DIR = CAL_PATHS.out_figures
const PLOT_FONT = "Computer Modern"
const TICK_FONT_SIZE = 16
const GUIDE_FONT_SIZE = 20
const TITLE_FONT_SIZE = 19
const LEGEND_FONT_SIZE = 16
const MARKER_SCALE = 10.0

const DATA = DataFrame(XLSX.readtable(
	CAL_PATHS.data_moments,
	"DataMoments",
))
const DATA_YEARS = Int.(DATA[!, "Year"])

"Convert a spreadsheet column to Float64, representing missing values by NaN."
column(name::String) = Float64.(coalesce.(DATA[!, name], NaN))

"Return the observed years and finite observations for a raw data series."
function observations(values::AbstractVector{<:Real})
	keep = isfinite.(values) .& (DATA_YEARS .>= first(YEARS)) .&
		   (DATA_YEARS .<= last(YEARS))
	return DATA_YEARS[keep], Float64.(values[keep])
end

"Apply the same smoothing spline used when constructing the SMM inputs."
function smooth(values::AbstractVector{<:Real})
	# Fit on every observation, pre-1950 years included, as _spline_fill does for the
	# model inputs; restricting the fit to 1950-2023 bends the curve near 1950.
	keep = isfinite.(values)
	years, observed = DATA_YEARS[keep], Float64.(values[keep])
	length(observed) >= 2 || error("A smoothing spline requires at least two observations")
	spline = fit(SmoothingSpline, Float64.(years), observed, SPLINE_LAMBDA)
	return predict(spline, Float64.(YEARS))
end

"Construct a unit-normalized constant-growth path from its two endpoint levels."
function endpoint_growth_path(smoothed_levels::AbstractVector{<:Real})
	growth = log(smoothed_levels[end] / smoothed_levels[1]) / (length(YEARS) - 1)
	return exp.(growth .* (0:length(YEARS)-1))
end

const SECTORS = (
	(key=:m,  label="Goods",        color="#0072E8", marker=:diamond),
	(key=:rs, label="RI-Services",  color="#EF233C", marker=:square),
	(key=:ns, label="NRI-Services", color="#00A676", marker=:utriangle),
)

# Raw source series.
const PROFIT_SHARE = Dict(
	:m  => column("profit_share_GOODS"),
	:rs => column("profit_share_RI"),
	:ns => column("profit_share_NRI"),
)
const LABOR_VA_SHARE = Dict(
	:m  => column("LABVASHARE_GOODS"),
	:rs => column("LABVASHARE_RI"),
	:ns => column("LABVASHARE_NRI"),
)
const CAPITAL_PRICE = Dict(
	:m  => column("CAPPRICE_GOODS"),
	:rs => column("CAPPRICE_RI"),
	:ns => column("CAPPRICE_NRI"),
)
const LABOR_QUALITY = Dict(
	:m  => column("LABQLTY_GOODS"),
	:rs => column("LABQLTY_RI"),
	:ns => column("LABQLTY_NRI"),
)
const POPULATION = column("pop")
const SCIENTIST_SHARE = column("SCI_SH")
const VA_SHARE_SMOOTH = let
	goods = clamp.(smooth(column("VASHARE_GOODS")), 0.0, 1.0)
	ri_services = clamp.(smooth(column("VASHARE_RI")), 0.0, 1.0)
	nri_services = clamp.(1.0 .- goods .- ri_services, 0.0, 1.0)
	Dict(:m => goods, :rs => ri_services, :ns => nri_services)
end

# Load the calibrated baseline in an isolated module so its constants do not
# conflict with the plotting constants above.
module CalibrationPlotBaseline
	include(joinpath(@__DIR__, "baseline.jl"))
end
const MODEL_EQ = CalibrationPlotBaseline.eq_base
# The model's alpha and labor shares (net of land and inventory income under LS_NET), so the
# raw final-markup points use the same map as the model input.
const ALPHA = CalibrationPlotBaseline.alpha_fixed
const LABOR_VA_SHARE_MODEL = Dict(k => CalibrationPlotBaseline._ls_net(LABOR_VA_SHARE[k], CalibrationPlotBaseline._LS_NET[i])
                                  for (i, k) in enumerate((:m, :rs, :ns)))

"Normalize a positive series to one in 1950."
normalize_1950(values::AbstractVector{<:Real}) = Float64.(values) ./ values[1]

# Smoothed data and model paths used in the four calibration-fit panels.
const SCIENTIST_SECTOR_SHARE_SMOOTH = let
	goods = clamp.(smooth(column("RDSH_GOODS")), 1.0e-4, 0.98)
	ri_services = clamp.(smooth(column("RDSH_RI")), 1.0e-4, 0.98)
	scale = min.(1.0, 0.98 ./ max.(goods .+ ri_services, 1.0e-6))
	goods = clamp.(goods .* scale, 1.0e-4, 0.98)
	ri_services = clamp.(ri_services .* scale, 1.0e-4, 0.98)
	nri_services = clamp.(1.0 .- goods .- ri_services, 1.0e-4, 0.98)
	Dict(:m => goods, :rs => ri_services, :ns => nri_services)
end
const RELATIVE_PRICE_SMOOTH = Dict(
	:m  => normalize_1950(smooth(column("VAPRICE_GOODS"))),
	:rs => normalize_1950(smooth(column("VAPRICE_RI"))),
	:ns => normalize_1950(smooth(column("VAPRICE_NRI"))),
)
"Raw observations on YEARS (NaN where missing), indexed to the 1950 observation."
function raw_index_1950(values::AbstractVector{<:Real})
	pos = Dict(y => i for (i, y) in enumerate(DATA_YEARS))
	v = [haskey(pos, y) ? Float64(values[pos[y]]) : NaN for y in YEARS]
	isfinite(v[1]) || error("No 1950 observation to index to")
	return v ./ v[1]
end
# Panel (d) shows raw labor productivity as markers and, as thin solid lines, the
# smoothing-spline trend of ln LP that the calibration targets (1950 = 1, like the model).
const LABOR_PRODUCTIVITY_DATA = Dict(
	:m  => raw_index_1950(column("LP_GOODS")),
	:rs => raw_index_1950(column("LP_RI")),
	:ns => raw_index_1950(column("LP_NRI")),
)
_lp_trend_plot(c) = (s = smooth(log.(column(c))); exp.(s .- s[1]))
const LABOR_PRODUCTIVITY_TREND = Dict(:m  => _lp_trend_plot("LP_GOODS"),
                                      :rs => _lp_trend_plot("LP_RI"),
                                      :ns => _lp_trend_plot("LP_NRI"))
# Value-added shares, omega_j s_j / sum_i omega_i s_i; the solver's s_j are shares of
# gross revenue.
const MODEL_VA_SHARE = let va = CalibrationPlotBaseline._va_shares(MODEL_EQ; T=length(YEARS))
	Dict(:m => va.m, :rs => va.rs, :ns => va.ns)
end
const MODEL_SCIENTIST_SHARE = let
	total = Float64.(MODEL_EQ["Lst"][1:length(YEARS)])
	Dict(
		:m  => Float64.(MODEL_EQ["Lsmt"][1:length(YEARS)]) ./ total,
		:rs => Float64.(MODEL_EQ["Lsrst"][1:length(YEARS)]) ./ total,
		:ns => Float64.(MODEL_EQ["Lsnst"][1:length(YEARS)]) ./ total,
	)
end
# Model prices relative to a Fisher chain of the three sector prices with value-added
# weights, the counterpart of the private VA deflator that deflates the data.
const MODEL_PRICE_DEFLATOR = let n = length(YEARS)
	P = (MODEL_EQ["pmt"], MODEL_EQ["prst"], MODEL_EQ["pnst"])
	W = (MODEL_VA_SHARE[:m], MODEL_VA_SHARE[:rs], MODEL_VA_SHARE[:ns])
	d = ones(n)
	for t in 2:n
		rel = [P[j][t] / P[j][t-1] for j in 1:3]
		lasp = sum(W[j][t-1] * rel[j] for j in 1:3)
		paas = 1 / sum(W[j][t] / rel[j] for j in 1:3)
		d[t] = d[t-1] * sqrt(lasp * paas)
	end
	d
end
const MODEL_RELATIVE_PRICE = Dict(
	:m  => normalize_1950(MODEL_EQ["pmt"][1:length(YEARS)] ./ MODEL_PRICE_DEFLATOR),
	:rs => normalize_1950(MODEL_EQ["prst"][1:length(YEARS)] ./ MODEL_PRICE_DEFLATOR),
	:ns => normalize_1950(MODEL_EQ["pnst"][1:length(YEARS)] ./ MODEL_PRICE_DEFLATOR),
)
const MODEL_LABOR_PRODUCTIVITY = Dict(
	:m  => normalize_1950(MODEL_EQ["Ymt"][1:length(YEARS)] ./ MODEL_EQ["Lpmt"][1:length(YEARS)]),
	:rs => normalize_1950(MODEL_EQ["Yrst"][1:length(YEARS)] ./ MODEL_EQ["Lprst"][1:length(YEARS)]),
	:ns => normalize_1950(MODEL_EQ["Ynst"][1:length(YEARS)] ./ MODEL_EQ["Lpnst"][1:length(YEARS)]),
)

# The model inputs μˣ and μʸ are inverse markups. Plot their reciprocals so the
# first two panels report conventional gross markups greater than one.
const MUX_RAW = Dict(j.key => 1.0 .- PROFIT_SHARE[j.key] for j in SECTORS)
# Input panels plot the model's own input paths (INPUT_TRENDS in smm_ge_unified.jl), so the
# dashed lines are exactly what the solver uses. Target panels keep the spline.
const _EXOG = CalibrationPlotBaseline.baseline_exog
_input(name, key) = Float64.(getfield(_EXOG, Symbol(name, key))[1:length(YEARS)])
const MUX_SMOOTH = Dict(j.key => _input("mux", j.key) for j in SECTORS)
# The labor share of VALUE ADDED is alpha mu^y/(1-(1-alpha) mu^y mu^x), so mu^y needs
# mu^x as well as the labor share. The smoothed series is the model input.
muy_from_labva(ls, mux) = ls ./ (ALPHA .+ (1.0 - ALPHA) .* ls .* mux)
# The model's mu^x path (the quadratic trend) on the data years, NaN outside 1950-2023. The raw
# final-markup and machine-productivity points use it, as the model inputs do, so the measured
# markup swings of panel (b) do not enter panels (a) and (c).
_on_data_years(v) = [y in YEARS ? v[y - YEARS[1] + 1] : NaN for y in DATA_YEARS]
const MUX_TREND_DATA = Dict(j.key => _on_data_years(MUX_SMOOTH[j.key]) for j in SECTORS)
const MUY_RAW = Dict(j.key => muy_from_labva(LABOR_VA_SHARE_MODEL[j.key], MUX_TREND_DATA[j.key]) for j in SECTORS)
const MUY_SMOOTH = Dict(j.key => _input("muy", j.key) for j in SECTORS)
for j in SECTORS
	# The model's mu^y is built from the trend labor share and mu^x through the same map.
	all(isfinite, MUY_SMOOTH[j.key]) || error("mu^y input not finite for $(j.label)")
end
const MARKUP_X_RAW = Dict(j.key => 1.0 ./ MUX_RAW[j.key] for j in SECTORS)
const MARKUP_X_SMOOTH = Dict(j.key => 1.0 ./ MUX_SMOOTH[j.key] for j in SECTORS)
const MARKUP_Y_RAW = Dict(j.key => 1.0 ./ MUY_RAW[j.key] for j in SECTORS)
const MARKUP_Y_SMOOTH = Dict(j.key => 1.0 ./ MUY_SMOOTH[j.key] for j in SECTORS)

# Machine productivity: z_j = 1/(p^K_j * mu^x_j), at the model's mu^x path. The model uses
# constant growth between the 1950 and 2023 observations.
const Z_RAW = Dict{Symbol,Vector{Float64}}()
const Z_SMOOTH = Dict{Symbol,Vector{Float64}}()
for j in SECTORS
	# Raw z is normalized by its first observation below, so the price level cancels.
	pcapital_raw_index = CAPITAL_PRICE[j.key]
	Z_RAW[j.key] = 1.0 ./ (pcapital_raw_index .* MUX_TREND_DATA[j.key])
	first_observation = findfirst(isfinite, Z_RAW[j.key])
	first_observation === nothing && error("No machine-productivity data for $(j.label)")
	Z_RAW[j.key] ./= Z_RAW[j.key][first_observation]

	Z_SMOOTH[j.key] = _input("z", j.key) ./ _input("z", j.key)[1]
end

# Labor quality: constant growth between the 1950 and 2023 observations, unit-normalized.
const HBAR_RAW = Dict{Symbol,Vector{Float64}}()
const HBAR_SMOOTH = Dict{Symbol,Vector{Float64}}()
for j in SECTORS
	first_quality = LABOR_QUALITY[j.key][findfirst(isfinite, LABOR_QUALITY[j.key])]
	HBAR_RAW[j.key] = LABOR_QUALITY[j.key] ./ first_quality
	HBAR_SMOOTH[j.key] = _input("hbar", j.key) ./ _input("hbar", j.key)[1]
end

function parameter_panel(raw::Dict, fitted::Dict, title::String; legend=:best)
	panel = plot(
		xlabel="Year",
		title=title,
		xlims=(first(YEARS), last(YEARS)),
		xticks=1950:10:2020,
		legend=legend,
		grid=false,
		framestyle=:box,
		fontfamily=PLOT_FONT,
		tickfontsize=TICK_FONT_SIZE,
		guidefontsize=GUIDE_FONT_SIZE,
		titlefontsize=TITLE_FONT_SIZE,
		legendfontsize=LEGEND_FONT_SIZE,
		bottom_margin=11 * Plots.mm,
	)

	for j in SECTORS
		raw_years, raw_values = observations(raw[j.key])
		year_indices = raw_years .- first(YEARS) .+ 1
		# Plots interprets markersize as a radius-like quantity, so using the
		# square root makes marker area proportional to the sector's VA share.
		# A restrained scale and thin outline preserve this encoding without
		# allowing large-sector observations to obscure neighboring years.
		marker_sizes = MARKER_SCALE .* sqrt.(VA_SHARE_SMOOTH[j.key][year_indices])
		scatter!(
			panel,
			raw_years,
			raw_values;
			label=j.label,
			color=j.color,
			markershape=j.marker,
			markersize=marker_sizes,
			markerstrokewidth=0.4,
			markercolor=:white,
			markerstrokecolor=j.color,
		)
		plot!(
			panel,
			YEARS,
			fitted[j.key];
			label=false,
			color=j.color,
			linewidth=2.0,
			linestyle=:dash,
		)
	end
	return panel
end

"Plot one aggregate data series and the corresponding smoothed model input."
function aggregate_panel(raw_years, raw_values, fitted_years, fitted_values,
					 title::String; legend=:best, markershape=:diamond,
					 scale=1.0)
	panel = plot(
		xlabel="Year",
		title=title,
		xlims=(first(YEARS), last(YEARS)),
		xticks=1950:10:2020,
		legend=legend,
		grid=false,
		framestyle=:box,
		fontfamily=PLOT_FONT,
		tickfontsize=TICK_FONT_SIZE,
		guidefontsize=GUIDE_FONT_SIZE,
		titlefontsize=TITLE_FONT_SIZE,
		legendfontsize=LEGEND_FONT_SIZE,
		bottom_margin=11 * Plots.mm,
	)
	scatter!(
		panel,
		raw_years,
		scale .* raw_values;
		label="Data",
		color=:black,
		markershape=markershape,
		markersize=7.0,
		markerstrokewidth=0.4,
		markercolor=:white,
		markerstrokecolor=:black,
	)
	plot!(
		panel,
		fitted_years,
		scale .* fitted_values;
		label=false,
		color=:black,
		linewidth=2.0,
		linestyle=:dash,
	)
	return panel
end

"Plot sectoral data as markers, calibrated model paths as dashed lines and, when given, a data trend as thin solid lines."
function calibration_fit_panel(data::Dict, model::Dict, title::String;
							   legend=false, scale=1.0, trend=nothing)
	panel = plot(
		xlabel="Year",
		title=title,
		xlims=(first(YEARS), last(YEARS)),
		xticks=1950:10:2020,
		legend=legend,
		grid=false,
		framestyle=:box,
		fontfamily=PLOT_FONT,
		tickfontsize=TICK_FONT_SIZE,
		guidefontsize=GUIDE_FONT_SIZE,
		titlefontsize=TITLE_FONT_SIZE,
		legendfontsize=LEGEND_FONT_SIZE,
		bottom_margin=11 * Plots.mm,
	)
	for j in SECTORS
		scatter!(panel, YEARS, scale .* data[j.key];
			label=j.label,
			color=j.color,
			markershape=j.marker,
			markersize=MARKER_SCALE .* sqrt.(VA_SHARE_SMOOTH[j.key]),
			markerstrokewidth=0.4,
			markercolor=:white,
			markerstrokecolor=j.color)
		plot!(panel, YEARS, scale .* model[j.key];
			label=false, color=j.color,
			linewidth=2.2, linestyle=:dash)
		trend === nothing || plot!(panel, YEARS, scale .* trend[j.key];
			label=false, color=j.color, linewidth=1.2, linestyle=:solid)
	end
	return panel
end

# Population growth is the annual log change used by the model, shown in percent.
# Raw growth is reported only where adjacent calendar-year population observations
# are both available.
const POPULATION_SMOOTH_LEVEL = Float64.(_EXOG.Lt[1:length(YEARS)])
const POPULATION_GROWTH_SMOOTH = 100.0 .* diff(log.(POPULATION_SMOOTH_LEVEL))
const POPULATION_GROWTH_YEARS, POPULATION_GROWTH_RAW = let
	years = Int[]
	values = Float64[]
	for i in 2:length(DATA_YEARS)
		if DATA_YEARS[i] == DATA_YEARS[i-1] + 1 &&
		   first(YEARS) < DATA_YEARS[i] <= last(YEARS) &&
		   isfinite(POPULATION[i-1]) && isfinite(POPULATION[i]) &&
		   POPULATION[i-1] > 0.0 && POPULATION[i] > 0.0
			push!(years, DATA_YEARS[i])
			push!(values, 100.0 * log(POPULATION[i] / POPULATION[i-1]))
		end
	end
	years, values
end

# This is the same clamped spline fed to the baseline as the exogenous aggregate
# research-employment share.
const SCIENTIST_SHARE_SMOOTH = Float64.(_EXOG.ls[1:length(YEARS)])
const SCIENTIST_SHARE_YEARS, SCIENTIST_SHARE_RAW = observations(SCIENTIST_SHARE)

p_mux = parameter_panel(
	MARKUP_X_RAW,
	MARKUP_X_SMOOTH,
	"(b) Machine Producer Markups",
	legend=false,
)
p_muy = parameter_panel(
	MARKUP_Y_RAW,
	MARKUP_Y_SMOOTH,
	"(a) Final Producer Markups",
	legend=false,
)
p_z = parameter_panel(
	Z_RAW,
	Z_SMOOTH,
	"(c) Machine Producer Productivity (1950 = 1)",
	legend=:topleft,
)
p_hbar = parameter_panel(
	HBAR_RAW,
	HBAR_SMOOTH,
	"(d) Labor-Enhancing Productivity (1950 = 1)",
	legend=false,
)
p_population = aggregate_panel(
	POPULATION_GROWTH_YEARS,
	POPULATION_GROWTH_RAW,
	YEARS[2:end],
	POPULATION_GROWTH_SMOOTH,
	"(e) Population Growth Rate (%)",
	legend=:topright,
	markershape=:hexagon,
)
p_scientists = aggregate_panel(
	SCIENTIST_SHARE_YEARS,
	SCIENTIST_SHARE_RAW,
	YEARS,
	SCIENTIST_SHARE_SMOOTH,
	"(f) Economy-Wide R&D Worker Share (%)",
	legend=false,
	markershape=:hexagon,
	scale=100.0,
)

figure = plot(
	p_muy,
	p_mux,
	p_z,
	p_hbar,
	p_population,
	p_scientists;
	layout=(3, 2),
	size=(1400, 1300),
	margin=5 * Plots.mm,
)

mkpath(OUTPUT_DIR)
png_path = joinpath(OUTPUT_DIR, "Figure_6.png")
pdf_path = joinpath(OUTPUT_DIR, "Figure_6.pdf")
savefig(figure, png_path)
savefig(figure, pdf_path)

println("Saved: $png_path")
println("Saved: $pdf_path")

p_va_fit = calibration_fit_panel(
	VA_SHARE_SMOOTH,
	MODEL_VA_SHARE,
	"(a) Value Added Shares (%)";
	legend=(0.68, 0.68),
	scale=100.0,
)
p_scientist_fit = calibration_fit_panel(
	SCIENTIST_SECTOR_SHARE_SMOOTH,
	MODEL_SCIENTIST_SHARE,
	"(b) R&D Worker Shares (%)";
	legend=false,
	scale=100.0,
)
p_price_fit = calibration_fit_panel(
	RELATIVE_PRICE_SMOOTH,
	MODEL_RELATIVE_PRICE,
	"(c) Relative Prices (1950 = 1)";
	legend=false,
)
p_lp_fit = calibration_fit_panel(
	LABOR_PRODUCTIVITY_DATA,
	MODEL_LABOR_PRODUCTIVITY,
	"(d) Labor Productivity (1950 = 1)";
	legend=false,
	trend=LABOR_PRODUCTIVITY_TREND,
)

fit_figure = plot(
	p_va_fit,
	p_scientist_fit,
	p_price_fit,
	p_lp_fit;
	layout=(2, 2),
	size=(1400, 950),
	margin=5 * Plots.mm,
)
fit_png_path = joinpath(OUTPUT_DIR, "Figure_7.png")
fit_pdf_path = joinpath(OUTPUT_DIR, "Figure_7.pdf")
savefig(fit_figure, fit_png_path)
savefig(fit_figure, fit_pdf_path)

println("Saved: $fit_png_path")
println("Saved: $fit_pdf_path")

display(figure)
display(fit_figure)
