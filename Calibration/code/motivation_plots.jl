#==========================================================================================
 An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
 ------------------------------------------------------------------------------------------
 File:    motivation_plots.jl
 Purpose: Draw the data figures of Section 2 and Appendix B. Uses the data only; no
          model solve.
 ------------------------------------------------------------------------------------------
 Usage:   julia --project=. code/motivation_plots.jl
 ------------------------------------------------------------------------------------------
 Inputs:  ../Empirical/output/dataset/data_moments.xlsx
          ../Empirical/input/manual/NSF/nsf_rd_benchmarks.csv
          ../Cross_Citations/output/citations/spillover_backward1.csv
 Outputs: output/figures/Figure_1.pdf    Structural change in innovation
          output/figures/Figure_2.pdf    Heterogeneity within services
          output/figures/Figure_3.pdf    Citation-based knowledge linkages
          output/figures/Figure_4.pdf    Relative prices and value-added shares
          output/figures/Figure_5.pdf    Labor productivity across sectors
          output/figures/Figure_B1.pdf   R&D employment: IPUMS measure and NSF surveys
==========================================================================================#

import Pkg
include(joinpath(@__DIR__, "paths.jl"))
Pkg.activate(CAL_PATHS.root)

using DataFrames
using DelimitedFiles
using Plots
using Printf
using SmoothingSplines
using XLSX

const DATA_PATH = CAL_PATHS.data_moments
const SPILLOVER_PATH = CAL_PATHS.spillovers
# The README in Empirical/input/manual/NSF documents every cell.
const NSF_PATH = CAL_PATHS.nsf_benchmarks
const OUTPUT_DIR = CAL_PATHS.out_figures
const SPLINE_LAMBDA = 1.0e4
const ANNUAL_YEARS = collect(1947:2023)
const SHARE_YEARS = collect(1950:2023)
const PLOT_FONT = "Computer Modern"
const TICK_FONT_SIZE = 16
const GUIDE_FONT_SIZE = 20
const TITLE_FONT_SIZE = 19
const LEGEND_FONT_SIZE = 16
const MARKER_SCALE = 10.0

const DATA = DataFrame(XLSX.readtable(DATA_PATH, "DataMoments"))
const DATA_YEARS = Int.(DATA[!, "Year"])

const SECTORS = (
	(key=:m,  label="Goods",     color="#0072E8", marker=:diamond),
	(key=:rs, label="RI-Serv.",  color="#EF233C", marker=:square),
	(key=:ns, label="NRI-Serv.", color="#00A676", marker=:utriangle),
)
const SERVICES_COLOR = "#C000FF"

"Read one spreadsheet column, replacing missing entries with NaN."
column(name::String) = Float64.(coalesce.(DATA[!, name], NaN))

"Return finite observations, optionally restricted to an inclusive year range."
function observations(values::AbstractVector{<:Real};
					  first_year=minimum(DATA_YEARS), last_year=maximum(DATA_YEARS))
	keep = isfinite.(values) .& (DATA_YEARS .>= first_year) .&
		   (DATA_YEARS .<= last_year)
	return DATA_YEARS[keep], Float64.(values[keep])
end

"Fit the same cubic smoothing spline and penalty used in calibration_plots.jl."
function smooth(values::AbstractVector{<:Real}, prediction_years::AbstractVector{<:Real})
	years, observed = observations(values)
	length(observed) >= 2 || error("A smoothing spline requires at least two observations")
	spline = fit(SmoothingSpline, Float64.(years), observed, SPLINE_LAMBDA)
	return predict(spline, Float64.(prediction_years))
end

const VA_SHARE_RAW = Dict(
	:m  => column("VASHARE_GOODS"),
	:rs => column("VASHARE_RI"),
	:ns => column("VASHARE_NRI"),
)
const VA_SHARE_SMOOTH = Dict(
	j.key => smooth(VA_SHARE_RAW[j.key], ANNUAL_YEARS) for j in SECTORS
)

const CAPITAL_SHARE_RAW = Dict(
	:m  => 100.0 .* (1.0 .- column("LABVASHARE_GOODS")),
	:rs => 100.0 .* (1.0 .- column("LABVASHARE_RI")),
	:ns => 100.0 .* (1.0 .- column("LABVASHARE_NRI")),
)
const VALUE_ADDED_SHARE_RAW = Dict(
	j.key => 100.0 .* VA_SHARE_RAW[j.key] for j in SECTORS
)
const LABOR_PRODUCTIVITY_RAW = Dict(
	:m  => column("LP_GOODS"),
	:rs => column("LP_RI"),
	:ns => column("LP_NRI"),
)
const RELATIVE_PRICE_RAW = Dict(
	:m  => column("VAPRICE_GOODS"),
	:rs => column("VAPRICE_RI"),
	:ns => column("VAPRICE_NRI"),
)
const SCIENTIST_SHARE_RAW = Dict(
	:m  => 100.0 .* column("RDSH_GOODS"),
	:rs => 100.0 .* column("RDSH_RI"),
	:ns => 100.0 .* column("RDSH_NRI"),
)
const PATENT_SHARE_RAW = Dict(
	:m  => 100.0 .* column("PATSH_GOODS"),
	:rs => 100.0 .* column("PATSH_RI"),
	:ns => 100.0 .* column("PATSH_NRI"),
)
const AGGREGATE_LABOR_PRODUCTIVITY_RAW = column("LP_AGGREGATE")
const RD_EMPLOYMENT_SHARE_RAW = 100.0 .* column("SCI_SH")

"Normalize a series by its observed value in a specified base year."
function normalize_at_year(values::AbstractVector{<:Real}, base_year::Int)
	base_index = findfirst(==(base_year), DATA_YEARS)
	base_index === nothing && error("Base year $base_year is absent from the data")
	base_value = values[base_index]
	isfinite(base_value) && base_value > 0.0 ||
		error("Series has no positive finite value in $base_year")
	return Float64.(values) ./ base_value
end

"Continuous average annual log growth between two observed endpoints, in percent."
function average_growth_percent(values::AbstractVector{<:Real}, first_year::Int, last_year::Int)
	first_index = findfirst(==(first_year), DATA_YEARS)
	last_index = findfirst(==(last_year), DATA_YEARS)
	(first_index !== nothing && last_index !== nothing) ||
		error("Growth endpoints must be present in the data")
	return 100.0 * log(values[last_index] / values[first_index]) /
		   (last_year - first_year)
end

"Legend text reporting the indexed endpoint and calibration-target average growth."
growth_legend(label::String, normalized_values, raw_values) = @sprintf(
	"%s: 2023 = %.2f; avg. growth = %.1f%%",
	label,
	normalized_values[findfirst(==(2023), DATA_YEARS)],
	average_growth_percent(raw_values, 1950, 2023),
)

# Knowledge-spillover weights. Rows index the receiving sector and columns the
# source sector, following the convention in smm_ge_unified.jl.
const SPILLOVER_TABLE, SPILLOVER_HEADER = readdlm(SPILLOVER_PATH, ',', header=true)
const SPILLOVER_COLUMNS = String.(vec(SPILLOVER_HEADER))
# Citing years 1978-2020, as in the calibration: 1977 rests on 28 citations and citing years
# after 2020 are incomplete. The splines are fitted and drawn over the same years.
const SPILLOVER_FIRST, SPILLOVER_LAST = 1978, 2020
const _SPILL_ROWS = let y = round.(Int, Float64.(SPILLOVER_TABLE[:, findfirst(==("start"), SPILLOVER_COLUMNS)]))
	findall(t -> SPILLOVER_FIRST <= t <= SPILLOVER_LAST, y)
end
spillover_column(name::String) =
	Float64.(SPILLOVER_TABLE[_SPILL_ROWS, findfirst(==(name), SPILLOVER_COLUMNS)])
const SPILLOVER_YEARS = Int.(round.(spillover_column("start")))
const SPILLOVER_RAW = Dict(
	:m => Dict(
		:m => spillover_column("delta_mm"),
		:rs => spillover_column("delta_mr"),
		:ns => spillover_column("delta_mn"),
	),
	:rs => Dict(
		:m => spillover_column("delta_rm"),
		:rs => spillover_column("delta_rr"),
		:ns => spillover_column("delta_rn"),
	),
	:ns => Dict(
		:m => spillover_column("delta_nm"),
		:rs => spillover_column("delta_nr"),
		:ns => spillover_column("delta_nn"),
	),
)

function smooth_spillover(values::AbstractVector{<:Real})
	spline = fit(
		SmoothingSpline,
		Float64.(SPILLOVER_YEARS),
		Float64.(values),
		SPLINE_LAMBDA,
	)
	return predict(spline, Float64.(SPILLOVER_YEARS))
end

"Smooth and renormalize a receiving-sector row for display; the calibration holds each probability at its 1978-2020 mean."
function normalized_spillover_row(receiving::Symbol)
	fitted = Dict(j.key => smooth_spillover(SPILLOVER_RAW[receiving][j.key]) for j in SECTORS)
	row_total = fitted[:m] .+ fitted[:rs] .+ fitted[:ns]
	return Dict(j.key => fitted[j.key] ./ row_total for j in SECTORS)
end

"Marker sizes with area proportional to the sector's smoothed value-added share."
function sector_marker_sizes(key::Symbol, years::AbstractVector{<:Integer})
	indices = years .- first(ANNUAL_YEARS) .+ 1
	return MARKER_SCALE .* sqrt.(clamp.(VA_SHARE_SMOOTH[key][indices], 0.0, 1.0))
end

function base_panel(ylabel::String, years::AbstractVector{<:Integer}; legend=false, title="")
	return plot(
		xlabel="Year",
		ylabel=ylabel,
		title=title,
		xlims=(first(years), last(years)),
		xticks=first(years):10:(last(years) - 3),
		legend=legend,
		grid=false,
		framestyle=:box,
		fontfamily=PLOT_FONT,
		tickfontsize=TICK_FONT_SIZE,
		guidefontsize=GUIDE_FONT_SIZE,
		titlefontsize=TITLE_FONT_SIZE,
		legendfontsize=LEGEND_FONT_SIZE,
		left_margin=13 * Plots.mm,
		bottom_margin=12 * Plots.mm,
		top_margin=7 * Plots.mm,
	)
end

"Plot the three source-sector spillover weights for one receiving sector."
function spillover_panel(receiving::Symbol, title::String; legend=false)
	panel = base_panel(
		"Share of Citations (%)",
		SHARE_YEARS;
		legend=legend,
		title=title,
	)
	plot!(panel; left_margin=15 * Plots.mm,
		xlims=(SPILLOVER_FIRST - 2, SPILLOVER_LAST + 2), xticks=1980:10:2020)
	fitted = normalized_spillover_row(receiving)
	for source in SECTORS
		raw_values = SPILLOVER_RAW[receiving][source.key]
		scatter!(panel, SPILLOVER_YEARS, 100.0 .* raw_values;
			label=source.label,
			color=source.color,
			markershape=source.marker,
			markersize=sector_marker_sizes(source.key, SPILLOVER_YEARS),
			markerstrokewidth=0.4,
			markercolor=:white,
			markerstrokecolor=source.color,
		)
		plot!(panel, SPILLOVER_YEARS, 100.0 .* fitted[source.key];
			label=false,
			color=source.color,
			linewidth=2.0,
			linestyle=:dash,
		)
	end
	return panel
end

"Plot sector observations as distinct hollow markers and smoothing splines as dashed lines."
function sector_panel(raw::Dict, years::Vector{Int}, ylabel::String;
					  legend=false, title="", normalization_year=nothing)
	panel = base_panel(ylabel, years; legend=legend, title=title)
	for j in SECTORS
		plotted_raw = normalization_year === nothing ? raw[j.key] :
			normalize_at_year(raw[j.key], normalization_year)
		raw_years, raw_values = observations(
			plotted_raw; first_year=first(years), last_year=last(years)
		)
		fitted = smooth(raw[j.key], years)
		if normalization_year !== nothing
			base_index = findfirst(==(normalization_year), years)
			base_index === nothing && error("Normalization year $normalization_year is outside the plot range")
			fitted ./= fitted[base_index]
		end
		scatter!(panel, raw_years, raw_values;
			label=j.label,
			color=j.color,
			markershape=j.marker,
			markersize=sector_marker_sizes(j.key, raw_years),
			markerstrokewidth=0.4,
			markercolor=:white,
			markerstrokecolor=j.color,
		)
		plot!(panel, years, fitted;
			label=false,
			color=j.color,
			linewidth=2.0,
			linestyle=:dash,
		)
	end
	return panel
end

"Plot goods and total-services shares, with service shares summed across both types."
function aggregate_services_panel(raw::Dict, years::Vector{Int}, ylabel::String;
								  legend=false, title="")
	panel = base_panel(ylabel, years; legend=legend, title=title)
	series = (
		(label="Goods", color="#0072E8", marker=:diamond, raw=raw[:m],
		 smooth=smooth(raw[:m], years), va_key=:m),
		(label="Services", color=SERVICES_COLOR, marker=:star5, raw=raw[:rs] .+ raw[:ns],
		 smooth=smooth(raw[:rs], years) .+ smooth(raw[:ns], years), va_key=:services),
	)
	for item in series
		raw_years, raw_values = observations(
			item.raw; first_year=first(years), last_year=last(years)
		)
		marker_sizes = if item.va_key == :m
			sector_marker_sizes(:m, raw_years)
		else
			indices = raw_years .- first(ANNUAL_YEARS) .+ 1
			service_va = VA_SHARE_SMOOTH[:rs][indices] .+ VA_SHARE_SMOOTH[:ns][indices]
			MARKER_SCALE .* sqrt.(clamp.(service_va, 0.0, 1.0))
		end
		scatter!(panel, raw_years, raw_values;
			label=item.label,
			color=item.color,
			markershape=item.marker,
			markersize=marker_sizes,
			markerstrokewidth=0.4,
			markercolor=:white,
			markerstrokecolor=item.color,
		)
		plot!(panel, years, item.smooth;
			label=false,
			color=item.color,
			linewidth=2.0,
			linestyle=:dash,
		)
	end
	return panel
end

"Plot sectoral labor-productivity indices normalized to one in 1950."
function sectoral_lp_panel(; legend=:topleft)
	panel = base_panel(
		"Labor Productivity (1950 = 1)",
		SHARE_YEARS;
		legend=legend,
		title="(a) Sectoral Labor Productivity",
	)
	long_labels = Dict(:m => "Goods", :rs => "RI-Services", :ns => "NRI-Services")
	for j in SECTORS
		normalized_raw = normalize_at_year(LABOR_PRODUCTIVITY_RAW[j.key], 1950)
		raw_years, raw_values = observations(
			normalized_raw; first_year=first(SHARE_YEARS), last_year=last(SHARE_YEARS)
		)
		smoothed = smooth(LABOR_PRODUCTIVITY_RAW[j.key], SHARE_YEARS)
		smoothed ./= smoothed[1]
		scatter!(panel, raw_years, raw_values;
			label=growth_legend(long_labels[j.key], normalized_raw,
				LABOR_PRODUCTIVITY_RAW[j.key]),
			color=j.color,
			markershape=j.marker,
			markersize=sector_marker_sizes(j.key, raw_years),
			markerstrokewidth=0.4,
			markercolor=:white,
			markerstrokecolor=j.color,
		)
		plot!(panel, SHARE_YEARS, smoothed;
			label=false,
			color=j.color,
			linewidth=2.0,
			linestyle=:dash,
		)
	end
	return panel
end

"Plot the data's VA-Törnqvist aggregate labor-productivity index."
function aggregate_lp_panel(; legend=:topleft)
	panel = base_panel(
		"Labor Productivity (1950 = 1)",
		SHARE_YEARS;
		legend=legend,
		title="(b) Aggregate Labor Productivity",
	)
	normalized_raw = normalize_at_year(AGGREGATE_LABOR_PRODUCTIVITY_RAW, 1950)
	raw_years, raw_values = observations(
		normalized_raw; first_year=first(SHARE_YEARS), last_year=last(SHARE_YEARS)
	)
	smoothed = smooth(AGGREGATE_LABOR_PRODUCTIVITY_RAW, SHARE_YEARS)
	smoothed ./= smoothed[1]
	scatter!(panel, raw_years, raw_values;
		label=growth_legend("Aggregate", normalized_raw,
			AGGREGATE_LABOR_PRODUCTIVITY_RAW),
		color=:black,
		markershape=:hexagon,
		markersize=7.5,
		markerstrokewidth=0.4,
		markercolor=:white,
		markerstrokecolor=:black,
	)
	plot!(panel, SHARE_YEARS, smoothed;
		label=false,
		color=:black,
		linewidth=2.0,
		linestyle=:dash,
	)
	return panel
end

const NSF_TABLE, NSF_HEADER = readdlm(NSF_PATH, ',', header=true)
nsf_column(name::String) = Float64.(replace(
	NSF_TABLE[:, findfirst(==(name), String.(vec(NSF_HEADER)))], "" => NaN))
const NSF_YEARS = Int.(nsf_column("year"))
const NSF_SECTOR_SHARE = Dict(
	:m  => nsf_column("goods_share"),
	:rs => nsf_column("ri_share"),
	:ns => nsf_column("nri_share"),
)
const NSF_MANUFACTURING_SHARE = nsf_column("mfg_share")
const NSF_RD_EMPLOYMENT_SHARE = nsf_column("rd_pct_priv")

"Finite NSF observations of one column, optionally restricted to a set of years."
function nsf_observations(values; years=NSF_YEARS)
	keep = isfinite.(values) .& in.(NSF_YEARS, Ref(years))
	return NSF_YEARS[keep], values[keep]
end

"IPUMS sector shares as in Figure 2, with the NSF survey shares as filled markers."
function nsf_sector_panel()
	panel = sector_panel(SCIENTIST_SHARE_RAW, SHARE_YEARS, "Share of R&D Workers (%)";
		legend=(0.1, 0.47), title="(a) Sectoral Composition")
	# The NSF benchmarks fall on 2023, so the frame leaves room for the endpoint markers
	plot!(panel; ylims=(0, 100), xlims=(1948, 2025))
	mfg_years, mfg_values = nsf_observations(NSF_MANUFACTURING_SHARE)
	scatter!(panel, mfg_years, mfg_values;
		label="NSF, manufacturing",
		markershape=:circle,
		markersize=4.0,
		markercolor=SECTORS[1].color,
		markerstrokecolor=SECTORS[1].color,
		markerstrokewidth=0.4,
	)
	for j in SECTORS
		nsf_years, nsf_values = nsf_observations(NSF_SECTOR_SHARE[j.key])
		scatter!(panel, nsf_years, nsf_values;
			label=false,
			markershape=j.marker,
			markersize=8.0,
			markercolor=j.color,
			markerstrokecolor=j.color,
			markerstrokewidth=0.4,
		)
	end
	# Legend key only: filled markers denote the NSF three-sector benchmarks
	scatter!(panel, [NaN], [NaN];
		label="NSF, three sectors",
		markershape=:diamond,
		markersize=8.0,
		markercolor=:gray55,
		markerstrokecolor=:gray55,
	)
	return panel
end

"R&D workers per 100 private employees: IPUMS against the NSF surveys."
function nsf_total_panel()
	panel = base_panel("R&D Workers per 100\nPrivate Employees", SHARE_YEARS;
		legend=:topleft, title="(b) R&D Employment")
	plot!(panel; ylims=(0, 2), xlims=(1948, 2025))
	raw_years, raw_values = observations(RD_EMPLOYMENT_SHARE_RAW;
		first_year=first(SHARE_YEARS), last_year=last(SHARE_YEARS))
	scatter!(panel, raw_years, raw_values;
		label="IPUMS (this paper)",
		markershape=:hexagon,
		markersize=7.5,
		markerstrokewidth=0.4,
		markercolor=:white,
		markerstrokecolor=:black,
	)
	plot!(panel, SHARE_YEARS, smooth(RD_EMPLOYMENT_SHARE_RAW, SHARE_YEARS);
		label=false,
		color=:black,
		linewidth=2.0,
		linestyle=:dash,
	)
	# 2023 is BERD's count of all R&D employees, a broader concept than the FTE
	# scientists and engineers of the earlier surveys
	fte_years, fte_values = nsf_observations(NSF_RD_EMPLOYMENT_SHARE; years=1950:2022)
	scatter!(panel, fte_years, fte_values;
		label="NSF, FTE R&D scientists and engineers",
		markershape=:circle,
		markersize=4.0,
		markercolor=:black,
		markerstrokecolor=:black,
		markerstrokewidth=0.4,
	)
	emp_years, emp_values = nsf_observations(NSF_RD_EMPLOYMENT_SHARE; years=2023:2023)
	scatter!(panel, emp_years, emp_values;
		label="NSF, R&D employees",
		markershape=:utriangle,
		markersize=8.0,
		markercolor=:black,
		markerstrokecolor=:black,
		markerstrokewidth=0.4,
	)
	return panel
end

function save_two_panel_figure(left, right, filename::String)
	figure = plot(
		left,
		right;
		layout=(1, 2),
		size=(1900, 800),
		margin=5 * Plots.mm,
	)
	path = joinpath(OUTPUT_DIR, filename)
	savefig(figure, path)
	println("Saved: $path")
	return figure
end

function save_three_panel_figure(left, center, right, filename::String)
	figure = plot(
		left,
		center,
		right;
		layout=(1, 3),
		size=(2500, 800),
		margin=5 * Plots.mm,
	)
	path = joinpath(OUTPUT_DIR, filename)
	savefig(figure, path)
	println("Saved: $path")
	return figure
end

mkpath(OUTPUT_DIR)

relative_prices_and_va = save_two_panel_figure(
	sector_panel(RELATIVE_PRICE_RAW, SHARE_YEARS,
		"Price of Sectoral Goods\nRelative to Final Goods (1950 = 1)";
		legend=:topleft, title="(a) Relative Prices", normalization_year=1950),
	sector_panel(VALUE_ADDED_SHARE_RAW, SHARE_YEARS,
		"Share of Value Added (%)"; legend=false,
		title="(b) Value Added Shares"),
	"Figure_4.pdf",
)

labor_productivity = save_two_panel_figure(
	sectoral_lp_panel(),
	aggregate_lp_panel(),
	"Figure_5.pdf",
)

sectoral_innovation = save_two_panel_figure(
	sector_panel(SCIENTIST_SHARE_RAW, SHARE_YEARS,
		"Share of R&D Workers (%)"; legend=:left, title="(a) R&D Workers"),
	sector_panel(PATENT_SHARE_RAW, SHARE_YEARS,
		"Share of Patents (%)"; legend=false, title="(b) Patents"),
	"Figure_2.pdf",
)

goods_services_innovation = save_two_panel_figure(
	aggregate_services_panel(SCIENTIST_SHARE_RAW, SHARE_YEARS,
		"Share of R&D Workers (%)"; legend=:top, title="(a) R&D Workers"),
	aggregate_services_panel(PATENT_SHARE_RAW, SHARE_YEARS,
		"Share of Patents (%)"; legend=false, title="(b) Patents"),
	"Figure_1.pdf",
)

knowledge_spillovers = save_three_panel_figure(
	spillover_panel(:m, "(a) Citing Sector: Goods"; legend=:left),
	spillover_panel(:rs, "(b) Citing Sector: RI-Services"; legend=false),
	spillover_panel(:ns, "(c) Citing Sector: NRI-Services"; legend=false),
	"Figure_3.pdf",
)

rd_workers_nsf = save_two_panel_figure(
	nsf_sector_panel(),
	nsf_total_panel(),
	"Figure_B1.pdf",
)

display(relative_prices_and_va)
display(labor_productivity)
display(sectoral_innovation)
display(goods_services_innovation)
display(knowledge_spillovers)
