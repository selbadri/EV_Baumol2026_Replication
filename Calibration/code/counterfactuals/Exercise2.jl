#==========================================================================================
 An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
 ------------------------------------------------------------------------------------------
 File:    Exercise2.jl
 Purpose: Exercise 2. Decomposition of the growth in the flow of new ideas along the
          baseline path, 1950-2023 (Section 5.2).
 ------------------------------------------------------------------------------------------
 Usage:   julia --project=. code/counterfactuals/Exercise2.jl
 ------------------------------------------------------------------------------------------
 Inputs:  output/calibration/final_calibration.txt, calibration_results.toml (via ../baseline.jl)
 Outputs: output/counterfactuals/Table_4.tex
          output/counterfactuals/Exercise2_stock_appendix.tex   stock decomposition (not in the paper)
 ------------------------------------------------------------------------------------------
 Notes:
 Exact decomposition of growth in the flow of new varieties along
 the baseline equilibrium path over 1950--2023.

 Main-text identity:

   Delta ln Ndot_j
     = Delta ln zeta_j
     + phi Delta ln L_j^s
     + gamma delta_jj Delta ln N_j
     + gamma sum(i != j) delta_ji Delta ln N_i,

 exact in annual differences because the spillover elasticities delta_ji are
 constant. The stock-based decomposition of Delta ln N_j goes to an appendix
 table.
==========================================================================================#

include(joinpath(@__DIR__, "..", "baseline.jl"))

using Printf

const YEAR0 = 1950
const ND = baseline_params.T_data
const GAMMA = baseline_params.gamma
const PHI = baseline_params.phi
const OUTDIR = CAL_PATHS.out_counterfact

_slice(v) = collect(v[1:ND])
# Research productivity along the baseline: the calibrated paths (knots for the services when
# the baseline has them), so the identity holds exactly; exponential otherwise.
const _ZETA_PATHS = zeta_paths(baseline_params)
_zeta(z0, g) = [z0 * exp(g * (t - 1)) for t in 1:ND]
_zeta_j(j, z0, g) = _ZETA_PATHS === nothing ? _zeta(z0, g) : collect(_ZETA_PATHS[j][1:ND])

struct SectorInputs
    name  :: String
    own   :: Int
    drow  :: NTuple{3,Float64}
    N     :: Vector{Float64}
    Ls    :: Vector{Float64}
    H     :: Vector{Float64}
    Ndot  :: Vector{Float64}
    zeta  :: Vector{Float64}
end

const N_PATHS = (
    _slice(eq_base["Nmt"]),
    _slice(eq_base["Nrst"]),
    _slice(eq_base["Nnst"]),
)

const N_LAMBDA_PATHS = (
    _slice(eq_base["NΛmt"]),
    _slice(eq_base["NΛrst"]),
    _slice(eq_base["NΛnst"]),
)

# The solver takes the spillover elasticities as per-period arrays; they are constant parameters.
function _delta(v)
    all(==(v[1]), v) || error("Exercise 2 requires constant spillover elasticities")
    v[1]
end

const SECTORS = [
    SectorInputs("Goods", 1,
        (_delta(baseline_params.dMM), _delta(baseline_params.dMR), _delta(baseline_params.dMN)),
        N_PATHS[1], _slice(eq_base["Lsmt"]), _slice(eq_base["Hmt"]),
        _slice(eq_base["Nmdot"]), _zeta_j(1, baseline_params.zetam_ini, baseline_params.gzetam)),
    SectorInputs("RI-Services", 2,
        (_delta(baseline_params.dRM), _delta(baseline_params.dRR), _delta(baseline_params.dRN)),
        N_PATHS[2], _slice(eq_base["Lsrst"]), _slice(eq_base["Hrst"]),
        _slice(eq_base["Nrsdot"]), _zeta_j(2, baseline_params.zetars_ini, baseline_params.gzetars)),
    SectorInputs("NRI-Services", 3,
        (_delta(baseline_params.dNM), _delta(baseline_params.dNR), _delta(baseline_params.dNN)),
        N_PATHS[3], _slice(eq_base["Lsnst"]), _slice(eq_base["Hnst"]),
        _slice(eq_base["Nnsdot"]), _zeta_j(3, baseline_params.zetans_ini, baseline_params.gzetans)),
]

const Q_INITIAL = let q=baseline_params.lambda_N0
    (q.m,q.rs,q.ns)
end

struct FlowDecomp
    productivity     :: Vector{Float64}
    employment       :: Vector{Float64}
    own_ideas        :: Vector{Float64}
    other_ideas      :: Vector{Float64}
    total            :: Vector{Float64}
    direct           :: Vector{Float64}
    resid            :: Vector{Float64}
    level_resid      :: Vector{Float64}
end

function decompose_flow(s::SectorInputs)
    productivity=zeros(ND); employment=zeros(ND)
    own_ideas=zeros(ND); other_ideas=zeros(ND)
    total=zeros(ND); direct=zeros(ND); resid=zeros(ND)
    level_resid=zeros(ND)

    lnN=ntuple(k -> log.(N_PATHS[k]),3)
    lnNlambda=ntuple(k -> log.(N_LAMBDA_PATHS[k]),3)
    lnNdot=log.(s.Ndot)

    # The innovation technology in levels uses the share-seeded stocks N^Lambda = q(1950) N; the
    # constant ln q drops out of the annual differences.
    for t in 1:ND
        rhs=log(s.zeta[t])+PHI*log(s.Ls[t])
        rhs+=GAMMA*sum(s.drow[k]*lnNlambda[k][t] for k in 1:3)
        level_resid[t]=rhs-lnNdot[t]
    end

    for t in 2:ND
        dlnN=ntuple(k -> lnN[k][t]-lnN[k][t-1],3)
        productivity[t]=log(s.zeta[t])-log(s.zeta[t-1])
        employment[t]=PHI*(log(s.Ls[t])-log(s.Ls[t-1]))
        own_ideas[t]=GAMMA*s.drow[s.own]*dlnN[s.own]
        other_ideas[t]=GAMMA*sum(s.drow[k]*dlnN[k] for k in 1:3 if k != s.own)
        total[t]=productivity[t]+employment[t]+own_ideas[t]+other_ideas[t]
        direct[t]=lnNdot[t]-lnNdot[t-1]
        resid[t]=total[t]-direct[t]
    end

    FlowDecomp(productivity,employment,own_ideas,other_ideas,total,direct,resid,level_resid)
end

# Spillover multiplier 1/(1-gamma delta_jj) on own ideas.
multiplier(s::SectorInputs)=1/(1-GAMMA*s.drow[s.own])

struct StockDecomp
    productivity :: Vector{Float64}
    labor        :: Vector{Float64}
    external     :: Vector{Float64}
    transition   :: Vector{Float64}
    total        :: Vector{Float64}
    direct       :: Vector{Float64}
    resid        :: Vector{Float64}
end

function decompose_stock(s::SectorInputs)
    productivity=zeros(ND); labor=zeros(ND); external=zeros(ND)
    transition=zeros(ND); total=zeros(ND)
    direct=zeros(ND); resid=zeros(ND)
    lnN=ntuple(k -> log.(N_PATHS[k]),3)
    lng=log.(s.Ndot./s.N)
    denominator=1-GAMMA*s.drow[s.own]
    denominator>0 || error("Nonpositive dilution denominator for $(s.name)")

    for t in 2:ND
        productivity[t]=(log(s.zeta[t])-log(s.zeta[t-1]))/denominator
        labor[t]=PHI*(log(s.Ls[t])-log(s.Ls[t-1]))/denominator
        external[t]=GAMMA*sum(s.drow[k]*(lnN[k][t]-lnN[k][t-1])
                              for k in 1:3 if k != s.own)/denominator
        transition[t]=-(lng[t]-lng[t-1])/denominator
        total[t]=productivity[t]+labor[t]+external[t]+transition[t]
        direct[t]=lnN[s.own][t]-lnN[s.own][t-1]
        resid[t]=total[t]-direct[t]
    end
    StockDecomp(productivity,labor,external,transition,total,direct,resid)
end

const FLOW_DECOMPS=[decompose_flow(s) for s in SECTORS]
const STOCK_DECOMPS=[decompose_stock(s) for s in SECTORS]

const FLOW_CHANNELS=(
    ("I. Research productivity",:productivity),
    ("II. Research employment",:employment),
    ("III. Own accumulated ideas",:own_ideas),
    ("IV. Other sectors' accumulated ideas",:other_ideas),
)
const STOCK_CHANNELS=(
    ("I. Research productivity",:productivity),
    ("II. Research labor",:labor),
    ("III. External spillovers",:external),
    ("IV. Transition",:transition),
)

cumulative(d,field::Symbol)=sum(getfield(d,field)[2:ND])

function validate_inputs()
    for j in 1:3
        s=SECTORS[j]
        all(>(0),s.N) && all(>(0),N_LAMBDA_PATHS[j]) && all(>(0),s.Ls) &&
            all(>(0),s.H) && all(>(0),s.Ndot) && all(>(0),s.zeta) ||
            error("Exercise 2 requires positive paths for $(s.name)")
        isapprox(N_PATHS[j][1],1.0;atol=1e-12,rtol=1e-12) ||
            error("Exercise 2 requires normalized N_j(1950)=1")
        maximum(abs.(N_LAMBDA_PATHS[j].-Q_INITIAL[j].*N_PATHS[j]))<=1e-12 ||
            error("Exercise 2 requires N_j^Lambda=q_j(1950)N_j")
        maximum(abs.(s.H.-s.Ls))<=1e-12 ||
            error("Exercise 2 requires homogeneous researchers with H_j=L_j^s")
        abs(sum(s.drow)-1)<=1e-10 ||
            error("Exercise 2 spillover row does not sum to one for $(s.name)")
        1-GAMMA*s.drow[j]>0 ||
            error("Exercise 2 spillover multiplier is undefined for $(s.name)")
    end
end

function run_validation()
    validate_inputs()
    println("\n"*"="^90)
    println("  Exercise 2 — exact flow and stock decomposition validation")
    println("="^90)
    for j in 1:3
        s=SECTORS[j]; f=FLOW_DECOMPS[j]; a=STOCK_DECOMPS[j]
        levelgap=maximum(abs.(f.level_resid))
        annualgap=maximum(abs.(f.resid[2:end]))
        cumgap=abs(cumulative(f,:total)-cumulative(f,:direct))
        endpointgap=abs(cumulative(f,:direct)-log(s.Ndot[end]/s.Ndot[1]))
        stockannual=maximum(abs.(a.resid[2:end]))
        stockcum=abs(cumulative(a,:total)-cumulative(a,:direct))
        lng=log.(s.Ndot./s.N)
        denominator=1-GAMMA*s.drow[j]
        bridgegap=maximum(2:ND) do t
            maximum(abs.((
                a.productivity[t]-f.productivity[t]/denominator,
                a.labor[t]-f.employment[t]/denominator,
                a.external[t]-f.other_ideas[t]/denominator,
                a.transition[t]+(lng[t]-lng[t-1])/denominator,
                f.direct[t]-a.direct[t]-(lng[t]-lng[t-1]),
            )))
        end
        @printf("  %-14s flow level=%7.1e annual=%7.1e cumulative=%7.1e endpoint=%7.1e\n",
                s.name,levelgap,annualgap,cumgap,endpointgap)
        @printf("  %-14s appendix annual=%7.1e cumulative=%7.1e bridge=%7.1e\n",
            "",stockannual,stockcum,bridgegap)
        levelgap<=1e-10 || error("Flow level identity failed for $(s.name): $levelgap")
        annualgap<=1e-10 || error("Flow annual identity failed for $(s.name): $annualgap")
        cumgap<=1e-12 || error("Flow cumulative identity failed for $(s.name): $cumgap")
        endpointgap<=1e-12 || error("Flow endpoint identity failed for $(s.name): $endpointgap")
        stockannual<=1e-10 || error("Stock annual identity failed for $(s.name): $stockannual")
        stockcum<=1e-10 || error("Stock cumulative identity failed for $(s.name): $stockcum")
        bridgegap<=1e-12 || error("Flow-to-stock bridge failed for $(s.name): $bridgegap")
    end
end

function print_console_table(title,channels,decomps,totalfield)
    println("\n"*"="^82)
    println("  $title")
    println("  Entries are cumulative contributions in log points")
    println("="^82)
    @printf("  %-47s %9s %12s %13s\n","Channel","Goods","RI-Services","NRI-Services")
    for (label,field) in channels
        @printf("  %-47s %9.4f %12.4f %13.4f\n",label,
                (cumulative(d,field) for d in decomps)...)
    end
    @printf("  %-47s %9.4f %12.4f %13.4f\n",totalfield,
            (cumulative(d,:direct) for d in decomps)...)
end

function latex_helpers(io)
    fmt(x,digits)=begin
        r=round(x,digits=digits)
        abs(r)<0.5*10.0^(-digits) ? @sprintf("%.*f",digits,0.0) : @sprintf("%.*f",digits,r)
    end
    row(label,vals;digits=2)=println(io,"    "*label*" & "*
        join((fmt(x,digits) for x in vals)," & ")*" \\\\")
    memopair(label,vals;digits=2)=println(io,"    "*label*" & "*
        join(("\\multicolumn{1}{c}{"*fmt(v[1],digits)*" / "*fmt(v[2],digits)*"}" for v in vals)," & ")*" \\\\")
    row,memopair
end

# Channel rows rounded so that each column adds to its printed total (largest remainder).
function printed_channels(channels,decomps;digits=2)
    tot=[round(cumulative(d,:direct),digits=digits) for d in decomps]
    vals=zeros(length(channels),length(decomps))
    for (j,d) in enumerate(decomps)
        vals[:,j]=round_to_total([cumulative(d,field) for (_,field) in channels],tot[j];digits=digits)
    end
    tot,vals
end

function write_flow_latex(path::String)
    io=IOBuffer(); row,memopair=latex_helpers(io)
    println(io,"% Auto-generated by Exercise2.jl — do not edit by hand.")
    println(io,"% Requires: booktabs, threeparttable, siunitx, ragged2e, graphicx")
    println(io,"\\begin{table}[htbp]")
    println(io,"  \\centering")
    println(io,"  \\caption{Decomposition of the Growth in the Flow of New Ideas, 1950--2023}")
    println(io,"  \\label{tab:ex2_flow}")
    println(io,"  \\begin{threeparttable}")
    println(io,"  \\centering")
    println(io,"  \\resizebox{\\textwidth}{!}{%")
    println(io,"  \\begin{tabular}{l *{3}{S[table-format=-1.2]}}")
    println(io,"    \\toprule")
    println(io,"     & \\multicolumn{3}{c}{Sector} \\\\")
    println(io,"    \\cmidrule(lr){2-4}")
    println(io,"     & {Goods} & {RI-Services} & {NRI-Services} \\\\")
    println(io,"    \\midrule")
    channels=FLOW_CHANNELS
    tot,vals=printed_channels(channels,FLOW_DECOMPS)
    for (i,(label,_)) in enumerate(channels)
        row(replace(label,". "=>".\\ "),vals[i,:])
    end
    println(io,"    \\midrule")
    row("Total",tot)
    println(io,"    \\addlinespace")
    println(io,"    \\multicolumn{4}{l}{\\textit{Memo:}} \\\\")
    row("Flow of new ideas, 1950",[s.Ndot[1] for s in SECTORS];digits=3)
    row("Flow of new ideas, 2023",[s.Ndot[end] for s in SECTORS];digits=3)
    row("Cumulative stock growth, \$\\ln N_j(2023)-\\ln N_j(1950)\$",
        [cumulative(d,:direct) for d in STOCK_DECOMPS])
    println(io,"    \\bottomrule")
    println(io,"  \\end{tabular}}")
    println(io,"  \\end{threeparttable}")
    println(io,"  \\begin{minipage}{\\textwidth}")
    println(io,"    \\footnotesize")
    println(io,"    \\RaggedRight")
    println(io,"    \\noindent\\textit{Notes:} Rows I--IV report contributions, in log points, to the change")
    println(io,"    in the log flow of new varieties, \$\\ln\\dot N_j(2023)-\\ln\\dot N_j(1950)\$. The decomposition")
    println(io,"    follows Equation (20) and is exact before rounding; entries are rounded so that they add to the")
    println(io,"    printed totals. The memo rows report the initial and final flow of new varieties and cumulative")
    println(io,"    growth in the stock of varieties.")
    println(io,"  \\end{minipage}")
    println(io,"\\end{table}")
    open(path,"w") do f; write(f,String(take!(io))); end
    println("  Wrote main Exercise 2 table -> $path")
end

function write_stock_latex(path::String)
    io=IOBuffer(); row,memopair=latex_helpers(io)
    println(io,"% Auto-generated by Exercise2.jl — do not edit by hand.")
    println(io,"\\begin{table}[htbp]")
    println(io,"  \\centering")
    println(io,"  \\caption{Decomposition of Cumulative Idea-Stock Growth, 1950--2023}")
    println(io,"  \\label{tab:ex2_stock_appendix}")
    println(io,"  \\begin{threeparttable}")
    println(io,"  \\centering")
    println(io,"  \\begin{tabular}{l *{3}{S[table-format=-1.2]}}")
    println(io,"    \\toprule")
    println(io,"     & \\multicolumn{3}{c}{Sector} \\\\")
    println(io,"    \\cmidrule(lr){2-4}")
    println(io,"     & {Goods} & {RI-Services} & {NRI-Services} \\\\")
    println(io,"    \\midrule")
    channels=STOCK_CHANNELS
    tot,vals=printed_channels(channels,STOCK_DECOMPS)
    for (i,(label,_)) in enumerate(channels)
        row(replace(label,". "=>".\\ "),vals[i,:])
    end
    println(io,"    \\midrule")
    row("Total (\$\\Delta\\ln N_j\$)",tot)
    println(io,"    \\addlinespace")
    row("\\textit{Memo:} \$1/(1-\\gamma\\delta_{jj})\$",[multiplier(s) for s in SECTORS])
    memopair("\\textit{Memo:} \$g_{N_j}\$ (percent), 1950 / 2023",
             [(100s.Ndot[1]/s.N[1],100s.Ndot[end]/s.N[end]) for s in SECTORS])
    println(io,"    \\bottomrule")
    println(io,"  \\end{tabular}")
    println(io,"  \\begin{tablenotes}[flushleft]")
    println(io,"    \\footnotesize")
    println(io,"    \\RaggedRight")
    println(io,"    \\item[] \\textit{Notes:} Stock-based decomposition of Exercise 2. Rows I--IV sum exactly to \$\\Delta\\ln N_j\$")
    println(io,"    before rounding. The transition row captures changes in \$g_{N_j}=\\dot N_j/N_j\$.")
    println(io,"    Entries are rounded so that printed rows add to the printed totals.")
    println(io,"    RI- and NRI-Services denote research- and non-research-intensive services, respectively.")
    println(io,"  \\end{tablenotes}")
    println(io,"  \\end{threeparttable}")
    println(io,"\\end{table}")
    open(path,"w") do f; write(f,String(take!(io))); end
    println("  Wrote Exercise 2 appendix table -> $path")
end

function main()
    mkpath(OUTDIR)
    run_validation()
    print_console_table("Exercise 2 — growth in the flow of new varieties, 1950--2023",
                        FLOW_CHANNELS,FLOW_DECOMPS,"Total (Delta ln Ndot_j)")
    print_console_table("Exercise 2 appendix — cumulative idea-stock growth, 1950--2023",
                        STOCK_CHANNELS,STOCK_DECOMPS,"Total (Delta ln N_j)")
    write_flow_latex(joinpath(OUTDIR,"Table_4.tex"))
    write_stock_latex(joinpath(OUTDIR,"Exercise2_stock_appendix.tex"))
    println("\nExercise 2 complete.")
end

main()
