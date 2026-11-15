#==========================================================================================
 An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
 ------------------------------------------------------------------------------------------
 File:    Exercise4.jl
 Purpose: Exercise 4. Homogenizing the innovation technology: three-player {D,L,G} Shapley
          decomposition (Section 5.4).
 ------------------------------------------------------------------------------------------
 Usage:   julia --project=. code/counterfactuals/Exercise4.jl
 ------------------------------------------------------------------------------------------
 Inputs:  output/calibration/final_calibration.txt, calibration_results.toml (via ../baseline.jl)
 Outputs: output/counterfactuals/Table_7.tex
          output/counterfactuals/Exercise4_shapley.tex   Shapley decomposition (not in the paper)
          output/counterfactuals/Exercise4_cells.csv, _shapley.csv, _spread.csv
          output/counterfactuals/Exercise4_research_shares.pdf   (not in the paper)
 ------------------------------------------------------------------------------------------
 Notes:
 Three-player {D,L,G} Shapley. D: symmetric spillover backdrop.
 L: common research-productivity LEVEL in 1950; each sector keeps its own growth path.
 G: common constant research-productivity GROWTH; each sector keeps its 1950 level.
    (All three paths are exponential in the baseline; ex4_L and ex4_G also handle knots.)
 In each L- or G-containing cell the shared scalar is solved so the aggregate
 idea ratio N(2023)/N(1950) equals baseline R*. Both-common cells {L,G} and
 {D,L,G} are order dependent: the last-homogenized channel carries the R*
 solve while the first is frozen at its single-channel value; the two orders
 are reported as their average with an endpoint spread. Option (b): idea
 growth is matched only in L- or G-containing cells; the empty and D-only
 cells keep baseline research productivity and let D move aggregate ideas.
 Baseline production wedges are fixed and theta is recovered to preserve
 baseline initial VA shares in every equilibrium.
==========================================================================================#
using Printf
using Statistics
import Plots

if !isdefined(@__MODULE__, :baseline_context)
    include(joinpath(@__DIR__, "..", "baseline.jl"))
end

const EX4_ND = baseline_params.T_data
const EX4_OUTDIR = CAL_PATHS.out_counterfact
const EX4_TOL = 2e-5
const EX4_CACHE = Ref{Any}(nothing)
_ex4sl(v)=collect(v[1:EX4_ND])
_ex4seed(eq)=(copy(eq["Vmt"]),copy(eq["Vrst"]),copy(eq["Vnst"]))

function ex4_symmetric_params(p=baseline_params)
    p.lambda_N0 == baseline_params.lambda_N0 ||
        error("Exercise 4 must hold the 1950 KPSS Lambda seeds fixed")
    diag=(p.dMM.+p.dRR.+p.dNN)./3
    off=(1 .- diag)./2
    with_params(p; dMM=copy(diag),dMR=copy(off),dMN=copy(off),
                   dRM=copy(off),dRR=copy(diag),dRN=copy(off),
                   dNM=copy(off),dNR=copy(off),dNN=copy(diag))
end

function ex4_va_shares(eq,e=baseline_exog)
    a=baseline_params.alpha
    omega=(1 .-(1-a).*_ex4sl(e.muym).*_ex4sl(e.muxm),
           1 .-(1-a).*_ex4sl(e.muyrs).*_ex4sl(e.muxrs),
           1 .-(1-a).*_ex4sl(e.muyns).*_ex4sl(e.muxns))
        Y=_ex4sl(eq["Yt"])
        nom=(omega[1].*_ex4sl(eq["smt"]).*Y,
            omega[2].*_ex4sl(eq["srst"]).*Y,
            omega[3].*_ex4sl(eq["snst"]).*Y)
    tot=nom[1].+nom[2].+nom[3]
    ntuple(j->nom[j]./tot,3)
end

"""Value-added Tornqvist aggregate normalized idea-stock ratio N(2023)/N(1950)."""
function ex4_idea_ratio(eq,e=baseline_exog)
    s=ex4_va_shares(eq,e)
    N=(_ex4sl(eq["Nmt"]),_ex4sl(eq["Nrst"]),_ex4sl(eq["Nnst"]))
    total=0.0
    for t in 2:EX4_ND
        total += sum(0.5*(s[j][t]+s[j][t-1])*log(N[j][t]/N[j][t-1]) for j in 1:3)
    end
    exp(total)
end
const EX4_TARGET = ex4_idea_ratio(eq_base)

# Baseline research-productivity anchors (Goods) and the sectoral extremes used
# as natural homogenization ranges. Root brackets are widened for robustness.
const EX4_ZETA0 = (baseline_params.zetam_ini, baseline_params.zetars_ini, baseline_params.zetans_ini)
const EX4_GZETA = (baseline_params.gzetam, baseline_params.gzetars, baseline_params.gzetans)

# Homogenization operators. In the exponential baseline (6 October 2026) L sets a common
# zeta_j,0 and G a common g_zeta,j. With services knots (BLOCKB_MODE=knots):
# L(c): common 1950 level c; each sector keeps its own growth path (Goods exponential, the
#   services' knot paths shifted so that zeta_j(1950) = c).
# G(c): common constant growth rate c; each sector keeps its 1950 level and its path becomes
#   exponential, which removes both the differences in mean growth and the curvature of the
#   services' paths. Applying both gives a common level and growth in either order.
function ex4_L(p, c)
    kn = p.zeta_knots === nothing ? nothing :
         (p.zeta_knots[1] .- p.zeta_knots[1][1] .+ log(c), p.zeta_knots[2] .- p.zeta_knots[2][1] .+ log(c))
    with_params(p; zetam_ini=c, zetars_ini=c, zetans_ini=c, zeta_knots=kn)
end
function ex4_G(p, c)
    p.zeta_knots === nothing && return with_params(p; gzetam=c, gzetars=c, gzetans=c)
    with_params(p; zetars_ini=exp(p.zeta_knots[1][1]), zetans_ini=exp(p.zeta_knots[2][1]),
                gzetam=c, gzetars=c, gzetans=c, zeta_knots=nothing)
end

function ex4_solve(p,label; seed=baseline_context.V_seed)
    run_counterfactual(p,baseline_exog;V_seed=seed,label="Exercise 4 $label")
end

"""Bracket and solve one scalar against the aggregate idea-stock target."""
function ex4_root(builder,bounds,label; seed=baseline_context.V_seed)
    cache=Dict{Float64,Tuple{Float64,Any}}()
    function evaluate(x)
        haskey(cache,x) && return cache[x]
        try
            eq=ex4_solve(builder(x),label;seed=seed)
            r=log(ex4_idea_ratio(eq)/EX4_TARGET)
            cache[x]=(r,eq)
        catch err
            @warn "Exercise 4 root evaluation failed" label x exception=(err,catch_backtrace())
            cache[x]=(NaN,nothing)
        end
    end
    grid=collect(range(bounds[1],bounds[2];length=21))
    vals=[evaluate(x) for x in grid]
    bracket=nothing
    for i in 1:length(grid)-1
        isfinite(vals[i][1])&&isfinite(vals[i+1][1])&&signbit(vals[i][1])!=signbit(vals[i+1][1]) &&
            (bracket=(grid[i],grid[i+1]);break)
    end
    # A root can sit between the last grid point that solves and the first that does not.
    # Subdivide each such interval before giving up; the ordinary search is unaffected.
    if bracket===nothing
        for i in 1:length(grid)-1
            isfinite(vals[i][1])==isfinite(vals[i+1][1]) && continue
            a=isfinite(vals[i][1]) ? grid[i] : grid[i+1]; fa=evaluate(a)[1]
            for x in range(a,isfinite(vals[i][1]) ? grid[i+1] : grid[i];length=12)[2:end-1]
                fx=evaluate(x)[1]
                isfinite(fx) || break
                if signbit(fx)!=signbit(fa)
                    bracket=a<x ? (a,x) : (x,a); break
                end
                a,fa=x,fx
            end
            bracket===nothing || break
        end
    end
    if bracket===nothing
        feasible=[(abs(vals[i][1]),grid[i],vals[i]) for i in eachindex(grid) if isfinite(vals[i][1])]
        isempty(feasible)&&error("Exercise 4 $label has no feasible root evaluations")
        best=feasible[argmin(first.(feasible))]
        if best[1]>EX4_TOL
            shape=join((@sprintf("%.4g:%+.4g",grid[i],vals[i][1]) for i in eachindex(grid)),"  ")
            error("Exercise 4 $label could not bracket target; best |log gap|=$(best[1])\n  grid x:logr = $shape")
        end
        return best[2],best[3][2]
    end
    lo,hi=bracket; flo,elo=evaluate(lo); fhi,ehi=evaluate(hi)
    best=abs(flo)<abs(fhi) ? (lo,flo,elo) : (hi,fhi,ehi)
    for _ in 1:40
        mid=(lo+hi)/2; fm,em=evaluate(mid)
        if isfinite(fm)
            abs(fm)<abs(best[2])&&(best=(mid,fm,em))
            abs(fm)<=EX4_TOL&&break
            if signbit(fm)==signbit(flo); lo,flo=mid,fm; else; hi,fhi=mid,fm; end
        else
            hi=mid
        end
    end
    abs(best[2])<=EX4_TOL || error("Exercise 4 $label root gap=$(best[2])")
    best[1],best[3]
end

function ex4_scientist_shares(eq)
    L=(_ex4sl(eq["Lsmt"]),_ex4sl(eq["Lsrst"]),_ex4sl(eq["Lsnst"]))
    tot=L[1].+L[2].+L[3]
    ntuple(j->L[j]./tot,3)
end

function ex4_outcomes(eq)
    pi=ex4_scientist_shares(eq)
    N=(_ex4sl(eq["Nmt"]),_ex4sl(eq["Nrst"]),_ex4sl(eq["Nnst"]))
    A=(_ex4sl(eq["Ymt"])./_ex4sl(eq["Lpmt"]),
       _ex4sl(eq["Yrst"])./_ex4sl(eq["Lprst"]),
       _ex4sl(eq["Ynst"])./_ex4sl(eq["Lpnst"]))
    gr(v)=100log(v[end]/v[1])/(EX4_ND-1)
    (;pi0=ntuple(j->pi[j][1],3),piT=ntuple(j->pi[j][end],3),
      gN=ntuple(j->gr(N[j]),3),gA=ntuple(j->gr(A[j]),3),ratio=ex4_idea_ratio(eq))
end

# Root brackets. On the symmetric backdrop aggregate ideas fall, so the common
# zeta0 that restores R* can be well above the sectoral extremes; the growth
# bracket stays modest since the solver diverges for large common g.
const EX4_Z0B = (0.02, 1.5)
const EX4_GB  = (-0.05, 0.045)

"""Average two outcome tuples field by field (order-averaged both-common cells)."""
_ex4avg(a,b)=(;pi0=ntuple(j->0.5*(a.pi0[j]+b.pi0[j]),3),
              piT=ntuple(j->0.5*(a.piT[j]+b.piT[j]),3),
              gN=ntuple(j->0.5*(a.gN[j]+b.gN[j]),3),
              gA=ntuple(j->0.5*(a.gA[j]+b.gA[j]),3),
              ratio=0.5*(a.ratio+b.ratio))

function ex4_validate(r)
    # Every coalition containing L or G must hit the aggregate idea target.
    rstar=(r.L.eq,r.G.eq,r.DL.eq,r.DG.eq,
           r.LG.LtoG.eq,r.LG.GtoL.eq,r.DLG.LtoG.eq,r.DLG.GtoL.eq)
    for eq in rstar
        gap=abs(log(ex4_idea_ratio(eq)/EX4_TARGET))
        gap<=EX4_TOL || error("Exercise 4 R* target gap=$gap")
    end
    q=baseline_params.lambda_N0; seeds=(q.m,q.rs,q.ns)
    for eq in (r.D.eq,rstar...)
        va0=ntuple(j->ex4_va_shares(eq)[j][1],3)
        maximum(abs.(collect(va0).-collect(baseline_identified.initial_va_shares)))<=1e-8 ||
            error("Exercise 4 theta recovery did not preserve baseline initial VA shares")
        maximum(abs.(eq["τprst"].-baseline_identified.tau_p_rs))<=1e-12 || error("RI wedge changed")
        maximum(abs.(eq["τpnst"].-baseline_identified.tau_p_ns))<=1e-12 || error("NRI wedge changed")
        N=(eq["Nmt"],eq["Nrst"],eq["Nnst"])
        Nlambda=(eq["NΛmt"],eq["NΛrst"],eq["NΛnst"])
        all(isapprox(N[j][1],1.0;atol=1e-12,rtol=1e-12) for j in 1:3) ||
            error("Exercise 4 normalized initial N changed")
        maximum(maximum(abs.(Nlambda[j].-seeds[j].*N[j])) for j in 1:3)<=1e-12 ||
            error("Exercise 4 Lambda seeds changed across scenarios")
    end
end

# Three-player {D,L,G} Shapley. D: symmetric spillover backdrop. L: common
# initial zeta level. G: common zeta growth. In every coalition containing L or
# G the last-homogenized research channel carries one free scalar solved so that
# N(2023)/N(1950) equals baseline; the empty and D-only cells keep baseline
# research productivity (option b), so D is free to move aggregate ideas.
function solve_exercise4(;force=false)
    !force&&EX4_CACHE[]!==nothing&&return EX4_CACHE[]
    pbase=baseline_params; psym=ex4_symmetric_params()
    Deq=ex4_solve(psym,"D symmetric spillovers")
    # Single-channel cells: homogenize one zeta primitive to a common value that hits R*.
    cL,eqL  =ex4_root(c->ex4_L(pbase,c),EX4_Z0B,"L";seed=baseline_context.V_seed)
    cG,eqG  =ex4_root(c->ex4_G(pbase,c),EX4_GB,"G";seed=baseline_context.V_seed)
    cDL,eqDL=ex4_root(c->ex4_L(psym,c),EX4_Z0B,"D,L";seed=_ex4seed(Deq))
    cDG,eqDG=ex4_root(c->ex4_G(psym,c),EX4_GB,"D,G";seed=_ex4seed(Deq))
    # Both-common cells: freeze the first channel's solved value, solve the second (both orders).
    gLG,eqLG1  =ex4_root(c->ex4_G(ex4_L(pbase,cL),c),EX4_GB,"L,G (L->G)";seed=_ex4seed(eqL))
    zLG,eqLG2  =ex4_root(c->ex4_L(ex4_G(pbase,cG),c),EX4_Z0B,"L,G (G->L)";seed=_ex4seed(eqG))
    gDLG,eqDLG1=ex4_root(c->ex4_G(ex4_L(psym,cDL),c),EX4_GB,"D,L,G (L->G)";seed=_ex4seed(eqDL))
    zDLG,eqDLG2=ex4_root(c->ex4_L(ex4_G(psym,cDG),c),EX4_Z0B,"D,L,G (G->L)";seed=_ex4seed(eqDG))
    result=(D=(eq=Deq,),
            L=(eq=eqL,z0=cL), G=(eq=eqG,g=cG),
            DL=(eq=eqDL,z0=cDL), DG=(eq=eqDG,g=cDG),
            LG=(LtoG=(eq=eqLG1,z0=cL,g=gLG), GtoL=(eq=eqLG2,g=cG,z0=zLG)),
            DLG=(LtoG=(eq=eqDLG1,z0=cDL,g=gDLG), GtoL=(eq=eqDLG2,g=cDG,z0=zDLG)))
    ex4_validate(result); EX4_CACHE[]=result
end

"""Coalition outcome vectors for the three-player Shapley (both-common cells order-averaged)."""
function ex4_coalitions(r)
    (empty=ex4_outcomes(eq_base), D=ex4_outcomes(r.D.eq),
     L=ex4_outcomes(r.L.eq), G=ex4_outcomes(r.G.eq),
     DL=ex4_outcomes(r.DL.eq), DG=ex4_outcomes(r.DG.eq),
     LG=_ex4avg(ex4_outcomes(r.LG.LtoG.eq),ex4_outcomes(r.LG.GtoL.eq)),
     DLG=_ex4avg(ex4_outcomes(r.DLG.LtoG.eq),ex4_outcomes(r.DLG.GtoL.eq)))
end

"""Three-player Shapley values of one reported object (field `sym`), by sector."""
function ex4_shapley(v,sym)
    f(x)=getfield(x,sym)
    v0=f(v.empty);vD=f(v.D);vL=f(v.L);vG=f(v.G);vDL=f(v.DL);vDG=f(v.DG);vLG=f(v.LG);vDLG=f(v.DLG)
    phiD=ntuple(j->(vD[j]-v0[j])/3+(vDL[j]-vL[j])/6+(vDG[j]-vG[j])/6+(vDLG[j]-vLG[j])/3,3)
    phiL=ntuple(j->(vL[j]-v0[j])/3+(vDL[j]-vD[j])/6+(vLG[j]-vG[j])/6+(vDLG[j]-vDG[j])/3,3)
    phiG=ntuple(j->(vG[j]-v0[j])/3+(vDG[j]-vD[j])/6+(vLG[j]-vL[j])/6+(vDLG[j]-vDL[j])/3,3)
    total=ntuple(j->vDLG[j]-v0[j],3)
    (;phiD,phiL,phiG,total)
end

"""Table 7: outcomes in the baseline and in each homogenization experiment."""
function ex4_write_table7(path,v)
    cols=(v.empty,v.D,v.L,v.G,v.DLG)
    sectors=("Goods","RI-Services","NRI-Services")
    io=IOBuffer()
    row(label,f,digits)=println(io,"    $label & "*
        join((Printf.format(Printf.Format("%.$(digits)f"),f(c)) for c in cols)," & ")*" \\\\")
    panel(title)=println(io,"    \\multicolumn{6}{l}{\\textit{$title}} \\\\")
    println(io,"% Auto-generated by Exercise4.jl — do not edit by hand.")
    println(io,"\\begin{table}[htbp]")
    println(io,"  \\centering")
    println(io,"  \\caption{The Direction of Research and the Innovation Technology}")
    println(io,"  \\label{tab:ex4_cells}")
    println(io,"  \\begin{threeparttable}")
    println(io,"  \\begin{tabular}{l *{5}{S[table-format=3.2]}}")
    println(io,"    \\toprule")
    println(io,"     & {Baseline} & {Symmetric} & {Common} & {Common} & {All} \\\\")
    println(io,"     & & {spillovers} & {\$\\zeta_0\$} & {\$g_\\zeta\$} & {three} \\\\")
    println(io,"    \\midrule")
    panel("Panel A. Share of researchers in 2023 (\\%)")
    for j in 1:3; row(sectors[j],c->100c.piT[j],1); end
    println(io,"    \\addlinespace")
    panel("Panel B. Average annual idea-stock growth (\\%)")
    for j in 1:3; row(sectors[j],c->c.gN[j],2); end
    println(io,"    \\addlinespace")
    panel("Panel C. Average annual labor productivity growth (\\%)")
    for j in 1:3; row(sectors[j],c->c.gA[j],2); end
    println(io,"    \\addlinespace")
    panel("Memo: Aggregate idea index, 2023 (1950 = 100)")
    row("",c->100c.ratio,1)
    println(io,"    \\bottomrule")
    println(io,"  \\end{tabular}")
    println(io,"  \\begin{tablenotes}[flushleft]")
    println(io,"    \\footnotesize")
    println(io,"    \\RaggedRight")
    println(io,"    \\item[] \\textit{Notes:} Growth rates are annualized log changes over 1950--2023. Each")
    println(io,"  single-channel experiment homogenizes the indicated feature of the innovation technology while")
    println(io,"  retaining the other two at baseline. In the common-\$\\zeta_0\$ and common-\$g_\\zeta\$ experiments,")
    println(io,"  the common parameter is chosen to match the baseline 2023 aggregate idea index. The")
    println(io,"  symmetric-spillover experiment does not impose this normalization. The final column imposes all")
    println(io,"  three restrictions and reports average outcomes across the two orders of homogenizing the initial")
    println(io,"  level and growth rate of research productivity. Research shares may not sum to 100 because of")
    println(io,"  rounding.")
    println(io,"  \\end{tablenotes}")
    println(io,"  \\end{threeparttable}")
    println(io,"\\end{table}")
    open(f->write(f,String(take!(io))),path,"w")
    println("  Wrote LaTeX table -> $path")
end

"""Shapley decomposition table (not in the paper)."""
function ex4_write_latex(path,v; goods_anchored=false)
    io=IOBuffer()
    fmt3(x)=fmt_fixed(x,3)
    # The three Shapley terms are rounded to add to the printed total (largest remainder).
    row(label,s,j)=begin
        tot=round(s.total[j],digits=3)
        parts=round_to_total([s.phiD[j],s.phiL[j],s.phiG[j]],tot;digits=3)
        println(io,"    "*label*" & "*join((fmt3(x) for x in (parts...,tot))," & ")*" \\\\")
    end
    sectors=("Goods","RI-Services","NRI-Services")
    spi=ex4_shapley(v,:piT); sN=ex4_shapley(v,:gN); sA=ex4_shapley(v,:gA)
    caption=goods_anchored ?
        "Homogenizing the Innovation Technology: Goods-Anchored Shapley Decomposition" :
        "Homogenizing the Innovation Technology: Shapley Decomposition"

    println(io,"% Auto-generated by Exercise4.jl — do not edit by hand.")
    println(io,"\\begin{table}[htbp]")
    println(io,"  \\centering")
    println(io,"  \\caption{"*caption*"}")
    println(io,"  \\label{tab:"*(goods_anchored ? "ex4b" : "ex4")*"_shapley}")
    println(io,"  \\begin{threeparttable}")
    println(io,"  \\centering")
    println(io,"  \\begin{tabular}{l *{4}{S[table-format=-1.3]}}")
    println(io,"    \\toprule")
    println(io,"     & {Spillovers} & {Initial productivity} & {Productivity growth} & {Total} \\\\")
    println(io,"     & {\$D\$} & {\$L\$} & {\$G\$} & {} \\\\")
    println(io,"    \\midrule")
    println(io,"    \\multicolumn{5}{l}{\\textit{Panel A. Terminal innovation employment shares}} \\\\")
    for j in 1:3; row(sectors[j],spi,j); end
    println(io,"    \\addlinespace")
    println(io,"    \\multicolumn{5}{l}{\\textit{Panel B. Idea growth}} \\\\")
    for j in 1:3; row(sectors[j],sN,j); end
    println(io,"    \\addlinespace")
    println(io,"    \\multicolumn{5}{l}{\\textit{Panel C. Labor productivity growth}} \\\\")
    for j in 1:3; row(sectors[j],sA,j); end
    println(io,"    \\bottomrule")
    println(io,"  \\end{tabular}")
    println(io,"  \\begin{tablenotes}[flushleft]")
    println(io,"    \\footnotesize")
    println(io,"    \\RaggedRight")
    println(io,"    \\item[] \\textit{Notes:} Entries are three-player Shapley contributions from homogenizing")
    println(io,"    the innovation technology. \$D\$ imposes symmetric knowledge spillovers;")
    println(io,"    \$L\$ imposes a common initial research-productivity level; and \$G\$ imposes a common")
    println(io,"    research-productivity growth rate. Panel A reports changes in terminal employment-share")
    println(io,"    fractions. Panels B and C report contributions to annualized growth rates in percentage")
    println(io,"    points over 1950--2023. Total is the fully homogenized outcome minus the baseline outcome")
    println(io,"    and equals the sum of the three Shapley contributions before rounding.")
    if goods_anchored
        println(io,"    In coalitions containing both \$L\$ and \$G\$, the first homogenized channel is fixed at")
        println(io,"    its Goods-sector baseline value and the second is chosen to match the baseline aggregate")
        println(io,"    idea-stock ratio; the two order-specific outcomes are averaged.")
    else
        println(io,"    In coalitions containing both \$L\$ and \$G\$, the first homogenized channel is fixed at")
        println(io,"    its corresponding single-channel solution and the second is chosen to match the baseline")
        println(io,"    aggregate idea-stock ratio; the two order-specific outcomes are averaged.")
    end
    println(io,"    RI- and NRI-Services denote research- and non-research-intensive services, respectively.")
    println(io,"  \\end{tablenotes}")
    println(io,"  \\end{threeparttable}")
    println(io,"\\end{table}")
    open(path,"w") do f
        write(f,String(take!(io)))
    end
    println("  Wrote LaTeX table -> $path")
end

function ex4_write(r; prefix="Exercise4")
    mkpath(EX4_OUTDIR)
    v=ex4_coalitions(r)
    fmt(x)=isnan(x) ? "" : @sprintf("%.10f",x)
    theta(eq)=(eq["θm"][1],eq["θrs"][1],eq["θns"][1])
    nanθ=(NaN,NaN,NaN)
    function row(io,cell,g,z0,o,th)
        vals=(g,z0,o.ratio,th...,o.pi0...,o.piT...,o.gN...,o.gA...)
        println(io,cell*","*join((fmt(x) for x in vals),","))
    end
    open(joinpath(EX4_OUTDIR,"$(prefix)_cells.csv"),"w") do io
        println(io,"cell,common_g,common_zeta0,idea_ratio,theta_m,theta_rs,theta_ns,pi0_m,pi0_rs,pi0_ns,piT_m,piT_rs,piT_ns,gN_m,gN_rs,gN_ns,gA_m,gA_rs,gA_ns")
        row(io,"empty",NaN,NaN,v.empty,theta(eq_base))
        row(io,"D",NaN,NaN,v.D,theta(r.D.eq))
        row(io,"L",NaN,r.L.z0,v.L,theta(r.L.eq))
        row(io,"G",r.G.g,NaN,v.G,theta(r.G.eq))
        row(io,"DL",NaN,r.DL.z0,v.DL,theta(r.DL.eq))
        row(io,"DG",r.DG.g,NaN,v.DG,theta(r.DG.eq))
        row(io,"LG_LtoG",r.LG.LtoG.g,r.LG.LtoG.z0,ex4_outcomes(r.LG.LtoG.eq),theta(r.LG.LtoG.eq))
        row(io,"LG_GtoL",r.LG.GtoL.g,r.LG.GtoL.z0,ex4_outcomes(r.LG.GtoL.eq),theta(r.LG.GtoL.eq))
        row(io,"LG_avg",NaN,NaN,v.LG,nanθ)
        row(io,"DLG_LtoG",r.DLG.LtoG.g,r.DLG.LtoG.z0,ex4_outcomes(r.DLG.LtoG.eq),theta(r.DLG.LtoG.eq))
        row(io,"DLG_GtoL",r.DLG.GtoL.g,r.DLG.GtoL.z0,ex4_outcomes(r.DLG.GtoL.eq),theta(r.DLG.GtoL.eq))
        row(io,"DLG_avg",NaN,NaN,v.DLG,nanθ)
    end
    open(joinpath(EX4_OUTDIR,"$(prefix)_shapley.csv"),"w") do io
        println(io,"object,sector,phi_D,phi_L,phi_G,total,addup_gap")
        sl=("m","rs","ns")
        for (sym,name) in ((:piT,"piT"),(:gN,"gN"),(:gA,"gA"))
            s=ex4_shapley(v,sym)
            for j in 1:3
                gap=s.phiD[j]+s.phiL[j]+s.phiG[j]-s.total[j]
                println(io,join((name,sl[j],fmt(s.phiD[j]),fmt(s.phiL[j]),fmt(s.phiG[j]),fmt(s.total[j]),fmt(gap)),","))
            end
        end
    end
    open(joinpath(EX4_OUTDIR,"$(prefix)_spread.csv"),"w") do io
        println(io,"cell,object,spread_m,spread_rs,spread_ns")
        for (cell,a,b) in (("LG",r.LG.LtoG.eq,r.LG.GtoL.eq),("DLG",r.DLG.LtoG.eq,r.DLG.GtoL.eq))
            oa=ex4_outcomes(a); ob=ex4_outcomes(b)
            for (sym,name) in ((:piT,"piT"),(:gN,"gN"),(:gA,"gA"))
                d=ntuple(j->getfield(oa,sym)[j]-getfield(ob,sym)[j],3)
                println(io,join((cell,name,fmt(d[1]),fmt(d[2]),fmt(d[3])),","))
            end
        end
    end
    prefix=="Exercise4" && ex4_write_table7(joinpath(EX4_OUTDIR,"Table_7.tex"),v)
    ex4_write_latex(joinpath(EX4_OUTDIR,"$(prefix)_shapley.tex"),v;
                    goods_anchored=(prefix=="Exercise4b"))
    years=1950:(1950+EX4_ND-1); panels=Plots.Plot[]
    avgpath(a,b,j)=0.5.*(ex4_scientist_shares(a)[j].+ex4_scientist_shares(b)[j])
    lines=[("Baseline",j->ex4_scientist_shares(eq_base)[j]),
           ("D",j->ex4_scientist_shares(r.D.eq)[j]),
           ("L",j->ex4_scientist_shares(r.L.eq)[j]),
           ("G",j->ex4_scientist_shares(r.G.eq)[j]),
           ("D,L,G",j->avgpath(r.DLG.LtoG.eq,r.DLG.GtoL.eq,j))]
    for j in 1:3
        q=Plots.plot(title=("Goods","RI-Serv.","NRI-Serv.")[j],xlabel="year",ylabel="share")
        for (name,f) in lines
            Plots.plot!(q,years,f(j),label=name,lw=1.5)
        end
        Plots.plot!(q,years,ex4_va_shares(eq_base)[j],label="Baseline VA",ls=:dash,color=:black)
        push!(panels,q)
    end
    Plots.savefig(Plots.plot(panels...;layout=(1,3),size=(1300,400)),joinpath(EX4_OUTDIR,"$(prefix)_research_shares.pdf"))
end

function main_exercise4()
    println("="^72*"\n  Exercise 4 — three-player {D,L,G} Shapley homogenization\n"*"="^72)
    r=solve_exercise4();ex4_write(r)
    @printf("  target idea ratio %.8f\n",EX4_TARGET)
    @printf("  D-cell idea ratio %.8f (option b: D moves aggregate ideas)\n",ex4_outcomes(r.D.eq).ratio)
    @printf("  L: common zeta0=%.6f | G: common g=%+.6f (baseline spillovers)\n",r.L.z0,r.G.g)
    @printf("  D,L: common zeta0=%.6f | D,G: common g=%+.6f (symmetric spillovers)\n",r.DL.z0,r.DG.g)
    r
end

get(ENV,"EX4_NORUN","0")!="1"&&main_exercise4()
