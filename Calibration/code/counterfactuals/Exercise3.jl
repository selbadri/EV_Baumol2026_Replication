#==========================================================================================
 An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
 ------------------------------------------------------------------------------------------
 File:    Exercise3.jl
 Purpose: Exercise 3. Structural change and aggregate labor productivity: remove income
          effects, then relative-price effects (Section 5.3, Appendix H.3).
 ------------------------------------------------------------------------------------------
 Usage:   julia --project=. code/counterfactuals/Exercise3.jl
 ------------------------------------------------------------------------------------------
 Inputs:  output/calibration/final_calibration.txt, calibration_results.toml (via ../baseline.jl)
 Outputs: output/counterfactuals/Table_5.tex, Table_6.tex
          output/counterfactuals/Table_H3.tex, Table_H4.tex
          output/counterfactuals/Exercise3_subperiods.tex   growth by subperiod (not in the paper)
          output/counterfactuals/Exercise3_results.csv, Exercise3_objects.csv
 ------------------------------------------------------------------------------------------
 Notes:
 Three equilibria are solved: the baseline, income effects removed (E: eta_j = 1) and the
 full shutdown (E+S: eta_j = 1 and sigma = 1). Setting sigma = 1 also removes income effects
 from expenditure shares, so E comes first and the order is unique. Income contribution =
 base - E; relative-price contribution = E - E+S; they add exactly to the total.
==========================================================================================#
using Printf
using Statistics

if !isdefined(@__MODULE__, :baseline_context)
    include(joinpath(@__DIR__, "..", "baseline.jl"))
end

const EX3_ND = baseline_params.T_data
const EX3_I1980 = 1980 - 1950 + 1
const EX3_OUTDIR = CAL_PATHS.out_counterfact
const EX3_PERIODS = [(1, EX3_ND), (1, EX3_I1980), (EX3_I1980, EX3_ND)]
const EX3_CELLS = [:base, :E, :ES]
const EX3_LABEL = Dict(:base=>"Baseline", :E=>"E", :ES=>"E+S")
const EX3_CACHE = Ref{Any}(nothing)

_ex3sl(v) = collect(v[1:EX3_ND])
_ex3seed(eq) = (copy(eq["Vmt"]), copy(eq["Vrst"]), copy(eq["Vnst"]))

function ex3_inputs(cell)
    p = cell in (:E,:ES) ? with_params(baseline_params; etam=1.0, etars=1.0, etans=1.0) : baseline_params
    p = cell == :ES ? with_params(p; sigma=1.0) : p
    return p, baseline_exog
end

"""Solve and cache the two counterfactual equilibria."""
function solve_all_cells(; force=false)
    !force && EX3_CACHE[] !== nothing && return EX3_CACHE[]
    eqs = Dict{Symbol,Any}(:base=>eq_base)
    exs = Dict{Symbol,Any}(:base=>baseline_exog)
    parent = Dict(:E=>:base, :ES=>:E)
    for cell in EX3_CELLS[2:end]
        p, e = ex3_inputs(cell)
        eqs[cell] = run_counterfactual(p, e; V_seed=_ex3seed(eqs[parent[cell]]),
                                       label="Exercise 3 $(EX3_LABEL[cell])")
        exs[cell] = e
        @printf("  %-6s converged (diff_V=%.3e)\n", EX3_LABEL[cell], eqs[cell]["diag_diff_V"][1])
    end
    EX3_CACHE[] = (; eqs, exs)
end

"""Consolidated sectoral real VA and nominal VA shares."""
function ex3_va(eq, e)
    a = baseline_params.alpha
    omega = (1 .- (1-a).*_ex3sl(e.muym).*_ex3sl(e.muxm),
             1 .- (1-a).*_ex3sl(e.muyrs).*_ex3sl(e.muxrs),
             1 .- (1-a).*_ex3sl(e.muyns).*_ex3sl(e.muxns))
    Y = (_ex3sl(eq["Ymt"]), _ex3sl(eq["Yrst"]), _ex3sl(eq["Ynst"]))
    expenditure = (_ex3sl(eq["smt"]), _ex3sl(eq["srst"]), _ex3sl(eq["snst"]))
    aggregate_Y = _ex3sl(eq["Yt"])
    real = ntuple(j -> omega[j].*Y[j], 3)
    nominal = ntuple(j -> omega[j].*expenditure[j].*aggregate_Y, 3)
    total = nominal[1] .+ nominal[2] .+ nominal[3]
    shares = ntuple(j -> nominal[j]./total, 3)
    (; omega, real, nominal, shares)
end

function ex3_aggregate_lp(eq, e, i0, i1)
    va = ex3_va(eq, e)
    Lp = _ex3sl(eq["Lpmt"]).+_ex3sl(eq["Lprst"]).+_ex3sl(eq["Lpnst"])
    total = 0.0
    for t in i0+1:i1
        total += sum(0.5*(va.shares[j][t]+va.shares[j][t-1]) *
                     log(va.real[j][t]/va.real[j][t-1]) for j in 1:3)
        total -= log(Lp[t]/Lp[t-1])
    end
    100total/(i1-i0)
end

function ex3_benchmark(i0, i1)
    va = ex3_va(eq_base, baseline_exog)
    A = (_ex3sl(eq_base["Ymt"])./_ex3sl(eq_base["Lpmt"]),
         _ex3sl(eq_base["Yrst"])./_ex3sl(eq_base["Lprst"]),
         _ex3sl(eq_base["Ynst"])./_ex3sl(eq_base["Lpnst"]))
    100sum(va.shares[j][1]*log(A[j][i1]/A[j][i0]) for j in 1:3)/(i1-i0)
end

function ex3_objects(eq, e)
    va = ex3_va(eq,e)
    Lp = (_ex3sl(eq["Lpmt"]),_ex3sl(eq["Lprst"]),_ex3sl(eq["Lpnst"]))
    Ls = (_ex3sl(eq["Lsmt"]),_ex3sl(eq["Lsrst"]),_ex3sl(eq["Lsnst"]))
    sp = Lp[1].+Lp[2].+Lp[3]; ss = Ls[1].+Ls[2].+Ls[3]
    A = (_ex3sl(eq["Ymt"])./Lp[1],_ex3sl(eq["Yrst"])./Lp[2],_ex3sl(eq["Ynst"])./Lp[3])
    # Production, LP, and reported idea growth use normalized N_j(1950)=1.
    N = (_ex3sl(eq["Nmt"]),_ex3sl(eq["Nrst"]),_ex3sl(eq["Nnst"]))
    gr(v)=100log(v[end]/v[1])/(EX3_ND-1)
    (; va0=ntuple(j->va.shares[j][1],3), vaT=ntuple(j->va.shares[j][end],3),
       emp0=ntuple(j->Lp[j][1]/sp[1],3), empT=ntuple(j->Lp[j][end]/sp[end],3),
       sci0=ntuple(j->Ls[j][1]/ss[1],3), sciT=ntuple(j->Ls[j][end]/ss[end],3),
       gA=ntuple(j->gr(A[j]),3), gN=ntuple(j->gr(N[j]),3),
       theta=(eq["θm"][1],eq["θrs"][1],eq["θns"][1]))
end

ex3_growth(eqs,exs) = Dict(c=>ntuple(k->ex3_aggregate_lp(eqs[c],exs[c],EX3_PERIODS[k]...),3)
                               for c in EX3_CELLS)

"""Income effects removed first, then relative prices."""
ex3_sequential(g) = (; E=collect(g[:base]).-collect(g[:E]), S=collect(g[:E]).-collect(g[:ES]))

function ex3_validate(eqs,exs)
    b0=ex3_va(eqs[:base],exs[:base]).shares
    q=baseline_params.lambda_N0; seeds=(q.m,q.rs,q.ns)
    for c in EX3_CELLS[2:end]
        s=ex3_va(eqs[c],exs[c]).shares
        maximum(abs(s[j][1]-b0[j][1]) for j in 1:3)<=1e-8 || error("theta recovery failed in $c")
        maximum(abs.(eqs[c]["τprst"].-baseline_identified.tau_p_rs))<=1e-12 || error("RI wedge changed")
        maximum(abs.(eqs[c]["τpnst"].-baseline_identified.tau_p_ns))<=1e-12 || error("NRI wedge changed")
        N=(eqs[c]["Nmt"],eqs[c]["Nrst"],eqs[c]["Nnst"])
        Nlambda=(eqs[c]["NΛmt"],eqs[c]["NΛrst"],eqs[c]["NΛnst"])
        all(isapprox(N[j][1],1.0;atol=1e-12,rtol=1e-12) for j in 1:3) ||
            error("normalized initial N changed in $c")
        maximum(maximum(abs.(Nlambda[j].-seeds[j].*N[j])) for j in 1:3)<=1e-12 ||
            error("Lambda-stock convention changed in $c")
    end
    # With eta_j = sigma = 1 and constant machine markups (MUX_MODE=mean) the value-added shares
    # stay at their 1950 values. The baseline markups vary and move the value-added factors
    # omega_j, so the shares move with them and the check does not apply.
    full=ex3_va(eqs[:ES],exs[:ES]).shares
    gap=maximum(maximum(abs.(full[j].-full[j][1])) for j in 1:3)
    MUX_MODE=="mean" && (gap<=1e-10 || error("full-shutdown VA shares move by $gap"))
    @printf("  full-shutdown VA-share change: %.3e\n",gap)
end

function ex3_write(g,b,eqs,exs)
    mkpath(EX3_OUTDIR)
    open(joinpath(EX3_OUTDIR,"Exercise3_results.csv"),"w") do io
        println(io,"cell,LP_1950_2023,LP_1950_1980,LP_1980_2023")
        @printf(io,"Arithmetic benchmark,%.8f,%.8f,%.8f\n",b...)
        for c in EX3_CELLS; @printf(io,"%s,%.8f,%.8f,%.8f\n",EX3_LABEL[c],g[c]...); end
        sq=ex3_sequential(g)
        @printf(io,"Income effects (base - E),%.8f,%.8f,%.8f\n",sq.E...)
        @printf(io,"Relative-price effects (E - E+S),%.8f,%.8f,%.8f\n",sq.S...)
    end
    open(joinpath(EX3_OUTDIR,"Exercise3_objects.csv"),"w") do io
        println(io,"cell,theta_m,theta_rs,theta_ns,VA0_m,VA0_rs,VA0_ns,VAT_m,VAT_rs,VAT_ns,gA_m,gA_rs,gA_ns,gN_m,gN_rs,gN_ns")
        for c in EX3_CELLS
            o=ex3_objects(eqs[c],exs[c]); vals=(o.theta...,o.va0...,o.vaT...,o.gA...,o.gN...)
            println(io,EX3_LABEL[c]*","*join((@sprintf("%.10f",x) for x in vals),","))
        end
    end
end

function ex3_write_latex(path, g, b)
    io = IOBuffer()
    fmt3(x) = begin
        r = round(x, digits=3)
        abs(r) < 0.0005 ? "0.000" : @sprintf("%.3f", r)
    end
    row(label, values) = println(io, "    " * label * " & " *
                                     join((fmt3(x) for x in values), " & ") * " \\\\")

    println(io, "% Auto-generated by Exercise3.jl — do not edit by hand.")
    println(io, "\\begin{table}[htbp]")
    println(io, "  \\centering")
    println(io, "  \\caption{Structural Change and Aggregate Labor Productivity Growth}")
    println(io, "  \\label{tab:ex3_channels_growth}")
    println(io, "  \\begin{threeparttable}")
    println(io, "  \\centering")
    println(io, "  \\begin{tabular}{l S[table-format=-1.3] S[table-format=-1.3] S[table-format=-1.3]}")
    println(io, "    \\toprule")
        println(io, "     & {1950--2023} & {1950--1980} & {1980--2023} \\\\")
    println(io, "    \\midrule")
    println(io, "    \\multicolumn{4}{l}{\\textit{Panel A. Annualized growth rates}} \\\\")
    row("Fixed-shares benchmark", b)
    row("Baseline model", g[:base])
    row("Income effects removed", g[:E])
    row("Full shutdown", g[:ES])
    println(io, "    \\addlinespace")
    println(io, "    \\multicolumn{4}{l}{\\textit{Panel B. Contributions to baseline growth}} \\\\")
    sq = ex3_sequential(g)
    tot = [round(x, digits=3) for x in sq.E .+ sq.S]
    parts = [round_to_total([sq.E[k], sq.S[k]], tot[k]; digits=3) for k in 1:3]
    row("I.\\ Income effects", [p[1] for p in parts])
    row("II.\\ Relative-price effects", [p[2] for p in parts])
    println(io, "    \\midrule")
    row("Total", tot)
    println(io, "    \\bottomrule")
    println(io, "  \\end{tabular}")
    println(io, "  \\begin{tablenotes}[flushleft]")
    println(io, "    \\footnotesize")
    println(io, "    \\RaggedRight")
    println(io, "    \\item[] \\textit{Notes:} Panel A reports annualized log growth rates of aggregate value added")
    println(io, "  per production worker, \$VA/L^p\$, in percentage points. The fixed-shares benchmark weights")
    println(io, "  baseline sectoral labor productivity growth by 1950 value-added shares. Income effects are")
    println(io, "  removed by setting \$\\eta_j=1\$; the full shutdown also sets \$\\sigma=1\$. Each row of")
    println(io, "  Panel A is a solved equilibrium. Panel B reports baseline growth minus growth with income effects removed")
    println(io, "  (I) and growth with income effects removed minus full-shutdown growth (II); a negative entry")
    println(io, "  means that the mechanism lowers aggregate labor productivity growth. Entries are rounded so")
    println(io, "  that they add to the printed totals.")
    println(io, "  \\end{tablenotes}")
    println(io, "  \\end{threeparttable}")
    println(io, "\\end{table}")

    open(path, "w") do f
        write(f, String(take!(io)))
    end
    println("  Wrote LaTeX table -> $path")
end

"""2023 and path objects of one solved equilibrium, as reported in Tables 6, H.3 and H.4."""
function ex3_report_objects(eqs,exs,g,c)
    eq=eqs[c];e=exs[c];va=ex3_va(eq,e)
    Lp=(_ex3sl(eq["Lpmt"]),_ex3sl(eq["Lprst"]),_ex3sl(eq["Lpnst"]))
    Ls=(_ex3sl(eq["Lsmt"]),_ex3sl(eq["Lsrst"]),_ex3sl(eq["Lsnst"]))
    p=(_ex3sl(eq["pmt"]),_ex3sl(eq["prst"]),_ex3sl(eq["pnst"]))
    N=(_ex3sl(eq["Nmt"]),_ex3sl(eq["Nrst"]),_ex3sl(eq["Nnst"]))
    A=(_ex3sl(eq["Ymt"])./Lp[1],_ex3sl(eq["Yrst"])./Lp[2],_ex3sl(eq["Ynst"])./Lp[3])
    lp_tot=Lp[1].+Lp[2].+Lp[3];ls_tot=Ls[1].+Ls[2].+Ls[3]
    ratio(v)=v[end]/v[1]
    growth(v)=100log(ratio(v))/(EX3_ND-1)
    (;va,prod=ntuple(j->Lp[j]./lp_tot,3),sci=ntuple(j->Ls[j]./ls_tot,3),
      agg_sci=ls_tot./_ex3sl(e.Lt),A,N,p,
      aggA=exp(g[c][1]*(EX3_ND-1)/100),
      Y=ratio(_ex3sl(eq["Yt"])),gY=growth(_ex3sl(eq["Yt"])))
end

function ex3_write_counterfactual_latex(paths,eqs,exs,g)
    cells=EX3_CELLS
    sectors=("Goods","RI-Services","NRI-Services")
    br="\\\\"
    fmt3(x)=begin
        r=round(x,digits=3)
        abs(r)<0.0005 ? "0.000" : @sprintf("%.3f",r)
    end
    fmt1(x)=begin
        r=round(x,digits=1)
        abs(r)<0.05 ? "0.0" : @sprintf("%.1f",r)
    end
    o=Dict(c=>ex3_report_objects(eqs,exs,g,c) for c in cells)
    reported_values(getter)=(getter(o[:base]),getter(o[:E]),getter(o[:ES]))
    function pairrow(io,label,getter)
        pairs=map((1,EX3_ND)) do t
            reported_values(x->getter(x,t))
        end
        values=reduce(vcat,([pairs[1][k],pairs[2][k]] for k in 1:3))
        println(io,label*" & "*join((fmt3(x) for x in values)," & ")*br)
    end
    indexrow(io,label,getter)=println(io,label*" & "*
        join((fmt3(x) for x in reported_values(getter))," & ")*br)
    growthrow(io,getter)=println(io,"\\hspace{1em}\\textit{Average annual growth (\\%)} & "*
        join((fmt1(x) for x in reported_values(getter))," & ")*br)

    io=IOBuffer()
    println(io,"% Auto-generated by Exercise3.jl — do not edit by hand.")
    println(io,"% Requires: booktabs, threeparttable, siunitx, ragged2e, graphicx, float")
    println(io,"\\begin{table}[H]")
    println(io,"\\centering")
    println(io,"\\caption{Structural Change Counterfactuals: Sectoral Allocations}")
    println(io,"\\label{tab:ex3_counterfactuals}")
    println(io,"\\begin{threeparttable}")
    println(io,"\\centering")
    println(io,"\\resizebox{\\textwidth}{!}{%")
    println(io,"\\begin{tabular}{l *{6}{S[table-format=-1.3]}}")
    println(io,"\\toprule")
    println(io," & \\multicolumn{2}{c}{Baseline} & \\multicolumn{2}{c}{Income effects removed} & \\multicolumn{2}{c}{Full shutdown} "*br)
    println(io,"\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7}")
    println(io," & {1950} & {2023} & {1950} & {2023} & {1950} & {2023} "*br)
    println(io,"\\midrule")
    println(io,"\\multicolumn{7}{l}{\\textit{Panel A. Sectoral value-added shares}} "*br)
    for j in 1:3;pairrow(io,sectors[j],(x,t)->x.va.shares[j][t]);end
    println(io,"\\addlinespace")
    println(io,"\\multicolumn{7}{l}{\\textit{Panel B. Sectoral innovation employment shares}} "*br)
    for j in 1:3;pairrow(io,sectors[j],(x,t)->x.sci[j][t]);end
    println(io,"\\addlinespace")
    println(io,"\\multicolumn{7}{l}{\\textit{Panel C. Sectoral production employment shares}} "*br)
    for j in 1:3;pairrow(io,sectors[j],(x,t)->x.prod[j][t]);end
    println(io,"\\bottomrule")
    println(io,"\\end{tabular}}")
    println(io,"\\end{threeparttable}")
    println(io,"\\begin{minipage}{\\textwidth}")
    println(io,"\\footnotesize")
    println(io,"\\RaggedRight")
    println(io,"\\noindent\\textit{Notes:} Each pair of columns reports a solved equilibrium in 1950 and 2023.")
    println(io,"Income effects are removed by setting \$\\eta_j=1\$; the full shutdown also sets \$\\sigma=1\$.")
    println(io,"In each counterfactual the expenditure weights are re-solved to reproduce the baseline 1950 value-added")
    println(io,"shares. Panel A reports sectoral value-added shares; Panels B and C report sectoral shares of")
    println(io,"innovation and production employment. The aggregate share of employment devoted to innovation is")
    println(io,"omitted because it follows the same externally calibrated path in all exercises.")
    println(io,"RI- and NRI-Services denote research- and non-research-intensive services, respectively.")
    println(io,"\\end{minipage}")
    println(io,"\\end{table}")
    open(paths[1],"w") do f
        write(f,String(take!(io)))
    end

    println(io,"% Auto-generated by Exercise3.jl — do not edit by hand.")
    println(io,"% Requires: booktabs, threeparttable, siunitx, ragged2e, graphicx, float")

    println(io,"\\begin{table}[H]")
    println(io,"\\centering")
    println(io,"\\caption{Structural Change Counterfactuals: Productivity, Innovation, Prices, and Output}")
    println(io,"\\label{tab:ex3_counterfactuals_cont}")
    println(io,"\\begin{threeparttable}")
    println(io,"\\centering")
    println(io,"\\begin{tabular}{l *{3}{S[table-format=-2.3]}}")
    println(io,"\\toprule")
    println(io," & {Baseline} & {Income effects removed} & {Full shutdown} "*br)
    println(io,"\\midrule")
    println(io,"\\multicolumn{4}{l}{\\textit{Panel D. Labor productivity indices (1950 = 1)}} "*br)
    for j in 1:3
        indexrow(io,sectors[j],x->x.A[j][end]/x.A[j][1]);growthrow(io,x->100log(x.A[j][end]/x.A[j][1])/(EX3_ND-1))
    end
    indexrow(io,"Aggregate",x->x.aggA);growthrow(io,x->100log(x.aggA)/(EX3_ND-1))
    println(io,"\\addlinespace")
    println(io,"\\multicolumn{4}{l}{\\textit{Panel E. Machine-variety indices (1950 = 1)}} "*br)
    for j in 1:3;indexrow(io,sectors[j],x->x.N[j][end]/x.N[j][1]);end
    println(io,"\\addlinespace")
    println(io,"\\multicolumn{4}{l}{\\textit{Panel F. Sectoral price indices (1950 = 1)}} "*br)
    for j in 1:3;indexrow(io,sectors[j],x->x.p[j][end]/x.p[j][1]);end
    println(io,"\\addlinespace")
    println(io,"\\multicolumn{4}{l}{\\textit{Panel G. Aggregate output index (1950 = 1)}} "*br)
    indexrow(io,"Output",x->x.Y);growthrow(io,x->x.gY)
    println(io,"\\bottomrule")
    println(io,"\\end{tabular}")
    println(io,"\\begin{tablenotes}[flushleft]")
    println(io,"\\footnotesize")
    println(io,"\\RaggedRight")
    println(io,"\\item[] \\textit{Notes:} Panels D--G report 2023 index levels, normalized to one in 1950, in")
    println(io,"the solved equilibria of the previous table. Average annual growth rates are annualized log")
    println(io,"changes in percent. Aggregate labor")
    println(io,"productivity is value added per production worker, \$VA/L^p\$, constructed using")
    println(io,"T\\\"{o}rnqvist value-added weights. Machine varieties are normalized knowledge stocks")
    println(io,"\$N_j\$. Sectoral prices and aggregate output are normalized by their 1950 values.")
    println(io,"RI- and NRI-Services denote research- and non-research-intensive services, respectively.")
    println(io,"\\end{tablenotes}")
    println(io,"\\end{threeparttable}")
    println(io,"\\end{table}")

    open(paths[2],"w") do f
        write(f,String(take!(io)))
    end
    println("  Wrote counterfactual LaTeX tables -> $(join(paths, ", "))")
end

"""Table 5: aggregate labor productivity growth, 1950-2023, and its 2023 index."""
function ex3_write_table5(path, g, b)
    rows = [("Baseline", g[:base][1]), ("Fixed-shares benchmark", b[1]), ("Full shutdown", g[:ES][1])]
    index(x) = 100exp(x*(EX3_ND-1)/100)
    io = IOBuffer()
    println(io, "% Auto-generated by Exercise3.jl — do not edit by hand.")
    println(io, "\\begin{table}[htbp]")
    println(io, "  \\centering")
    println(io, "  \\caption{Structural Change and Aggregate Labor Productivity}")
    println(io, "  \\label{tab:ex3_aggregate}")
    println(io, "  \\begin{threeparttable}")
    println(io, "  \\begin{tabular}{l S[table-format=1.2] S[table-format=3.1]}")
    println(io, "    \\toprule")
    println(io, "     & {Avg. annual growth (\\%)} & {Productivity index, 2023} \\\\")
    println(io, "     & {1950--2023} & {1950 = 100} \\\\")
    println(io, "    \\midrule")
    for (label, x) in rows
        @printf(io, "    %s & %.2f & %.1f \\\\\n", label, x, index(x))
    end
    println(io, "    \\bottomrule")
    println(io, "  \\end{tabular}")
    println(io, "  \\begin{tablenotes}[flushleft]")
    println(io, "    \\footnotesize")
    println(io, "    \\RaggedRight")
    println(io, "    \\item[] \\textit{Notes:} Aggregate labor productivity is real value added per production worker.")
    println(io, "  Average annual growth is the annualized log change over 1950--2023. Productivity indices are")
    println(io, "  normalized to 100 in 1950. The fixed-shares benchmark weights baseline sectoral labor")
    println(io, "  productivity growth by 1950 value added shares; its endpoint index is approximated from the")
    println(io, "  reported cumulative log growth. The full shutdown sets \$\\eta_j = 1\$ and \$\\sigma = 1\$, with the")
    println(io, "  equilibrium re-solved.")
    println(io, "  \\end{tablenotes}")
    println(io, "  \\end{threeparttable}")
    println(io, "\\end{table}")
    open(f -> write(f, String(take!(io))), path, "w")
    println("  Wrote LaTeX table -> $path")
end

"""Table 6: 2023 outcomes in the three solved equilibria."""
function ex3_write_table6(path, eqs, exs, g)
    o = Dict(c => ex3_report_objects(eqs, exs, g, c) for c in EX3_CELLS)
    sectors = ("Goods", "RI-Serv.", "NRI-Serv.")
    io = IOBuffer()
    row(label, f, fmt) = println(io, "    $label & " *
        join((Printf.format(Printf.Format(fmt), f(o[c])) for c in EX3_CELLS), " & ") * " \\\\")
    panel(title) = println(io, "    \\multicolumn{4}{l}{\\textit{$title}} \\\\")
    println(io, "% Auto-generated by Exercise3.jl — do not edit by hand.")
    println(io, "\\begin{table}[htbp]")
    println(io, "  \\centering")
    println(io, "  \\caption{Mechanisms Behind Structural Change and Labor Productivity}")
    println(io, "  \\label{tab:ex3_mechanisms}")
    println(io, "  \\begin{threeparttable}")
    println(io, "  \\begin{tabular}{l *{3}{S[table-format=3.3]}}")
    println(io, "    \\toprule")
    println(io, "     & {Baseline} & {Income effects removed} & {Full shutdown} \\\\")
    println(io, "    \\midrule")
    panel("Panel A. Aggregate labor productivity, 2023 (1950 = 100)")
    row("Aggregate", x -> 100x.aggA, "%.1f")
    println(io, "    \\addlinespace")
    panel("Panel B. Sectoral value-added shares, 2023 (\\%)")
    for j in 1:3; row(sectors[j], x -> 100x.va.shares[j][end], "%.1f"); end
    println(io, "    \\addlinespace")
    panel("Panel C. Sectoral scientist shares, 2023 (\\%)")
    for j in 1:3; row(sectors[j], x -> 100x.sci[j][end], "%.1f"); end
    println(io, "    \\addlinespace")
    panel("Panel D. Machine-variety indices, 2023 (1950 = 1)")
    for j in 1:3; row(sectors[j], x -> x.N[j][end]/x.N[j][1], "%.3f"); end
    println(io, "    \\addlinespace")
    panel("Panel E. Labor-productivity indices, 2023 (1950 = 1)")
    for j in 1:3; row(sectors[j], x -> x.A[j][end]/x.A[j][1], "%.3f"); end
    println(io, "    \\bottomrule")
    println(io, "  \\end{tabular}")
    println(io, "  \\begin{tablenotes}[flushleft]")
    println(io, "    \\footnotesize")
    println(io, "    \\RaggedRight")
    println(io, "    \\item[] \\textit{Notes:} The table reports 2023 outcomes in three solved equilibria. Income")
    println(io, "  effects are removed by setting \$\\eta_j = 1\$; the full shutdown also sets \$\\sigma = 1\$.")
    println(io, "  Appendix H.3 reports the sectoral outcomes in 1950 and 2023.")
    println(io, "  \\end{tablenotes}")
    println(io, "  \\end{threeparttable}")
    println(io, "\\end{table}")
    open(f -> write(f, String(take!(io))), path, "w")
    println("  Wrote LaTeX table -> $path")
end

function main_exercise3()
    println("="^72*"\n  Exercise 3 — income effects (E), then relative prices (S)\n"*"="^72)
    c=solve_all_cells(); g=ex3_growth(c.eqs,c.exs)
    b=ntuple(k->ex3_benchmark(EX3_PERIODS[k]...),3)
    ex3_validate(c.eqs,c.exs); ex3_write(g,b,c.eqs,c.exs)
    ex3_write_latex(joinpath(EX3_OUTDIR,"Exercise3_subperiods.tex"),g,b)
    ex3_write_table5(joinpath(EX3_OUTDIR,"Table_5.tex"),g,b)
    ex3_write_table6(joinpath(EX3_OUTDIR,"Table_6.tex"),c.eqs,c.exs,g)
    ex3_write_counterfactual_latex(joinpath.(EX3_OUTDIR,("Table_H3.tex","Table_H4.tex")),c.eqs,c.exs,g)
    @printf("  %-26s %10s %10s %10s\n","cell","1950-2023","1950-1980","1980-2023")
    @printf("  %-26s %10.4f %10.4f %10.4f\n","Arithmetic benchmark",b...)
    for x in EX3_CELLS; @printf("  %-26s %10.4f %10.4f %10.4f\n",EX3_LABEL[x],g[x]...); end
    sq=ex3_sequential(g)
    @printf("  %-26s %10.4f %10.4f %10.4f\n","Income effects",sq.E...)
    @printf("  %-26s %10.4f %10.4f %10.4f\n","Relative-price effects",sq.S...)
    (;c...,growth=g,benchmark=b,sequential=sq)
end

get(ENV,"EX3_NORUN","0")!="1" && main_exercise3()
