#==========================================================================================
 An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
 ------------------------------------------------------------------------------------------
 File:    Exercise5.jl
 Purpose: Exercise 5. Composition and aggregate labor productivity growth: aggregation and
          the innovation feedback (Section 5.5).
 ------------------------------------------------------------------------------------------
 Usage:   julia --project=. code/counterfactuals/Exercise5.jl
 ------------------------------------------------------------------------------------------
 Inputs:  output/calibration/final_calibration.txt, calibration_results.toml (via ../baseline.jl)
          Exercise3.jl and Exercise4.jl (included)
 Outputs: output/counterfactuals/Table_8.tex
          output/counterfactuals/Figure_H1.pdf
          output/counterfactuals/Exercise5_terms.csv, _feedback.csv, _perturbations.csv,
          _1980_dynamic_path.csv
 ------------------------------------------------------------------------------------------
 Notes:
 Compares the Exercise 3 baseline and full-shutdown endpoints, then performs
 dated demand perturbations under the baseline and symmetric-spillover specifications.
==========================================================================================#
let old=get(ENV,"EX3_NORUN",nothing)
    try
        ENV["EX3_NORUN"]="1"
        include(joinpath(@__DIR__,"Exercise3.jl"))
    finally
        old===nothing ? delete!(ENV,"EX3_NORUN") : (ENV["EX3_NORUN"]=old)
    end
end
let old=get(ENV,"EX4_NORUN",nothing)
    try
        ENV["EX4_NORUN"]="1"
        include(joinpath(@__DIR__,"Exercise4.jl"))
    finally
        old===nothing ? delete!(ENV,"EX4_NORUN") : (ENV["EX4_NORUN"]=old)
    end
end

using Printf
using Statistics
import Plots

const EX5_ND=EX3_ND
const EX5_OUTDIR=CAL_PATHS.out_counterfact
const EX5_YEARS=collect(1950:2023)
const EX5_DATES=(1950,1980,2000,2023)
const EX5_EPS=0.001
const EX5_ENDPOINTS=(:base,:ES)
_ex5sl(v)=collect(v[1:EX5_ND])

"""Solve only the full-shutdown endpoint; the baseline is already available."""
function ex5_solve_endpoints()
        p,e=ex3_inputs(:ES)
        full=run_counterfactual(p,e;V_seed=_ex3seed(eq_base),label="Exercise 5 full shutdown")
        (;eqs=Dict(:base=>eq_base,:ES=>full),
            exs=Dict(:base=>baseline_exog,:ES=>e))
end

function ex5_sector_rows(eq,e,j)
    keys=(("Nmt","pmt","Ymt","Lpmt","e_Y_mt",:zm,:hbarm,:muxm),
          ("Nrst","prst","Yrst","Lprst","e_Y_rst",:zrs,:hbarrs,:muxrs),
          ("Nnst","pnst","Ynst","Lpnst","e_Y_nst",:zns,:hbarns,:muxns))[j]
    N,p,Y,Lp,eY=(_ex5sl(eq[keys[k]]) for k in 1:5)
    z,h,mux=(_ex5sl(getproperty(e,keys[k])) for k in 6:8)
    c=(1-baseline_params.alpha)/baseline_params.alpha
    row1=zeros(EX5_ND);row2=zeros(EX5_ND);row3=zeros(EX5_ND);row4=zeros(EX5_ND);row5=zeros(EX5_ND)
    gA=zeros(EX5_ND)
    lnN=log.(N)
    for t in 2:EX5_ND
        row1[t]=eY[t]*(lnN[t]-lnN[t-1])
        row2[t]=c*log(z[t]/z[t-1])
        row3[t]=c*log(p[t]/p[t-1])
        row4[t]=log(h[t]/h[t-1])
        # Machine market power, including the drift of e_Y with mu^x; zero with constant markups.
        row5[t]=c*log(mux[t]/mux[t-1])+(eY[t]-eY[t-1])*lnN[t-1]
        gA[t]=log((Y[t]/Lp[t])/(Y[t-1]/Lp[t-1]))
    end
    (;row1,row2,row3,row4,row5,gA)
end

"""Tornqvist aggregation identity per production worker (within + reallocation + value-added
conversion), with Exercise 1 rows.

The conversion term is the share-weighted change in the value-added factors omega_j. It is nonzero
when the machine markups vary (the baseline) and zero when they are held at their means
(MUX_MODE=mean); the check total == direct in ex5_validate confirms the identity either way."""
function ex5_path_terms(eq,e)
    va=ex3_va(eq,e); rows=ntuple(j->ex5_sector_rows(eq,e,j),3)
    Lp=(_ex5sl(eq["Lpmt"]),_ex5sl(eq["Lprst"]),_ex5sl(eq["Lpnst"]))
    sumLp=Lp[1].+Lp[2].+Lp[3]
    weights=ntuple(_->zeros(EX5_ND),3)
    within=zeros(EX5_ND);realloc=zeros(EX5_ND);conv=zeros(EX5_ND)
    total=zeros(EX5_ND);direct=zeros(EX5_ND);benchmark=zeros(EX5_ND)
    sub=ntuple(_->zeros(EX5_ND),5)
    for t in 2:EX5_ND
        for j in 1:3; weights[j][t]=0.5*(va.shares[j][t]+va.shares[j][t-1]); end
        gLp=log(sumLp[t]/sumLp[t-1])
        within[t]=sum(weights[j][t]*rows[j].gA[t] for j in 1:3)
        realloc[t]=sum(weights[j][t]*(log(Lp[j][t]/Lp[j][t-1])-gLp) for j in 1:3)
        conv[t]=sum(weights[j][t]*log(va.omega[j][t]/va.omega[j][t-1]) for j in 1:3)
        total[t]=within[t]+realloc[t]+conv[t]
        direct[t]=sum(weights[j][t]*log(va.real[j][t]/va.real[j][t-1]) for j in 1:3)-gLp
        benchmark[t]=sum(va.shares[j][1]*rows[j].gA[t] for j in 1:3)
        for r in 1:5; sub[r][t]=sum(weights[j][t]*getfield(rows[j],Symbol("row$r"))[t] for j in 1:3); end
    end
    (;within,realloc,conv,total,direct,benchmark,sub,weights,rows)
end

ex5_avg(v,i0=1,i1=EX5_ND)=100sum(@view(v[i0+1:i1]))/(i1-i0)

"""Exact base-minus-counterfactual product split of the within term."""
function ex5_difference(base,cf)
    keep=zeros(EX5_ND);move=zeros(EX5_ND)
    rowdiff=ntuple(_->zeros(EX5_ND),5)
    for t in 2:EX5_ND
        for j in 1:3
            sbar=0.5*(base.weights[j][t]+cf.weights[j][t])
            gbar=0.5*(base.rows[j].gA[t]+cf.rows[j].gA[t])
            keep[t]+=sbar*(base.rows[j].gA[t]-cf.rows[j].gA[t])
            move[t]+=gbar*(base.weights[j][t]-cf.weights[j][t])
            for r in 1:5
                rowdiff[r][t]+=sbar*(getfield(base.rows[j],Symbol("row$r"))[t]-getfield(cf.rows[j],Symbol("row$r"))[t])
            end
        end
    end
        (;keep,move,rowdiff,dwithin=base.within.-cf.within,drealloc=base.realloc.-cf.realloc,
            dconv=base.conv.-cf.conv,dtotal=base.total.-cf.total)
end

function ex5_validate(terms,cache)
    for c in EX5_ENDPOINTS
        gap=maximum(abs.(terms[c].total.-terms[c].direct));gap<=1e-11||error("Exercise 5 identity failed in $c: $gap")
        d=ex5_difference(terms[:base],terms[c])
        maximum(abs.(d.dwithin.-d.keep.-d.move))<=1e-11||error("within product split failed in $c")
        g3=ex3_aggregate_lp(cache.eqs[c],cache.exs[c],1,EX5_ND)
        abs(ex5_avg(terms[c].direct)-g3)<=1e-10 ||
            error("Exercise 3/5 aggregate mismatch in $c")
    end
    println("  production-worker aggregation identity and product split validated")
end

function ex5_demand_paths(date,reference)
    idx=date-1950+1;T=baseline_params.T_ge
    tm=fill(reference["θm"][1],T);tr=fill(reference["θrs"][1],T)
    tm[idx:end].-=EX5_EPS
    # Locks apply to normalized N. solve_GE maps every locked state into
    # N_j^Lambda=q_j(1950)N_j before evaluating the knowledge aggregator.
    Nm=fill(NaN,T);Nr=fill(NaN,T);Nn=fill(NaN,T)
    Nm[1:idx].=reference["Nmt"][1:idx];Nr[1:idx].=reference["Nrst"][1:idx];Nn[1:idx].=reference["Nnst"][1:idx]
    (;idx,tm,tr,Nm,Nr,Nn)
end

function ex5_perturbation(date,symmetric,reference)
    q=ex5_demand_paths(date,reference)
    p=symmetric ? ex4_symmetric_params() : baseline_params
    eq=run_counterfactual(p,baseline_exog;V_seed=_ex3seed(reference),label="Exercise 5 demand shock $date",
        theta_m_path=q.tm,theta_rs_path=q.tr,N_m_lock=q.Nm,N_rs_lock=q.Nr,N_ns_lock=q.Nn)
    i=q.idx
    pi(x)=x["Lsnst"][i]/x["Lst"][i]
    gn(x)=x["Nnsdot"][i]/x["Nnst"][i]
    A(x)=x["Ynst"][i]/x["Lpnst"][i]
    impact=(pi_ns=pi(eq)-pi(reference),gN_ns=gn(eq)-gn(reference),lnA_ns=log(A(eq)/A(reference)))
    (;date,idx=i,reference,eq,impact)
end

function ex5_run_perturbations()
    symmetric_reference=ex4_solve(ex4_symmetric_params(),"Exercise 5 symmetric reference")
    base=Dict(d=>ex5_perturbation(d,false,eq_base) for d in EX5_DATES)
    sym=Dict(d=>ex5_perturbation(d,true,symmetric_reference) for d in EX5_DATES)
    (;base,sym)
end

function ex5_validate_perturbations(pert)
    q=baseline_params.lambda_N0; seeds=(q.m,q.rs,q.ns)
    for set in (pert.base,pert.sym), date in EX5_DATES
        item=set[date]; idx=item.idx
        N=(item.eq["Nmt"],item.eq["Nrst"],item.eq["Nnst"])
        Nref=(item.reference["Nmt"],item.reference["Nrst"],item.reference["Nnst"])
        Nlambda=(item.eq["NΛmt"],item.eq["NΛrst"],item.eq["NΛnst"])
        maximum(maximum(abs.(N[j][1:idx].-Nref[j][1:idx])) for j in 1:3)<=1e-12 ||
            error("Exercise 5 normalized-N lock failed for $date")
        maximum(maximum(abs.(Nlambda[j].-seeds[j].*N[j])) for j in 1:3)<=1e-12 ||
            error("Exercise 5 Lambda refresh failed for $date")
    end
    return nothing
end

function ex5_cum(v)
    out=ones(EX5_ND)
    for t in 2:EX5_ND;out[t]=out[t-1]*exp(v[t]);end
    out
end

function ex5_write(terms,pert)
    mkpath(EX5_OUTDIR)
    fields=(:within,:realloc,:conv,:total)
    open(joinpath(EX5_OUTDIR,"Exercise5_terms.csv"),"w") do io
        println(io,"term,"*join((EX3_LABEL[c] for c in EX5_ENDPOINTS),","))
        for f in fields;println(io,string(f)*","*join((@sprintf("%.8f",ex5_avg(getfield(terms[c],f))) for c in EX5_ENDPOINTS),","));end
        for r in 1:5;println(io,"within_row$r,"*join((@sprintf("%.8f",ex5_avg(terms[c].sub[r])) for c in EX5_ENDPOINTS),","));end
    end
    open(joinpath(EX5_OUTDIR,"Exercise5_perturbations.csv"),"w") do io
        println(io,"date,specification,delta_pi_ns,delta_gN_ns,delta_lnA_ns")
        for d in EX5_DATES, (name,set) in (("baseline",pert.base),("symmetric",pert.sym))
            x=set[d].impact;@printf(io,"%d,%s,%.10f,%.10f,%.10f\n",d,name,x.pi_ns,x.gN_ns,x.lnA_ns)
        end
    end
    d=pert.base[1980];open(joinpath(EX5_OUTDIR,"Exercise5_1980_dynamic_path.csv"),"w") do io
        println(io,"year,delta_pi_ns,delta_gN_ns,delta_lnA_ns")
        for i in d.idx:EX5_ND
            pi(x)=x["Lsnst"][i]/x["Lst"][i];gn(x)=x["Nnsdot"][i]/x["Nnst"][i];A(x)=x["Ynst"][i]/x["Lpnst"][i]
            @printf(io,"%d,%.10f,%.10f,%.10f\n",i+1949,pi(d.eq)-pi(d.reference),gn(d.eq)-gn(d.reference),log(A(d.eq)/A(d.reference)))
        end
    end
    b=ex5_cum(terms[:base].direct);arith=ex5_cum(terms[:base].benchmark);full=ex5_cum(terms[:ES].direct)
    plt=Plots.plot(EX5_YEARS,b,label="Baseline",lw=2,color=:black,
                   ylims=(1,5),yticks=(1:5,string.(1:5)),
                   xlims=(first(EX5_YEARS),last(EX5_YEARS)),xlabel="Year",
                   ylabel="Aggregate Labor Productivity Growth\n(1950 = 1.00)")
    Plots.plot!(plt,EX5_YEARS,arith,label="Fixed-share benchmark",ls=:dash,lw=2)
    Plots.plot!(plt,EX5_YEARS,full,label="Full shutdown",ls=:dashdot,lw=2)
    Plots.savefig(plt,joinpath(EX5_OUTDIR,"Figure_H1.pdf"))
end

function ex5_write_latex(path,terms,d,X,Y,Z)
    # Table 8 reports cumulative log points over 1950-2023: annual averages times 73 years.
    K=EX5_ND-1
    io=IOBuffer()
    fmt1(x)=begin
        r=round(x,digits=1)
        abs(r)<0.05 ? "0.0" : @sprintf("%.1f",r)
    end
    triple(label,b,f,delta)=println(io,"    "*label*" & "*fmt1(b)*" & "*fmt1(f)*" & "*fmt1(delta)*" "*"\\\\")
    effect(label,x)=println(io,"    "*label*" & {--} & {--} & "*fmt1(x)*" "*"\\\\")
    cum(x)=K*ex5_avg(x)

    base=terms[:base];full=terms[:ES]
    println(io,"% Auto-generated by Exercise5.jl — do not edit by hand.")
    println(io,"\\begin{table}[htbp]")
    println(io,"  \\centering")
    println(io,"  \\caption{Composition and Cumulative Aggregate Labor Productivity Growth}")
    println(io,"  \\label{tab:ex5_decomposition}")
    println(io,"  \\begin{threeparttable}")
    println(io,"  \\centering")
    println(io,"  \\begin{tabular}{l S[table-format=-3.1] S[table-format=-3.1] S[table-format=-2.1]}")
    println(io,"    \\toprule")
    println(io,"     & {Baseline} & {Full shutdown} & {Difference} "*"\\\\")
    println(io,"    \\midrule")
    # Panel C drops rows that are identically zero (embodied technical progress and labor
    # efficiency follow the same exogenous paths in both equilibria). Entries are rounded so that
    # rows add to the printed totals.
    nz(x)=abs(x)>5e-7
    rowsA=[("I.\\ Within-sector productivity",base.within,full.within,d.dwithin),
           ("II.\\ Labor reallocation",base.realloc,full.realloc,d.drealloc)]
    # The value-added conversion row appears only when the machine markups vary.
    any(nz(ex5_avg(x)) for x in (base.conv,full.conv)) &&
        push!(rowsA,("III.\\ Value-added conversion",base.conv,full.conv,d.dconv))
    totA=(round(cum(base.total),digits=1),round(cum(full.total),digits=1),round(K*X,digits=1))
    colA=[round_to_total([cum(r[k+1]) for r in rowsA],totA[k];digits=1) for k in 1:3]
    println(io,"    \\multicolumn{4}{l}{\\textit{Panel A. Cumulative aggregate productivity growth}} "*"\\\\")
    for (i,r) in enumerate(rowsA)
        triple(r[1],colA[1][i],colA[2][i],colA[3][i])
    end
    println(io,"    \\midrule")
    triple("Total",totA...)
    println(io,"    \\addlinespace")
    println(io,"    \\multicolumn{4}{l}{\\textit{Panel B. Decomposition of the total difference}} "*"\\\\")
    YZ=round_to_total([K*Y,K*Z],totA[3];digits=1)
    effect("I.\\ Accounting component",YZ[1])
    effect("II.\\ Productivity feedback",YZ[2])
    println(io,"    \\midrule")
    effect("Total",totA[3])
    println(io,"    \\addlinespace")
    println(io,"    \\multicolumn{4}{l}{\\textit{Panel C. Decomposition of productivity feedback}} "*"\\\\")
    labels=("Endogenous innovation","Embodied technical progress","Relative-price effect","Labor efficiency",
            "Machine market power")
    keepC=[r for r in 1:5 if nz(ex5_avg(d.rowdiff[r]))]
    valsC=round_to_total([cum(d.rowdiff[r]) for r in keepC],YZ[2];digits=1)
    for (i,r) in enumerate(keepC)
        effect("\\hspace{1em}($('a'+i-1)) $(labels[r])",valsC[i])
    end
    println(io,"    \\midrule")
    effect("Total",YZ[2])
    println(io,"    \\bottomrule")
    println(io,"  \\end{tabular}")
    println(io,"  \\begin{tablenotes}[flushleft]")
    println(io,"    \\footnotesize")
    println(io,"    \\RaggedRight")
    println(io,"    \\item[] \\textit{Notes:} Entries are cumulative log-point contributions over 1950--2023.")
    println(io,"  Difference is full shutdown minus baseline. Panel A follows Equation (23). Panel B separates")
    println(io,"  changes in sectoral productivity growth from reweighting using Equation (24); the accounting")
    println(io,"  component combines reweighting with the change in the labor-reallocation term. Panel C applies")
    println(io,"  the sectoral productivity decomposition in Equation (18) to the productivity-feedback component;")
    println(io,"  embodied technical progress and labor efficiency follow the same exogenous paths in both")
    println(io,"  equilibria and contribute nothing. Entries are rounded so that they add to the printed totals.")
    println(io,"  Differences and totals are calculated before rounding.")
    println(io,"  \\end{tablenotes}")
    println(io,"  \\end{threeparttable}")
    println(io,"\\end{table}")

    open(path,"w") do f
        write(f,String(take!(io)))
    end
    println("  Wrote LaTeX table -> $path")
end

function main_exercise5()
    println("="^72*"\n  Exercise 5 — aggregation and innovation feedback\n"*"="^72)
    cache=ex5_solve_endpoints();terms=Dict(c=>ex5_path_terms(cache.eqs[c],cache.exs[c]) for c in EX5_ENDPOINTS)
    ex5_validate(terms,cache)
    d=ex5_difference(terms[:ES],terms[:base])
    X=ex5_avg(d.dtotal);Z=ex5_avg(d.keep);Y=X-Z
    @printf("  composition effect X=%+.4f pp/yr; accounting Y=%+.4f; feedback Z=%+.4f\n",X,Y,Z)
    @printf("  feedback split: endogenous innovation=%+.6f, embodied=%+.6f, relative-price=%+.6f, labor=%+.6f, machine market power=%+.6f\n",
            (ex5_avg(d.rowdiff[r]) for r in 1:5)...)
    @printf("  value-added conversion: baseline=%+.6f, full shutdown=%+.6f\n",
            ex5_avg(terms[:base].conv),ex5_avg(terms[:ES].conv))
    pert=ex5_run_perturbations();ex5_validate_perturbations(pert);ex5_write(terms,pert)
    # Panel B of Table 8 in pp a year (read by robustness.jl for Table J.1).
    open(joinpath(EX5_OUTDIR,"Exercise5_feedback.csv"),"w") do io
        println(io,"gain,accounting,feedback");@printf(io,"%.10f,%.10f,%.10f\n",X,Y,Z)
    end
    ex5_write_latex(joinpath(EX5_OUTDIR,"Table_8.tex"),terms,d,X,Y,Z)
    (;cache,terms,perturbations=pert,X,Y,Z)
end

get(ENV,"EX5_NORUN","0")!="1"&&main_exercise5()
