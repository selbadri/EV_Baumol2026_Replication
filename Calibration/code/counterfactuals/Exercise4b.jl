#==========================================================================================
 An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
 ------------------------------------------------------------------------------------------
 File:    Exercise4b.jl
 Purpose: Exercise 4b. Goods-anchored variant of Exercise 4 (robustness; not in the paper).
 ------------------------------------------------------------------------------------------
 Usage:   julia --project=. code/counterfactuals/Exercise4b.jl
 ------------------------------------------------------------------------------------------
 Inputs:  output/calibration/final_calibration.txt, calibration_results.toml (via ../baseline.jl)
          Exercise4.jl (included)
 Outputs: output/counterfactuals/Exercise4b_shapley.tex, _cells.csv, _shapley.csv, _spread.csv
          output/counterfactuals/Exercise4b_research_shares.pdf
 ------------------------------------------------------------------------------------------
 Notes:
 Identical to the three-player {D,L,G} Shapley of Exercise 4 in every respect
 EXCEPT the two both-common coalitions {L,G} and {D,L,G}. There, the frozen
 channel is pinned to the GOODS (m-sector) primitive instead of its own
 single-channel R*-solving value:
   L->G : zeta0_j = zeta0_m (Goods level), solve common g for R*.
   G->L : g_j    = g_m      (Goods growth), solve common zeta0 for R*.
 Both orders are computed and averaged, with the endpoint spread reported.
 The empty, D, and single-channel {L}/{G}/{D,L}/{D,G} cells are unchanged
 from Exercise 4 (the single-channel cells have nothing to freeze, so their
 one active channel is the R* solver). All Exercise 4 machinery — outcomes,
 coalitions, Shapley weights, validation, CSV/figure writer — is reused; only
 the coalition-solving step differs.
==========================================================================================#
if !isdefined(@__MODULE__, :solve_exercise4)
    let old=get(ENV,"EX4_NORUN",nothing)
        try
            ENV["EX4_NORUN"]="1"
            include(joinpath(@__DIR__,"Exercise4.jl"))
        finally
            old===nothing ? delete!(ENV,"EX4_NORUN") : (ENV["EX4_NORUN"]=old)
        end
    end
end

using Printf

const EX4B_CACHE = Ref{Any}(nothing)
const EX4B_GOODS_Z0 = baseline_params.zetam_ini
const EX4B_GOODS_G  = baseline_params.gzetam

function solve_exercise4b(;force=false)
    !force && EX4B_CACHE[]!==nothing && return EX4B_CACHE[]
    pbase=baseline_params; psym=ex4_symmetric_params()
    Deq=ex4_solve(psym,"D symmetric spillovers")
    # Single-channel cells: identical to Exercise 4 (one channel, so it solves R*).
    cL,eqL  =ex4_root(c->ex4_L(pbase,c),EX4_Z0B,"L";seed=baseline_context.V_seed)
    cG,eqG  =ex4_root(c->ex4_G(pbase,c),EX4_GB,"G";seed=baseline_context.V_seed)
    cDL,eqDL=ex4_root(c->ex4_L(psym,c),EX4_Z0B,"D,L";seed=_ex4seed(Deq))
    cDG,eqDG=ex4_root(c->ex4_G(psym,c),EX4_GB,"D,G";seed=_ex4seed(Deq))
    # Both-common cells: freeze the FIRST channel at GOODS, solve the second (both orders).
    z0=EX4B_GOODS_Z0; gg=EX4B_GOODS_G
    gLG,eqLG1  =ex4_root(c->ex4_G(ex4_L(pbase,z0),c),EX4_GB,"L,G (L->G)";seed=_ex4seed(eqL))
    zLG,eqLG2  =ex4_root(c->ex4_L(ex4_G(pbase,gg),c),EX4_Z0B,"L,G (G->L)";seed=_ex4seed(eqG))
    gDLG,eqDLG1=ex4_root(c->ex4_G(ex4_L(psym,z0),c),EX4_GB,"D,L,G (L->G)";seed=_ex4seed(eqDL))
    zDLG,eqDLG2=ex4_root(c->ex4_L(ex4_G(psym,gg),c),EX4_Z0B,"D,L,G (G->L)";seed=_ex4seed(eqDG))
    result=(D=(eq=Deq,),
            L=(eq=eqL,z0=cL), G=(eq=eqG,g=cG),
            DL=(eq=eqDL,z0=cDL), DG=(eq=eqDG,g=cDG),
            LG=(LtoG=(eq=eqLG1,z0=z0,g=gLG), GtoL=(eq=eqLG2,g=gg,z0=zLG)),
            DLG=(LtoG=(eq=eqDLG1,z0=z0,g=gDLG), GtoL=(eq=eqDLG2,g=gg,z0=zDLG)))
    ex4_validate(result); EX4B_CACHE[]=result
end

function main_exercise4b()
    println("="^72*"\n  Exercise 4b — {D,L,G} Shapley, Goods-anchored both-common freeze\n"*"="^72)
    r=solve_exercise4b(); ex4_write(r; prefix="Exercise4b")
    @printf("  target idea ratio %.8f\n",EX4_TARGET)
    @printf("  Goods anchors: zeta0=%.6f  g=%+.6f\n",EX4B_GOODS_Z0,EX4B_GOODS_G)
    @printf("  D-cell idea ratio %.8f (option b: D moves aggregate ideas)\n",ex4_outcomes(r.D.eq).ratio)
    @printf("  single-channel  L: zeta0=%.6f | G: g=%+.6f (baseline spillovers)\n",r.L.z0,r.G.g)
    @printf("  single-channel  D,L: zeta0=%.6f | D,G: g=%+.6f (symmetric spillovers)\n",r.DL.z0,r.DG.g)
    @printf("  L,G   (L->G) freeze zeta0=Goods, solve g=%+.6f | (G->L) freeze g=Goods, solve zeta0=%.6f\n",r.LG.LtoG.g,r.LG.GtoL.z0)
    @printf("  D,L,G (L->G) freeze zeta0=Goods, solve g=%+.6f | (G->L) freeze g=Goods, solve zeta0=%.6f\n",r.DLG.LtoG.g,r.DLG.GtoL.z0)
    r
end

get(ENV,"EX4B_NORUN","0")!="1" && main_exercise4b()
