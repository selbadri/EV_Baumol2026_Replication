#==========================================================================================
 An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
 ------------------------------------------------------------------------------------------
 File:    equilibrium_path.jl
 Purpose: General-equilibrium solver. Solves the 100-period transition path of the
          three-sector model (goods, RI-services, NRI-services) and the patent-value fixed
          point for given parameters, in calibration or counterfactual mode. Included by
          smm_ge_unified.jl; not run on its own.
 ------------------------------------------------------------------------------------------
 Inputs:  v_seed.jl (warm-start patent values)
 Outputs: none
==========================================================================================#

using LinearAlgebra
using Printf

# Converged baseline V profiles, used as the default warm-start seed so that
# baseline and all counterfactuals begin near the correct equilibrium. θ is
# recovered in closed form at t=1 to match the data VA shares. See v_seed.jl.
include(joinpath(@__DIR__, "v_seed.jl"))

#------------------------------------------------------------------------------
# equilibrium_path.jl
#
# GE solver for the 3-sector (m = goods, rs = R&D-intensive services,
# ns = non-RI services) endogenous-growth, structural-change model.
# Two operating modes:
#
#   mode = :calibrate       ℓˢ(t), sectoral scientist shares and relative-price
#                           targets are DATA. The production wedges τ^p_rs, τ^p_ns
#                           are backed out as residuals here (goods = reference,
#                           τ_m ≡ 1). The scientist wedges τ^s_rs, τ^s_ns are NOT
#                           recovered inside solve_GE; they are passed in and are
#                           solved for by the SMM inner loop, jointly with ζ_j(0),
#                           to match the 1950 scientist shares.
#
#   mode = :counterfactual  the four wedges are FIXED (passed in). Sectoral
#                           scientist shares and relative prices are endogenous.
#
# In BOTH modes the aggregate research share ℓˢ(t) is exogenous:
#   L^p(t) = (1-ℓˢ)·L,   L^s(t) = ℓˢ·L.
# There is no aggregate research-share root-find; only the sectoral split of
# L^s responds to the model.
#
# Model blocks:
#   variety production  y_j = (h_j L_j^p)^α X_j^{1-α},  X_j = N_j^κ (∫ x^{μ^x})^{1/μ^x}
#   innovation          Ṅ_j = ζ_j Λ_j^γ (L_j^s)^φ,   Λ_j = ∏_i N_i^{δ_ji(t)}
#   demand              non-homothetic CES over (m, rs, ns)
#   machine markup      = 1/μ^x_j,   machine price = 1/(μ^x_j z_j)
#
# Solution: outer fixed point on patent values V_j(t) (PDV of per-variety
# profits, BGP tail horizon H = 500-T). Per period, forward-shooting:
#   Y such that Y = Y_supply(L^p, p_m(Y); μ̃^y)   [Brent]
#   p_m from the price-index normalization         [closed form]
# θ = (θ_m, θ_rs, θ_ns) is recovered in closed form at t=1 to match the data
# VA shares. δ_ji may be scalars or length-T paths.
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# 1. Brent root finder
#------------------------------------------------------------------------------

function brent_root(f::Function, a::Float64, b::Float64;
                    xtol::Float64=1e-8, ftol::Float64=1e-8, max_iter::Int=50)
    fa = f(a)
    fb = f(b)
    if fa * fb > 0.0
        return abs(fa) < abs(fb) ? (a, fa, 0) : (b, fb, 0)
    end
    if abs(fa) < abs(fb)
        a, b = b, a
        fa, fb = fb, fa
    end
    c, fc = a, fa
    mflag = true
    d = 0.0
    iter = 0
    while iter < max_iter
        iter += 1
        if abs(fb) < ftol || abs(b - a) < xtol
            return (b, fb, iter)
        end
        s = 0.0
        if fa != fc && fb != fc
            s = a*fb*fc/((fa-fb)*(fa-fc)) +
                b*fa*fc/((fb-fa)*(fb-fc)) +
                c*fa*fb/((fc-fa)*(fc-fb))
        else
            s = b - fb*(b-a)/(fb-fa)
        end
        cond1 = !(((3a + b)/4 < s < b) || (b < s < (3a + b)/4))
        cond2 = mflag  && abs(s-b) >= abs(b-c)/2
        cond3 = !mflag && abs(s-b) >= abs(c-d)/2
        cond4 = mflag  && abs(b-c) < xtol
        cond5 = !mflag && abs(c-d) < xtol
        if cond1 || cond2 || cond3 || cond4 || cond5
            s = (a + b)/2
            mflag = true
        else
            mflag = false
        end
        fs = f(s)
        d, c, fc = c, b, fb
        if fa * fs < 0
            b, fb = s, fs
        else
            a, fa = s, fs
        end
        if abs(fa) < abs(fb)
            a, b = b, a
            fa, fb = fb, fa
        end
    end
    return (b, fb, iter)
end

#------------------------------------------------------------------------------
# 2. Inline helpers — structural primitives (one-to-one with paper equations)
#------------------------------------------------------------------------------

# ----- Love-of-variety exponent (block-4 Eq. e_Y_def):
#   e_Y_j(t) = (κ + (1-μ^x_j(t)) / μ^x_j(t)) * (1-α) / α
@inline e_Y_(κ, μx, α) = (κ + (1.0 - μx)/μx) * (1.0 - α)/α

# ----- Exogenous productivity bundle (block-4 Eq. A_j_def):
#   A_j(t) = [μ^y_j · μ^x_j · (1-α) · z_j]^((1-α)/α) · h_j
@inline A_(μy, μx, z, α, hbar) = (μy * μx * (1.0 - α) * z)^((1.0 - α)/α) * hbar

# ----- Knowledge spillovers (model Eq. for Λ_j):
#   Λ_j = (N_m^Λ)^δ_jm · (N_rs^Λ)^δ_jrs · (N_ns^Λ)^δ_jns
@inline Λ_(Nm, Nrs, Nns, δm, δrs, δns) = (Nm^δm) * (Nrs^δrs) * (Nns^δns)

# ----- Relative prices from production-side free labor mobility (block-5):
#   p_k / p_m = (μ^y_m · N_m^e_Y_m · A_m / (μ^y_k · N_k^e_Y_k · A_k))^α
@inline rel_price_(μym, μyk, Nm, Nk, e_Y_m, e_Y_k, Am, Ak, α) =
    (μym * Nm^e_Y_m * Am / (μyk * Nk^e_Y_k * Ak))^α

# ----- Expenditure-share-weighted average variety markup (governs C):
#   μ̄^y = Σ_j s_j · μ^y_j
@inline μy_avg_(sm, srs, sns, μym, μyrs, μyns) = sm*μym + srs*μyrs + sns*μyns

# ----- Wedge-adjusted average markup (governs aggregate output and w^p):
#   μ̃^y = Σ_j s_j · μ^y_j · (τ^p_j)^{-1/α}   with τ^p_m ≡ 1
@inline μy_tilde_(sm, srs, sns, μym, μyrs, μyns, τp_rs, τp_ns, α) =
    sm*μym + srs*μyrs*τp_rs^(-1.0/α) + sns*μyns*τp_ns^(-1.0/α)

# ----- Expenditure-share-weighted average of μ^y·μ^x (block-6):
#   μ̄^y μ^x = Σ_j s_j · μ^y_j · μ^x_j
@inline μyμx_avg_(sm, srs, sns, μym, μyrs, μyns, μxm, μxrs, μxns) =
    sm*μym*μxm + srs*μyrs*μxrs + sns*μyns*μxns

# ----- Aggregate output anchored to sector m (block-6 Eq. aggregate_output):
#   Y(t) = (μ^y_m / μ̄^y) · p_m^(1/α) · N_m^e_Y_m · A_m · L^p
@inline Y_supply_(μym, μy_av, pm, Nm, e_Y_m, Am, Lp, α) =
    (μym/μy_av) * pm^(1.0/α) * Nm^e_Y_m * Am * Lp

# ----- Sectoral output from demand side (revenue / price):
@inline Y_demand_(s_share, Y, p) = s_share * Y / p

# ----- Conditional production labor from sector output identity:
#   L_j^p = (s_j · Y) / (p_j^(1/α) · N_j^e_Y_j · A_j)
@inline Lp_conditional_(s_share, Y, p, N, e_Y_j, A, α) =
    s_share * Y / max(1e-12, p^(1.0/α) * N^e_Y_j * A)

# ----- Per-variety machine quantity at symmetric eqm (block-6):
#   x_j(t) = μ^y_j · μ^x_j · (1-α) · z_j · p_j Y_j / N_j
@inline x_sym_(μy, μx, z, pj, Yj, Nj, α) =
    μy * μx * (1.0 - α) * z * pj * Yj / Nj

# ----- Per-variety machine-producer profit (block-5 pi_revenue):
#   π_j^x(t) = μ^y_j · (1-μ^x_j) · (1-α) · p_j Y_j / N_j
@inline π_x_(μy, μx, pj, Yj, Nj, α) =
    μy * (1.0 - μx) * (1.0 - α) * pj * Yj / Nj

# ----- Aggregate consumption (block-6 C_in_Y):
#   C(t) = [1 - (1-α) · μ̄^y μ^x] · Y
# (No equipment in lean spec.)
@inline C_(Y, μyμx_avg, α) = (1.0 - (1.0 - α) * μyμx_avg) * Y

# ----- Real interest rate (Euler, per-capita c growth; g_ℓ cancels under CRRA):
#   r(t-1) = ρ + χ · (c_t - c_{t-1}) / c_{t-1}
@inline r_euler_(ρ, χ, ct, ctm1) = ρ + χ * ((ct - ctm1) / ctm1)

# ----- Gamma function (Lanczos g=7, n=9); dependency-free. Used only for
# Γ(1-1/θ_F) in the Fréchet efficiency index, an overall constant absorbed by
# the ζ_j calibration, so full precision is not critical.
@inline function _gammafn(x::Float64)::Float64
    g = 7.0
    c = (0.99999999999980993, 676.5203681218851, -1259.1392167224028,
         771.32342877765313, -176.61502916214059, 12.507343278686905,
         -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7)
    if x < 0.5
        return π / (sin(π * x) * _gammafn(1.0 - x))
    end
    x -= 1.0
    a = c[1]
    tt = x + g + 0.5
    @inbounds for i in 2:9
        a += c[i] / (x + (i - 1))
    end
    return sqrt(2π) * tt^(x + 0.5) * exp(-tt) * a
end

#------------------------------------------------------------------------------
# 3. Innovation / scientist block (Roy / Fréchet research talent)
#
# Research talent is sector-specific: worker i draws an efficiency triplet
# {z_ij} from Fréchet F_j(z)=exp(-T_j z^{-θ}), θ>1, and freely selects the
# sector with the highest AFTER-WEDGE income max_j τ^s_j w_j z_ij. The wedge is
# applied to the allocation only, through X̃_j = τ^s_j X_j, while Ṅ_j keeps the true
# ζ_j. Equivalently τ^s_j discounts the private value of an innovation, since
# X̃_j = ζ_j Λ_j^γ (τ^s_j V_j). That asymmetry is what separates a wedge from a
# rescaling of ζ: same allocation, different idea flow, so the LP-growth moments
# tell them apart. With τ^s_j ≡ 1 the block reduces to the no-wedge case below.
#
#   Headcount shares:   π_j = T_j w_j^θ / Φ,             Φ ≡ Σ_k T_k w_k^θ
#   Selection:          max_j w_j z_ij ~ Fréchet(Φ, θ), independent of the sector
#                       chosen, so average INCOME G Φ^{1/θ} is equalized across
#                       sectors while average EFFICIENCY among those choosing j is
#                             ē_j = G Φ^{1/θ} / w_j,     G ≡ Γ(1-1/θ),
#                       which is DECREASING in w_j.
#   Efficiency units:   H_j = π_j L^s ē_j = G L^s T_j w_j^{θ-1} Φ^{1/θ-1}
#   Innovation:         Ṅ_j = ζ_j Λ_j^γ H_j^{φ}
#   Free entry:         w_j = X_j H_j^{φ-1},   X_j ≡ ζ_j Λ_j^γ V_j
#
# CORRECTION (2026-08-29). The previous implementation set H_j = ē π_j L^s with a
# COMMON ē = G Φ^{1/θ}, described as "average talent equalized". What selection
# equalizes is average income, so the previous H_j overstated efficiency units by
# a factor w_j in every sector and the innovation technology was effectively
# Ṅ_j = ζ_j Λ_j^γ (w_j H_j)^{φ}, ideas produced from the research wage bill.
#
# CLOSED FORM (requires a COMMON φ; see _roy_require_common_phi). Write
#       u = 1-φ,   D = 1 + (θ-1)u,   b = 1/D.
# Substituting H_j into free entry and solving,
#       w_j = { X_j (K T_j)^{-u} }^{b},        K ≡ G L^s Φ^{(1-θ)/θ},
#       Φ   = (G L^s)^{-θu} M^{D},             M ≡ Σ_k (T_k^{φ} X_k^{θ})^{b},
#       π_j = (T_j^{φ} X_j^{θ})^{b} / M,
#       H_j = G L^s T_j w_j^{θ-1} Φ^{1/θ-1}.
# No root find is required. The share elasticity to innovation value is
#       ε = θ b = θ / [1 + (θ-1)(1-φ)] ∈ ( 1 , 1/(1-φ) ),
# bounded BELOW by one for every θ>1 and φ∈(0,1). At θ=2, φ=0.5: ε = 4/3, against
# the ε = 1 delivered by the previous (headcount) implementation. The share
# exponent on T_j is φb = 1/3, against 1/2 previously, so holding the 1950 shares
# requires a WIDER dispersion of T_j than the superseded block reported.
#
# IDENTIFICATION WARNING. Under this block the map
#       (ζ_j, T_j) ⟼ (c_j ζ_j, T_j c_j^{-θ/φ}),  any c_j > 0,
# leaves Φ, every π_j, every w_j H_j and every Ṅ_j unchanged at every date, so the
# ζ_j levels and the T_j are jointly unidentified by any moment built from shares,
# ideas or productivity. Verified numerically in roy_invariance_check. Either
# normalize T_j ≡ 1 for all j and let ζ_j(0) match the 1950 shares, or bring in a
# moment on sectoral research wages w_j, which is the only observable that moves
# along the manifold.
#
# The aggregate research pool L^s = ℓ^s·L is fixed; the shares π_j split it. T_m ≡ 1
# (normalization); T_rs, T_ns are recovered in CLOSED FORM at t=1 to match the 1950
# scientist-share targets (see solve_period_inner!).
#------------------------------------------------------------------------------

# SUPERSEDED headcount exponent a_j = 1/(1+θ(1-φ_j)). Retained only so that stale
# call sites fail loudly rather than silently reproducing the old allocation.
@inline roy_a(θ::Float64, φ::Float64) = 1.0 / (1.0 + θ * (1.0 - φ))

# Roy exponents under a common φ.
@inline roy_D(θ::Float64, φ::Float64) = 1.0 + (θ - 1.0) * (1.0 - φ)
@inline roy_b(θ::Float64, φ::Float64) = 1.0 / roy_D(θ, φ)
# Elasticity of the headcount share to X_j.
@inline roy_eps(θ::Float64, φ::Float64) = θ * roy_b(θ, φ)
# Exponent on T_j inside the share index; equals φ·b.
@inline roy_eT(θ::Float64, φ::Float64) = φ * roy_b(θ, φ)

# The closed form above holds only for a common φ. With sector-specific φ_j the
# share carries extra factors (ē L^s)^{a_j-1} Φ^{-a_j} that no longer cancel, and
# π and Φ must be solved as a joint fixed point. Fail rather than approximate.
@inline function _roy_require_common_phi(φm::Float64, φrs::Float64, φns::Float64)::Float64
    (isapprox(φm, φrs; rtol=1e-12, atol=1e-12) &&
     isapprox(φm, φns; rtol=1e-12, atol=1e-12)) ||
        error("Roy block: the closed form requires a common φ; got " *
              "φm=$φm, φrs=$φrs, φns=$φns. Sector-specific φ_j needs the joint " *
              "fixed point on (π, Φ), which is not implemented.")
    return φm
end

# Innovation flow Ṅ_j = ζ_j Λ_j^γ H_j^{φ_j}, with H_j in efficiency units.
@inline function Ndot_sector(ζ::Float64, Λ::Float64, γ::Float64,
                             H::Float64, φ::Float64)::Float64
    if ζ <= 0.0 || Λ <= 0.0 || H <= 0.0
        return 0.0
    end
    return ζ * Λ^γ * H^φ
end

# Roy allocation headcount shares π_j = (T_j^φ X_j^θ)^b / Σ_k (T_k^φ X_k^θ)^b under a
# common φ, with b = 1/[1+(θ-1)(1-φ)]. X_j = ζ_j Λ_j^γ V_j predetermined. Inactive
# sectors (X_j ≤ 0 or T_j ≤ 0) receive zero share. Evaluated by log-sum-exp so that
# large X_j^θ does not overflow. Returns (π_m, π_rs, π_ns).
@inline function roy_shares(Tm::Float64, Trs::Float64, Tns::Float64,
                            Xm::Float64, Xrs::Float64, Xns::Float64,
                            φm::Float64, φrs::Float64, φns::Float64,
                            θ::Float64)
    φ = _roy_require_common_phi(φm, φrs, φns)
    b = roy_b(θ, φ)
    l1 = (Tm  > 0.0 && Xm  > 0.0) ? b * (φ * log(Tm)  + θ * log(Xm))  : -Inf
    l2 = (Trs > 0.0 && Xrs > 0.0) ? b * (φ * log(Trs) + θ * log(Xrs)) : -Inf
    l3 = (Tns > 0.0 && Xns > 0.0) ? b * (φ * log(Tns) + θ * log(Xns)) : -Inf
    mx = max(l1, max(l2, l3))
    isfinite(mx) || return (0.0, 0.0, 0.0)
    e1 = isfinite(l1) ? exp(l1 - mx) : 0.0
    e2 = isfinite(l2) ? exp(l2 - mx) : 0.0
    e3 = isfinite(l3) ? exp(l3 - mx) : 0.0
    S = e1 + e2 + e3
    S > 0.0 || return (0.0, 0.0, 0.0)
    return (e1 / S, e2 / S, e3 / S)
end

# Level of the Fréchet income index Φ = Σ_k T_k w_k^θ in closed form:
#   Φ = (G L^s)^{-θ(1-φ)} M^{D},  M = Σ_k (T_k^φ X_k^θ)^b,  D = 1/b.
# Returned in logs to keep M^D finite when D > 1 and M is large.
@inline function roy_logPhi(Tm::Float64, Trs::Float64, Tns::Float64,
                            Xm::Float64, Xrs::Float64, Xns::Float64,
                            φ::Float64, θ::Float64, G::Float64, Ls::Float64)::Float64
    b = roy_b(θ, φ); D = roy_D(θ, φ)
    l1 = (Tm  > 0.0 && Xm  > 0.0) ? b * (φ * log(Tm)  + θ * log(Xm))  : -Inf
    l2 = (Trs > 0.0 && Xrs > 0.0) ? b * (φ * log(Trs) + θ * log(Xrs)) : -Inf
    l3 = (Tns > 0.0 && Xns > 0.0) ? b * (φ * log(Tns) + θ * log(Xns)) : -Inf
    mx = max(l1, max(l2, l3))
    isfinite(mx) || return -Inf
    acc = (isfinite(l1) ? exp(l1 - mx) : 0.0) +
          (isfinite(l2) ? exp(l2 - mx) : 0.0) +
          (isfinite(l3) ? exp(l3 - mx) : 0.0)
    logM = mx + log(acc)
    return D * logM - θ * (1.0 - φ) * (log(G) + log(Ls))
end

# ---------------------------------------------------------------------------
# SECTOR-SPECIFIC φ_j.
#
# With a common φ the pair (π, Φ) separates and everything is closed form. With
# φ_j differing across sectors it does not, because each sector's free-entry
# condition carries its own exponent. The system still collapses to ONE scalar
# equation, however. Given Φ, Fréchet selection and free entry pin every wage
# explicitly,
#     log K   = log(G Lˢ) + ((1-θ)/θ) log Φ
#     log w_j = [ log X_j - (1-φ_j)(log K + log T_j) ] / D_j,  D_j = 1+(θ-1)(1-φ_j)
# and Φ must reproduce itself through Φ = Σ_k T_k w_k^θ. So the whole block is the
# root of
#     g(logΦ) = log Σ_k exp(log T_k + θ log w_k) − logΦ.
#
# g is strictly decreasing: dg/dlogΦ = (θ-1) Σ_k π_k (1-φ_k)/D_k − 1, and each
# term obeys (θ-1)(1-φ_k)/D_k < 1, so the weighted average is below one and the
# derivative is negative. The root is therefore unique. Bisection is used in
# preference to Newton because it cannot leave the bracket, and the block sits
# inside the SMM inner loop where a silent divergence would be costly.
#
# Setting φ_j equal in all three sectors reproduces the closed form exactly;
# test_roy_block.jl checks that against roy_solve.
# ---------------------------------------------------------------------------

@inline _roy_phi_is_common(φm::Float64, φrs::Float64, φns::Float64)::Bool =
    isapprox(φm, φrs; rtol=1e-12, atol=1e-12) &&
    isapprox(φm, φns; rtol=1e-12, atol=1e-12)

@inline function _roy_logw_j(logT::Float64, logX::Float64, φ::Float64,
                             θ::Float64, logK::Float64)::Float64
    u = 1.0 - φ
    return (logX - u * (logK + logT)) / (1.0 + (θ - 1.0) * u)
end

function roy_logPhi_het(Tm::Float64, Trs::Float64, Tns::Float64,
                        Xm::Float64, Xrs::Float64, Xns::Float64,
                        φm::Float64, φrs::Float64, φns::Float64,
                        θ::Float64, G::Float64, Ls::Float64)::Float64
    (G > 0.0 && Ls > 0.0) || return -Inf
    Ts = (Tm, Trs, Tns); Xs = (Xm, Xrs, Xns); φs = (φm, φrs, φns)
    act = (Ts[1] > 0.0 && Xs[1] > 0.0,
           Ts[2] > 0.0 && Xs[2] > 0.0,
           Ts[3] > 0.0 && Xs[3] > 0.0)
    (act[1] || act[2] || act[3]) || return -Inf
    logGLs = log(G) + log(Ls)
    @inline function g(lp::Float64)::Float64
        logK = logGLs + ((1.0 - θ) / θ) * lp
        l1 = act[1] ? log(Ts[1]) + θ*_roy_logw_j(log(Ts[1]), log(Xs[1]), φs[1], θ, logK) : -Inf
        l2 = act[2] ? log(Ts[2]) + θ*_roy_logw_j(log(Ts[2]), log(Xs[2]), φs[2], θ, logK) : -Inf
        l3 = act[3] ? log(Ts[3]) + θ*_roy_logw_j(log(Ts[3]), log(Xs[3]), φs[3], θ, logK) : -Inf
        mx = max(l1, max(l2, l3))
        isfinite(mx) || return NaN
        acc = (isfinite(l1) ? exp(l1 - mx) : 0.0) +
              (isfinite(l2) ? exp(l2 - mx) : 0.0) +
              (isfinite(l3) ? exp(l3 - mx) : 0.0)
        acc > 0.0 || return NaN
        return mx + log(acc) - lp
    end
    # Bracket around the common-φ answer at the mean φ, then widen geometrically.
    lp0 = roy_logPhi(Tm, Trs, Tns, Xm, Xrs, Xns, (φm + φrs + φns)/3.0, θ, G, Ls)
    isfinite(lp0) || (lp0 = 0.0)
    g0 = g(lp0)
    isfinite(g0) || return -Inf
    lo = lp0; hi = lp0
    if g0 > 0.0
        step = 1.0; okb = false
        for _ in 1:200
            hi = lp0 + step; gh = g(hi)
            if isfinite(gh) && gh <= 0.0; okb = true; break; end
            step *= 2.0
        end
        okb || return -Inf
    else
        step = 1.0; okb = false
        for _ in 1:200
            lo = lp0 - step; gl = g(lo)
            if isfinite(gl) && gl >= 0.0; okb = true; break; end
            step *= 2.0
        end
        okb || return -Inf
    end
    for _ in 1:300
        mid = 0.5*(lo + hi)
        gm  = g(mid)
        isfinite(gm) || return -Inf
        gm > 0.0 ? (lo = mid) : (hi = mid)
        (hi - lo) <= 1e-14 * max(1.0, abs(mid)) && break
    end
    return 0.5*(lo + hi)
end

# Full Roy allocation under sector-specific φ_j. Same return contract as roy_solve.
function _roy_solve_het(Tm::Float64, Trs::Float64, Tns::Float64,
                        Xm::Float64, Xrs::Float64, Xns::Float64,
                        φm::Float64, φrs::Float64, φns::Float64,
                        θ::Float64, G::Float64, Ls::Float64)
    logPhi = roy_logPhi_het(Tm, Trs, Tns, Xm, Xrs, Xns, φm, φrs, φns, θ, G, Ls)
    isfinite(logPhi) || return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    logK = log(G) + log(Ls) + ((1.0 - θ) / θ) * logPhi
    Ts = (Tm, Trs, Tns); Xs = (Xm, Xrs, Xns); φs = (φm, φrs, φns)
    π = zeros(3); H = zeros(3); w = zeros(3)
    @inbounds for j in 1:3
        (Ts[j] > 0.0 && Xs[j] > 0.0) || continue
        lT = log(Ts[j])
        lw = _roy_logw_j(lT, log(Xs[j]), φs[j], θ, logK)
        lπ = lT + θ*lw - logPhi
        lH = log(G) + log(Ls) + lT + (θ - 1.0)*lw + (1.0/θ - 1.0)*logPhi
        wj = exp(lw); Hj = exp(lH); πj = exp(lπ)
        w[j] = isfinite(wj) ? wj : 0.0
        H[j] = isfinite(Hj) ? Hj : 0.0
        π[j] = isfinite(πj) ? πj : 0.0
    end
    # Σ_j π_j is 1 at the root; renormalize only to absorb rounding.
    S = π[1] + π[2] + π[3]
    S > 0.0 && (π ./= S)
    return (π[1], π[2], π[3], H[1], H[2], H[3], w[1], w[2], w[3])
end

# ---------------------------------------------------------------------------
# ALTERNATIVE ALLOCATION: homogeneous researchers, free mobility, Stone-Geary
# subsistence in HEADCOUNT.
#
# Replaces the Roy/Fréchet block. Researchers are identical, so there is one
# research wage w and headcount equals efficiency units. L̄_j is committed to sector
# j the way subsistence consumption is committed to a good; only Lˢ - Σ_k L̄_k is
# allocated on the margin. Keeping the model's zero-profit convention,
#
#     w = X_j (Lˢ_j)^{φ_j-1}   =>   Lˢ_j = (X_j/w)^{1/(1-φ_j)}
#
# for any sector whose unconstrained demand exceeds its floor, and Lˢ_j = L̄_j for
# one that does not. The market clears at the w solving
#
#     Σ_j max( L̄_j , (X_j/w)^{1/(1-φ_j)} ) = Lˢ.
#
# The left side is non-increasing in w and strictly decreasing wherever a sector is
# interior, so the root is unique. Bisection is used rather than Newton because the
# max() kinks the residual at each corner and this sits inside the SMM inner loop.
#
# Setting every L̄_j = 0 gives the plain free-mobility allocation. There is no T_j
# and no θ_F here, so the (ζ,T) manifold does not arise.
# ---------------------------------------------------------------------------

@inline function _sg_demand(X::Float64, φ::Float64, logw::Float64)::Float64
    (X > 0.0) || return 0.0
    ex = 1.0 / (1.0 - φ)
    return exp(ex * (log(X) - logw))
end

function sg_solve(Xm::Float64, Xrs::Float64, Xns::Float64,
                  φm::Float64, φrs::Float64, φns::Float64,
                  Lbm::Float64, Lbrs::Float64, Lbns::Float64,
                  Ls::Float64)
    (Ls > 0.0) || return (0.0, 0.0, 0.0, 0.0)
    X = (Xm, Xrs, Xns); φ = (φm, φrs, φns); Lb = (max(Lbm,0.0), max(Lbrs,0.0), max(Lbns,0.0))
    all(f -> 0.0 < f < 1.0, φ) || error("Stone-Geary block needs φ_j in (0,1); got $φ")
    # Infeasible floors return zeros rather than throwing: inside the SMM this
    # is an ordinary bad draw and should score a penalty, not abort the run.
    sum(Lb) < Ls || return (0.0, 0.0, 0.0, 0.0)
    total(logw) = sum(max(Lb[j], _sg_demand(X[j], φ[j], logw)) for j in 1:3)
    # Bracket: total() falls as logw rises.
    lo = -40.0; hi = 40.0
    tl = total(lo); th = total(hi)
    (isfinite(tl) && isfinite(th)) || return (0.0, 0.0, 0.0, 0.0)
    tl >= Ls || return (Lb[1], Lb[2], Lb[3], exp(lo))   # even a tiny wage cannot absorb Lˢ
    th <= Ls || return (0.0, 0.0, 0.0, 0.0)
    for _ in 1:300
        mid = 0.5*(lo+hi)
        total(mid) > Ls ? (lo = mid) : (hi = mid)
        (hi - lo) <= 1e-14 * max(1.0, abs(mid)) && break
    end
    logw = 0.5*(lo+hi)
    L = ntuple(j -> max(Lb[j], _sg_demand(X[j], φ[j], logw)), 3)
    # Rounding can leave Σ L slightly off Lˢ; put the residual on the interior
    # sectors in proportion to their marginal share so the pool clears exactly.
    S = L[1]+L[2]+L[3]
    if S > 0.0 && abs(S - Ls) > 0.0
        free = ntuple(j -> max(L[j] - Lb[j], 0.0), 3)
        fs = free[1]+free[2]+free[3]
        L = fs > 0.0 ? ntuple(j -> Lb[j] + free[j]*(Ls - Lb[1]-Lb[2]-Lb[3])/fs, 3) :
                       ntuple(j -> Lb[j], 3)
    end
    return (L[1], L[2], L[3], exp(logw))
end

# Full Roy allocation under a common φ. Returns
# (π_m,π_rs,π_ns, H_m,H_rs,H_ns, w_m,w_rs,w_ns) with H_j in EFFICIENCY UNITS,
#   H_j = π_j L^s ē_j,  ē_j = G Φ^{1/θ}/w_j,  w_j = X_j H_j^{φ-1}.
# The aggregate research wage bill is Σ_j w_j H_j = G Φ^{1/θ} L^s, so the average
# research wage per worker is G Φ^{1/θ}, recoverable ex post from the stored arrays
# as (Σ_j w_j H_j)/L^s. It is not any w_j and not a share-weighted average of them.
function roy_solve(Tm::Float64, Trs::Float64, Tns::Float64,
                   Xm::Float64, Xrs::Float64, Xns::Float64,
                   φm::Float64, φrs::Float64, φns::Float64,
                   θ::Float64, G::Float64, Ls::Float64)
    (Ls > 0.0 && G > 0.0) ||
        return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    # Sector-specific φ_j has no closed form; use the scalar fixed point.
    _roy_phi_is_common(φm, φrs, φns) ||
        return _roy_solve_het(Tm, Trs, Tns, Xm, Xrs, Xns, φm, φrs, φns, θ, G, Ls)
    φ = φm
    πm, πrs, πns = roy_shares(Tm, Trs, Tns, Xm, Xrs, Xns, φm, φrs, φns, θ)
    (πm + πrs + πns) > 0.0 ||
        return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    logPhi = roy_logPhi(Tm, Trs, Tns, Xm, Xrs, Xns, φ, θ, G, Ls)
    isfinite(logPhi) ||
        return (πm, πrs, πns, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    b = roy_b(θ, φ)
    # log K = log(G L^s) + ((1-θ)/θ) log Φ ;  log w_j = b[ log X_j - (1-φ)(log K + log T_j) ]
    logK = log(G) + log(Ls) + ((1.0 - θ) / θ) * logPhi
    @inline function _wH(T::Float64, X::Float64)
        (T > 0.0 && X > 0.0) || return (0.0, 0.0)
        lw = b * (log(X) - (1.0 - φ) * (logK + log(T)))
        # H_j = G L^s T_j w_j^{θ-1} Φ^{1/θ-1}
        lH = log(G) + log(Ls) + log(T) + (θ - 1.0) * lw + (1.0 / θ - 1.0) * logPhi
        w = exp(lw); H = exp(lH)
        return (isfinite(w) ? w : 0.0, isfinite(H) ? H : 0.0)
    end
    w_m,  H_m  = _wH(Tm,  Xm)
    w_rs, H_rs = _wH(Trs, Xrs)
    w_ns, H_ns = _wH(Tns, Xns)
    return (πm, πrs, πns, H_m, H_rs, H_ns, w_m, w_rs, w_ns)
end

# Closed-form recovery of T_j (with T_m ≡ 1) from date-1 headcount share targets:
#   T_j = [ (π_j/π_m)^{D} (X_m/X_j)^{θ} ]^{1/φ},   D = 1 + (θ-1)(1-φ).
# The superseded block used the exponents (1/a_j, θ) with a_j = 1/(1+θ(1-φ_j)).
@inline function roy_recover_T(π_target::Float64, π_m::Float64,
                               X_j::Float64, X_m::Float64,
                               φ::Float64, θ::Float64)::Float64
    (π_target > 0.0 && π_m > 0.0 && X_j > 0.0 && X_m > 0.0) || return 1.0
    D = roy_D(θ, φ)
    lT = (D * (log(π_target) - log(π_m)) + θ * (log(X_m) - log(X_j))) / φ
    lT = clamp(lT, -300.0, 300.0)
    T = exp(lT)
    return (isfinite(T) && T > 0.0) ? T : 1.0
end

# Diagnostic: the (ζ_j, T_j) ⟼ (c_j ζ_j, T_j c_j^{-θ/φ}) manifold leaves the shares,
# Φ, the wage bills and the idea flows unchanged. Returns the largest relative
# deviation across those objects; a value at solver precision confirms that the ζ
# levels are not identified by any moment built from them.
function roy_invariance_check(Tm, Trs, Tns, Xm, Xrs, Xns, φ, θ, G, Ls,
                              c = (1.0, 0.7, 1.4))
    base = roy_solve(Tm, Trs, Tns, Xm, Xrs, Xns, φ, φ, φ, θ, G, Ls)
    e = -θ / φ
    alt = roy_solve(Tm * c[1]^e, Trs * c[2]^e, Tns * c[3]^e,
                    Xm * c[1], Xrs * c[2], Xns * c[3], φ, φ, φ, θ, G, Ls)
    dev = 0.0
    for j in 1:3                                    # π_j
        dev = max(dev, abs(alt[j] - base[j]))
    end
    for (j, cj) in enumerate(c)                     # Ṅ_j ∝ (X_j/V_j) H_j^φ
        b0 = base[3+j] > 0.0 ? base[3+j]^φ : 0.0
        b1 = alt[3+j]  > 0.0 ? alt[3+j]^φ  : 0.0
        b0 > 0.0 && (dev = max(dev, abs(cj * b1 / b0 - 1.0)))
    end
    return dev
end

#------------------------------------------------------------------------------
# 4. GEState: struct-of-arrays for all length-T equilibrium variables.
#
# The naming uses subscript "t" for time-series (so Yt is the path of Y).
# Sector labels: m = goods/MFG, rs = R&D-intensive services, ns = non-RI services.
#------------------------------------------------------------------------------

mutable struct GEState
    # Aggregate output, sectoral outputs
    Yt::Vector{Float64};   Ymt::Vector{Float64};   Yrst::Vector{Float64};   Ynst::Vector{Float64}

    # Sectoral prices (numeraire pinned at t=1)
    pmt::Vector{Float64};  prst::Vector{Float64};  pnst::Vector{Float64}

    # Sectoral expenditure shares of total Y
    smt::Vector{Float64};  srst::Vector{Float64};  snst::Vector{Float64}

    # Per-variety machine quantities and profits
    xmt::Vector{Float64};  xrst::Vector{Float64};  xnst::Vector{Float64}
    πxmt::Vector{Float64}; πxrst::Vector{Float64}; πxnst::Vector{Float64}

    # Exogenous productivity bundles A_j(t)
    Amt::Vector{Float64};  Arst::Vector{Float64};  Anst::Vector{Float64}

    # Love-of-variety exponents e_Y_j(t)
    e_Y_mt::Vector{Float64}; e_Y_rst::Vector{Float64}; e_Y_nst::Vector{Float64}

    # Consumption, interest rate
    Ct::Vector{Float64}
    rt::Vector{Float64}

    # Patent values
    Vmt::Vector{Float64};  Vrst::Vector{Float64};  Vnst::Vector{Float64}

    # Aggregate markup wedges
    μy_avt::Vector{Float64};   μyμx_avt::Vector{Float64}

    # Sectoral production labor
    Lpmt::Vector{Float64};  Lprst::Vector{Float64}; Lpnst::Vector{Float64}

    # Labor totals
    Lt::Vector{Float64}; Lst::Vector{Float64}; Lpt::Vector{Float64}

    # Knowledge spillovers Λ_j(t), Lambda-only stocks N_j^Λ(t), and normalized varieties N_j(t)
    Λmt::Vector{Float64};  Λrst::Vector{Float64};  Λnst::Vector{Float64}
    NΛmt::Vector{Float64}; NΛrst::Vector{Float64}; NΛnst::Vector{Float64}
    Nmt::Vector{Float64};  Nrst::Vector{Float64};  Nnst::Vector{Float64}

    # Sectoral scientist allocation
    Lsmt::Vector{Float64}; Lsrst::Vector{Float64}; Lsnst::Vector{Float64}

    # Sectoral researcher efficiency units used by the innovation technology
    Hmt::Vector{Float64}; Hrst::Vector{Float64}; Hnst::Vector{Float64}

    # Variety creation flows Ṅ_j
    Nmdot::Vector{Float64}; Nrsdot::Vector{Float64}; Nnsdot::Vector{Float64}

    # Production wages w_j^p (numeraire-relative; goods = w_m^p)
    wpmt::Vector{Float64}; wprst::Vector{Float64}; wpnst::Vector{Float64}

    # Scientist wages w_j^s
    wsmt::Vector{Float64}; wsrst::Vector{Float64}; wsnst::Vector{Float64}
end

function GEState(T::Int)
    z = () -> zeros(T)
    return GEState(
        z(), z(), z(), z(),    # Yt..Ynst
        z(), z(), z(),         # pmt..pnst
        z(), z(), z(),         # smt..snst
        z(), z(), z(),         # xmt..xnst
        z(), z(), z(),         # πxmt..πxnst
        z(), z(), z(),         # Amt..Anst
        z(), z(), z(),         # e_Y_mt..e_Y_nst
        z(),                   # Ct
        z(),                   # rt
        z(), z(), z(),         # Vmt..Vnst
        z(), z(),              # μy_avt, μyμx_avt
        z(), z(), z(),         # Lpmt..Lpnst
        z(), z(), z(),         # Lt, Lst, Lpt
        z(), z(), z(),         # Λmt..Λnst
        z(), z(), z(),         # NΛmt..NΛnst
        z(), z(), z(),         # Nmt..Nnst
        z(), z(), z(),         # Lsmt..Lsnst
        z(), z(), z(),         # Hmt..Hnst
        z(), z(), z(),         # Nmdot..Nnsdot
        z(), z(), z(),         # wpmt..wpnst
        z(), z(), z(),         # wsmt..wsnst
    )
end

# Pack GEState into a Dict matching v8's expected output keys (so callers
# downstream — smm_ge_unified.jl, plotting scripts — don't break).
function state_to_dict(s::GEState)::Dict{String, Vector{Float64}}
    d = Dict{String, Vector{Float64}}()
    d["Yt"]    = s.Yt;     d["Ymt"]   = s.Ymt;    d["Yrst"]  = s.Yrst;   d["Ynst"]  = s.Ynst
    d["pmt"]   = s.pmt;    d["prst"]  = s.prst;   d["pnst"]  = s.pnst
    d["smt"]   = s.smt;    d["srst"]  = s.srst;   d["snst"]  = s.snst
    d["xmt"]   = s.xmt;    d["xrst"]  = s.xrst;   d["xnst"]  = s.xnst
    d["πxmt"]  = s.πxmt;   d["πxrst"] = s.πxrst;  d["πxnst"] = s.πxnst
    d["Amt"]   = s.Amt;    d["Arst"]  = s.Arst;   d["Anst"]  = s.Anst
    d["e_Y_mt"]  = s.e_Y_mt;  d["e_Y_rst"] = s.e_Y_rst;  d["e_Y_nst"] = s.e_Y_nst
    d["Ct"]    = s.Ct
    d["rt"]    = s.rt
    d["Vmt"]   = s.Vmt;    d["Vrst"]  = s.Vrst;   d["Vnst"]  = s.Vnst
    d["μy_av"]   = s.μy_avt;   d["μyμx_av"] = s.μyμx_avt
    d["Lpmt"]  = s.Lpmt;   d["Lprst"] = s.Lprst;  d["Lpnst"] = s.Lpnst
    d["Lt"]    = s.Lt;     d["Lst"]   = s.Lst;    d["Lpt"]   = s.Lpt
    d["Λmt"]   = s.Λmt;    d["Λrst"]  = s.Λrst;   d["Λnst"]  = s.Λnst
    d["NΛmt"]  = s.NΛmt;   d["NΛrst"] = s.NΛrst;  d["NΛnst"] = s.NΛnst
    d["Nmt"]   = s.Nmt;    d["Nrst"]  = s.Nrst;   d["Nnst"]  = s.Nnst
    d["Lsmt"]  = s.Lsmt;   d["Lsrst"] = s.Lsrst;  d["Lsnst"] = s.Lsnst
    d["Hmt"]   = s.Hmt;    d["Hrst"]  = s.Hrst;   d["Hnst"]  = s.Hnst
    d["Nmdot"] = s.Nmdot;  d["Nrsdot"]= s.Nrsdot; d["Nnsdot"]= s.Nnsdot
    d["wpmt"]  = s.wpmt;   d["wprst"] = s.wprst;  d["wpnst"] = s.wpnst
    d["wsmt"]  = s.wsmt;   d["wsrst"] = s.wsrst;  d["wsnst"] = s.wsnst
    d["wpt"]   = s.wpmt    # back-compat alias for the goods production wage
    return d
end

#------------------------------------------------------------------------------
# 5. solve_GE — main entry point
#
# Positional inputs:
#   T                    : number of periods
#   χ, ρ                 : CRRA curvature, discount rate
#   δ_*                  : 9 knowledge-spillover intensities; each a scalar or a
#                          length-T path. Every receiving-sector row must sum to 1.
#   α, κ                 : labor share and love-of-variety scaling
#   σ, η_*               : non-homothetic CES preference parameters
#   γ, φ                 : innovation exponents (spillover, scientist labor)
#   θm_in, θrs_in        : θ seeds used when equilibrium θ is recovered at t=1
#   gζ_*, ζ*0            : innovation-productivity growth rates and t=1 levels
#   ls_path              : aggregate research share ℓ^s(t), length T (EXOGENOUS)
#   sci_sh_rs/ns_target  : t=1 sectoral scientist-share targets used to recover
#                          the two service-sector Fréchet scales when not fixed
#   p_ratio_rs/ns_target : relative-price targets p_rs/p_m, p_ns/p_m (calibrate);
#                          same length T_data ≤ T and same wedge-freeze rule.
#   smt_1/srst_1_target  : t=1 VA-share targets that pin θ
#   ini_N_*, gl          : initial varieties (normalized to 1) and population path
#   μx*, μy*, z*, hbar*  : exogenous length-T paths per sector
#   Lt                   : total employment path (length T)
#
# Keyword inputs:
#   mode                 : :calibrate (default) or :counterfactual
#   V_init_*             : warm-start patent values (length T)
#   v_damp, iter_max_v, tol_v: V fixed-point damping, iteration cap, and tolerance
#   θfr                  : Fréchet shape parameter (θ_F > 1); default 2.0
#   T_m                  : Fréchet goods scale (normalization); default 1.0
#   T_rs_fixed, T_ns_fixed: optional fixed Fréchet service-sector scales. When
#                          omitted, they are recovered from t=1 scientist shares.
#   theta_m/rs_fixed     : optional fixed CES weights. When omitted, θ is
#                          recovered from t=1 VA shares; both must be supplied.
#   φrs, φns             : scientist-labor elasticities (default = φm; Roy
#                          efficiency index assumes sector-neutral φ)
#   tau_p_rs/ns_path     : the two fixed production wedges; REQUIRED when
#                          mode=:counterfactual, ignored otherwise
#   lambda_N0            : 1950 patent shares (m, rs, ns), used exclusively to
#                          seed N_j^Λ inside the knowledge aggregator
#
# The Fréchet scales and CES weights are recovered independently in closed form
# at t=1 unless their paired fixed-value keywords are supplied.
#
# Returns Dict{String,Vector{Float64}} of all eqm paths (incl. w_j^p, w_j^s,
# π_j scientist shares, τ^p_j, recovered T_rs/T_ns) plus solver diagnostics.
#------------------------------------------------------------------------------

function solve_GE(
    T, χ, ρ,
    δmm, δmr, δmn, δnm, δnr, δnn, δrm, δrr, δrn,
    α, κ,
    σ, ηm, ηrs, ηns, γ, φm, θm_in, θrs_in,
    gζm, gζrs, gζns,
    ζm0, ζrs0, ζns0,
    ls_path,
    sci_sh_rs_target, sci_sh_ns_target,
    p_ratio_rs_target, p_ratio_ns_target,
    smt_1_target, srst_1_target,
    ini_Nm, ini_Nrs, ini_Nns, gl,
    μxmt, μxrst, μxnst, μymt, μyrst, μynst,
    zmt, zrst, znst, hbarmt, hbarrst, hbarnst, Lt;
    mode::Symbol=:calibrate,
    V_init_m::Union{Nothing,Vector{Float64}}=nothing,
    V_init_rs::Union{Nothing,Vector{Float64}}=nothing,
    V_init_ns::Union{Nothing,Vector{Float64}}=nothing,
    v_damp::Float64=0.9,
    iter_max_v::Int=45,
    tol_v::Float64=1e-3,
    θfr::Float64=2.0,
    T_m::Float64=1.0,
    T_rs_fixed::Union{Nothing,Float64}=nothing,
    T_ns_fixed::Union{Nothing,Float64}=nothing,
    # Scientist mobility wedges on the Roy selection margin, goods = reference
    # (τ^s_m ≡ 1). Researchers select on after-wedge income τ^s_j w_j z_ij, so the
    # allocation sees X̃_j = τ^s_j X_j while the idea flow keeps the true ζ_j. That
    # is what separates a wedge from a rescaling of ζ: same allocation, different
    # ideas, so the LP-growth moments tell them apart.
    tau_s_rs::Float64=1.0,
    tau_s_ns::Float64=1.0,
    # Allocation of researchers. :roy is the Fréchet/Roy block; :stonegeary is the
    # homogeneous free-mobility block with subsistence headcounts, which ignores
    # T_j and θ_F entirely. The L̄_j are CONSTANT IN LEVELS, so they bind early and
    # fade as Lˢ grows, which is the subsistence-consumption analogue.
    alloc::Symbol=:roy,
    Lbar_m::Float64=0.0,
    Lbar_rs::Float64=0.0,
    Lbar_ns::Float64=0.0,
    theta_m_fixed::Union{Nothing,Float64}=nothing,
    theta_rs_fixed::Union{Nothing,Float64}=nothing,
    theta_m_path::Union{Nothing,Vector{Float64}}=nothing,
    theta_rs_path::Union{Nothing,Vector{Float64}}=nothing,
    N_m_lock::Union{Nothing,Vector{Float64}}=nothing,
    N_rs_lock::Union{Nothing,Vector{Float64}}=nothing,
    N_ns_lock::Union{Nothing,Vector{Float64}}=nothing,
    lambda_N0::NamedTuple,
    tau_p_rs_path::Union{Nothing,Vector{Float64}}=nothing,
    tau_p_ns_path::Union{Nothing,Vector{Float64}}=nothing,
    φrs::Float64=φm,
    φns::Float64=φm,
    # Full research-productivity paths (m, rs, ns), length T; replaces ζ0 exp(g t) when given.
    ζ_paths::Union{Nothing,NTuple{3,Vector{Float64}}}=nothing,
)
    s = GEState(T)

    # ----- Iteration / tolerance parameters -----
    iter_max_y   = 100
    # tol_y feeds p_m (raised to 1/α), so keep it far tighter than Y itself.
    tol_y        = 1e-10
    isfinite(tol_v) && tol_v > 0.0 || error("tol_v must be finite and strictly positive")
    ω_v          = v_damp   # damping on V update
    FLOW_CAP     = 1e8
    N_CAP        = 1e10

    # ----- Diagnostics -----
    nfev_pm       = Ref(0)
    nfev_y        = Ref(0)
    iter_pm_total = Ref(0)
    iter_y_total  = Ref(0)
    n_v_iters     = Ref(0)
    T_rs_const    = Ref(something(T_rs_fixed, 1.0))
    T_ns_const    = Ref(something(T_ns_fixed, 1.0))

    all(isfinite, (T_m, T_rs_const[], T_ns_const[])) ||
        error("Fréchet scales must be finite")
    all(>(0.0), (T_m, T_rs_const[], T_ns_const[])) ||
        error("Fréchet scales must be strictly positive")

    fixed_theta = theta_m_fixed !== nothing && theta_rs_fixed !== nothing
    (theta_m_fixed === nothing) == (theta_rs_fixed === nothing) ||
        error("theta_m_fixed and theta_rs_fixed must be supplied together")
    if fixed_theta
        theta_ns_fixed = 1.0 - theta_m_fixed - theta_rs_fixed
        all(isfinite, (theta_m_fixed, theta_rs_fixed, theta_ns_fixed)) ||
            error("Fixed preference weights must be finite")
        all(>(0.0), (theta_m_fixed, theta_rs_fixed, theta_ns_fixed)) ||
            error("Fixed preference weights must be strictly positive and sum to one")
    end
    path_theta = theta_m_path !== nothing && theta_rs_path !== nothing
    (theta_m_path === nothing) == (theta_rs_path === nothing) ||
        error("theta_m_path and theta_rs_path must be supplied together")
    if path_theta
        length(theta_m_path) == T && length(theta_rs_path) == T ||
            error("time-varying theta paths must have length T")
        theta_ns_path = 1 .- theta_m_path .- theta_rs_path
        all(isfinite, theta_m_path) && all(isfinite, theta_rs_path) && all(isfinite, theta_ns_path) ||
            error("time-varying theta paths must be finite")
        all(>(0.0), theta_m_path) && all(>(0.0), theta_rs_path) && all(>(0.0), theta_ns_path) ||
            error("time-varying theta weights must be strictly positive and sum to one")
    end
    lock_N = N_m_lock !== nothing && N_rs_lock !== nothing && N_ns_lock !== nothing
    all(x -> x === nothing, (N_m_lock, N_rs_lock, N_ns_lock)) || lock_N ||
        error("all three knowledge-state lock paths must be supplied together")
    if lock_N
        all(length(x) == T for x in (N_m_lock, N_rs_lock, N_ns_lock)) ||
            error("knowledge-state lock paths must have length T")
    end

    (:m in keys(lambda_N0) && :rs in keys(lambda_N0) && :ns in keys(lambda_N0)) ||
        error("lambda_N0 must contain fields :m, :rs, and :ns")
    qλ_m, qλ_rs, qλ_ns = Float64(lambda_N0.m), Float64(lambda_N0.rs), Float64(lambda_N0.ns)
    all(isfinite, (qλ_m, qλ_rs, qλ_ns)) || error("lambda_N0 shares must be finite")
    all(>(0.0), (qλ_m, qλ_rs, qλ_ns)) || error("lambda_N0 shares must be strictly positive")
    isapprox(qλ_m + qλ_rs + qλ_ns, 1.0; atol=1e-12, rtol=1e-12) ||
        error("lambda_N0 shares must sum to one")
    # Two admissible initialisations. Either N_j(1950)=1, the normalised index, in
    # which case the 1950 cross-section enters ONLY through lambda_N0 in the
    # spillover aggregator; or N_j(1950)=q_j, the observed patent shares, in which
    # case N carries the level everywhere and the spillover must use N directly or
    # the shares would be applied twice.
    _N_LEVELS = !all(x -> isapprox(x, 1.0; atol=1e-12, rtol=1e-12),
                           (ini_Nm, ini_Nrs, ini_Nns))
    if _N_LEVELS
        all(isapprox(a, b; rtol=1e-9) for (a,b) in
            ((ini_Nm, qλ_m), (ini_Nrs, qλ_rs), (ini_Nns, qλ_ns))) ||
            error("N_j(1950) must equal either 1 or the lambda_N0 shares")
    end

    # N_j^Λ starts from the data patent share and inherits the endogenous gross
    # growth of normalized N_j. Since N_j(1950)=1, N_j^Λ(t)=q_j(1950)N_j(t).
    @inline function refresh_lambda!(t::Int)
        # When N_j(1950)=q_j the level is already in N, so scaling by q_j again
        # would double-count the 1950 cross-section.
        s.NΛmt[t]  = _N_LEVELS ? s.Nmt[t] : qλ_m  * s.Nmt[t]
        s.NΛrst[t] = _N_LEVELS ? s.Nrst[t] : qλ_rs * s.Nrst[t]
        s.NΛnst[t] = _N_LEVELS ? s.Nnst[t] : qλ_ns * s.Nnst[t]
        s.Λmt[t]  = Λ_(s.NΛmt[t], s.NΛrst[t], s.NΛnst[t], _δt(δmm,t), _δt(δmr,t), _δt(δmn,t))
        s.Λrst[t] = Λ_(s.NΛmt[t], s.NΛrst[t], s.NΛnst[t], _δt(δrm,t), _δt(δrr,t), _δt(δrn,t))
        s.Λnst[t] = Λ_(s.NΛmt[t], s.NΛrst[t], s.NΛnst[t], _δt(δnm,t), _δt(δnr,t), _δt(δnn,t))
        return nothing
    end

    # ----- Mode / wedge configuration -----
    counterfactual = (mode === :counterfactual)
    if counterfactual
        all(p -> p !== nothing && length(p) == T,
            (tau_p_rs_path, tau_p_ns_path)) ||
            error("mode=:counterfactual requires τ^p paths of length T")
    elseif mode !== :calibrate
        error("mode must be :calibrate or :counterfactual (got $mode)")
    end

    κ isa NamedTuple || error("κ must be a NamedTuple with fields (:m, :rs, :ns); got $(typeof(κ))")
    (:m in keys(κ) && :rs in keys(κ) && :ns in keys(κ)) ||
        error("κ NamedTuple must contain fields :m, :rs, :ns")
    κ_m = Float64(κ.m)
    κ_rs = Float64(κ.rs)
    κ_ns = Float64(κ.ns)
    all(isfinite, (κ_m, κ_rs, κ_ns)) || error("κ values must be finite")
    all(>(0.0), (κ_m, κ_rs, κ_ns)) || error("κ values must be strictly positive")

    # ----- Calibrate-mode wedge horizon -----
    # In :calibrate the production-price targets (p_ratio_rs/ns) may be shorter
    # than T. The production wedges τ^p_j are recovered to match those targets
    # for t ≤ T_data and held fixed at their t=T_data (2023) values thereafter.
    # The Fréchet scales T_rs, T_ns use only the t=1 scientist-share targets.
    T_data = T
    if !counterfactual
        T_data = length(p_ratio_rs_target)
        length(p_ratio_ns_target) == T_data ||
            error(":calibrate price targets must share one length T_data (got p_rs=$(length(p_ratio_rs_target)), " *
                  "p_ns=$(length(p_ratio_ns_target)))")
        1 <= T_data <= T || error("calibrate target length T_data=$T_data must satisfy 1 ≤ T_data ≤ T=$T")
        (length(sci_sh_rs_target) >= 1 && length(sci_sh_ns_target) >= 1) ||
            error(":calibrate requires scientist-share targets of length ≥ 1 (t=1 used for T_rs, T_ns)")
    end

    # δ_ji may be scalars or length-T paths; index safely and check row sums.
    _δt(v, t) = v isa AbstractVector ? @inbounds(v[t]) : v
    for t in 1:T
        @assert abs(_δt(δmm,t)+_δt(δmr,t)+_δt(δmn,t) - 1.0) < 1e-6 "δ_m row ≠ 1 at t=$t"
        @assert abs(_δt(δrm,t)+_δt(δrr,t)+_δt(δrn,t) - 1.0) < 1e-6 "δ_rs row ≠ 1 at t=$t"
        @assert abs(_δt(δnm,t)+_δt(δnr,t)+_δt(δnn,t) - 1.0) < 1e-6 "δ_ns row ≠ 1 at t=$t"
    end

    # ----- Exogenous ζ_j paths; price wedges recorded for reporting -----
    ζmt  = [ζm0  * exp(gζm  * (t - 1)) for t in 1:T]
    ζrst = [ζrs0 * exp(gζrs * (t - 1)) for t in 1:T]
    ζnst = [ζns0 * exp(gζns * (t - 1)) for t in 1:T]
    if ζ_paths !== nothing
        all(length(v) == T for v in ζ_paths) || error("ζ_paths must have length T=$T")
        ζmt, ζrst, ζnst = ζ_paths
    end
    τprst = zeros(T)     # τ^p_rs (report)
    τpnst = zeros(T)     # τ^p_ns (report)

    # ----- Targets used to recover theta internally -----
    # theta (θ_m, θ_rs, θ_ns) is recovered in closed form at t=1 so that the
    # model's t=1 expenditure shares equal the data VA-share targets
    # (smt_1_target, srst_1_target); θ_ns = 1 - θ_m - θ_rs. See the θ-recovery
    # block after the t=1 solve inside the V fixed-point loop.

    # ----- Initial varieties: taken as passed (either 1, or the 1950 shares) -----

    # ----- Precompute time-varying exponents and productivity scales -----
    # These depend only on time-varying μ^x_j(t) and exogenous {μ^y, z}, not on
    # equilibrium objects, so we compute once at solve_GE entry.
    @inbounds for t in 1:T
        s.e_Y_mt[t]  = e_Y_(κ_m,  μxmt[t],  α)
        s.e_Y_rst[t] = e_Y_(κ_rs, μxrst[t], α)
        s.e_Y_nst[t] = e_Y_(κ_ns, μxnst[t], α)
    end

    # ----------------------------------------------------------------------
    # solve_period_inner!:
    #   Given V_j(t), Λ_j(t), N_j(t) and θ, solve the within-period equilibrium.
    #   Lᵖ = (1-ℓˢ)L is known, so only Y is solved:
    #     Y such that Y = Y_supply(Lᵖ, p_m(Y); μ̃^y)   [Brent, p_m closed form].
    #   Then split the fixed pool Lˢ = ℓˢL across sectors and compute N flows and
    #   all sectoral wages. Wedges τ^p, τ^s are data (:counterfactual) or backed
    #   out from the price / scientist-share targets (:calibrate).
    #
    #   The closed-form p_m exploits p_j = (p_j/p_m)·p_m, so the price index
    #     Σ_j θ_j p_j^(1-σ) Y^((η_j-1)(1-σ)) = 1
    #   is linear in p_m^(1-σ).
    # ----------------------------------------------------------------------
    function solve_period_inner!(
        t::Int, V_m::Float64, V_rs::Float64, V_ns::Float64,
        θm::Float64, θrs::Float64, θns::Float64,
        Y_init::Float64;
        match_t1_shares::Bool=false
    )
        # Productivity bundles at time t
        s.Amt[t]  = A_(μymt[t],  μxmt[t],  zmt[t],  α, hbarmt[t])
        s.Arst[t] = A_(μyrst[t], μxrst[t], zrst[t], α, hbarrst[t])
        s.Anst[t] = A_(μynst[t], μxnst[t], znst[t], α, hbarnst[t])

        Y_ref  = Y_init
        e_Y_m  = s.e_Y_mt[t]
        e_Y_rs = s.e_Y_rst[t]
        e_Y_ns = s.e_Y_nst[t]

        # Fixed research share splits total labor (no L^p root-find).
        ls_t   = clamp(ls_path[t], 0.0, 1.0 - 1e-9)
        Lp_t   = (1.0 - ls_t) * s.Lt[t]
        Ls_tot = ls_t * s.Lt[t]
        s.Lpt[t] = Lp_t

        # Fundamental relative prices p_j/p_m (before mobility wedge).
        rp_rs = rel_price_(μymt[t], μyrst[t], s.Nmt[t], s.Nrst[t],
                           e_Y_m, e_Y_rs, s.Amt[t], s.Arst[t], α)
        rp_ns = rel_price_(μymt[t], μynst[t], s.Nmt[t], s.Nnst[t],
                           e_Y_m, e_Y_ns, s.Amt[t], s.Anst[t], α)

        # Production-mobility wedge τ^p_j: data in :counterfactual, else the
        # residual that reproduces the target relative price.
        if counterfactual
            τp_rs_t = max(tau_p_rs_path[t], 1e-8)
            τp_ns_t = max(tau_p_ns_path[t], 1e-8)
        elseif t <= T_data
            τp_rs_t = max(p_ratio_rs_target[t], 1e-8) / max(rp_rs, 1e-12)
            τp_ns_t = max(p_ratio_ns_target[t], 1e-8) / max(rp_ns, 1e-12)
        else
            # t > T_data: hold production wedges at their t=T_data (2023) values.
            τp_rs_t = τprst[T_data]
            τp_ns_t = τpnst[T_data]
        end
        p_rs_ratio = rp_rs * τp_rs_t      # observed relative price p_rs/p_m
        p_ns_ratio = rp_ns * τp_ns_t
        τprst[t] = τp_rs_t
        τpnst[t] = τp_ns_t

        ε       = 1.0 - σ
        ln_κ_rs = log(p_rs_ratio)
        ln_κ_ns = log(p_ns_ratio)          # ln(p_m/p_m) = 0

        # t=1 VA-share targets, used only when match_t1_shares=true to recover θ.
        sm_star  = clamp(smt_1_target,  1e-8, 1.0 - 2e-8)
        srs_star = clamp(srst_1_target, 1e-8, 1.0 - 2e-8)
        sns_star = max(1.0 - sm_star - srs_star, 1e-8)

        # ---- p_m closed form given Y (also sets the softmax shares) ----
        function solve_pm_given_Y!(Y_local::Float64)
            lnY  = log(Y_local)
            a_m  =           (ηm  - 1.0) * lnY
            a_rs = ln_κ_rs + (ηrs - 1.0) * lnY
            a_ns = ln_κ_ns + (ηns - 1.0) * lnY

            local pm::Float64
            if match_t1_shares
                # Impose shares s_j* and solve p_m from Σ_j θ_j = 1; θ recovered after.
                if ε == 0.0
                    pm = exp(-(sm_star*a_m + srs_star*a_rs + sns_star*a_ns))
                else
                    Z  = sm_star*exp(-ε*a_m) + srs_star*exp(-ε*a_rs) + sns_star*exp(-ε*a_ns)
                    pm = exp(log(Z) / ε)
                end
                s.smt[t]  = sm_star
                s.srst[t] = srs_star
                s.snst[t] = sns_star
            elseif ε == 0.0
                pm        = exp(-(θm*a_m + θrs*a_rs + θns*a_ns))
                s.smt[t]  = θm
                s.srst[t] = θrs
                s.snst[t] = θns
            else
                S   = θm*expm1(ε*a_m) + θrs*expm1(ε*a_rs) + θns*expm1(ε*a_ns)
                pm  = exp(-log1p(S) / ε)
                den = 1.0 + S
                s.smt[t]  = θm  * (1.0 + expm1(ε*a_m))  / den
                s.srst[t] = θrs * (1.0 + expm1(ε*a_rs)) / den
                s.snst[t] = θns * (1.0 + expm1(ε*a_ns)) / den
            end

            nfev_pm[]       += 1
            iter_pm_total[] += 1
            s.pmt[t]  = pm
            s.prst[t] = p_rs_ratio * pm
            s.pnst[t] = p_ns_ratio * pm
            return 0.0
        end

        # ---- Y such that Y = Y_supply(Lᵖ, p_m(Y); μ̃^y) ----
        f_Y = function (Y_try)
            solve_pm_given_Y!(Y_try)
            μ̃y = μy_tilde_(s.smt[t], s.srst[t], s.snst[t],
                           μymt[t], μyrst[t], μynst[t], τp_rs_t, τp_ns_t, α)
            s.μy_avt[t] = μ̃y
            Y_rhs = Y_supply_(μymt[t], μ̃y, s.pmt[t],
                              s.Nmt[t], e_Y_m, s.Amt[t], Lp_t, α)
            return Y_try - Y_rhs
        end

        # Log-symmetric bracket around the warm start, widened by doubling.
        lnY0 = log(max(Y_ref, 1e-8))
        w    = log(2.0)
        Y_lo = exp(lnY0 - w); Y_hi = exp(lnY0 + w)
        f_lo = f_Y(Y_lo);     f_hi = f_Y(Y_hi)
        expand_iter = 0
        while f_lo * f_hi > 0.0 && expand_iter < 12
            w   *= 2.0
            Y_lo = exp(lnY0 - w); Y_hi = exp(lnY0 + w)
            f_lo = f_Y(Y_lo);     f_hi = f_Y(Y_hi)
            expand_iter += 1
        end

        Y_root, fY, iters = brent_root(f_Y, Y_lo, Y_hi;
                                       xtol=tol_y, ftol=tol_y, max_iter=iter_max_y)
        nfev_y[]       += 1
        iter_y_total[] += iters
        f_Y(Y_root)            # restore state at the root
        s.Yt[t] = Y_root
        s.μyμx_avt[t] = μyμx_avg_(s.smt[t], s.srst[t], s.snst[t],
                                   μymt[t], μyrst[t], μynst[t],
                                   μxmt[t], μxrst[t], μxnst[t])

        # Production wages: w_m^p = α μ̃^y Y / Lᵖ, w_j^p = (τ^p_j)^{1/α} w_m^p.
        wp_m = α * s.μy_avt[t] * s.Yt[t] / Lp_t
        s.wpmt[t]  = wp_m
        s.wprst[t] = τp_rs_t^(1.0/α) * wp_m
        s.wpnst[t] = τp_ns_t^(1.0/α) * wp_m

        # ---- Roy / Fréchet split of the fixed scientist pool Lˢ = ℓˢ·L ----
        # Common φ: headcount shares π_j ∝ (T_j^φ X_j^θ)^b, b = 1/[1+(θ-1)(1-φ)];
        # efficiency units H_j = π_j Lˢ ē_j with ē_j = Γ(1-1/θ) Φ^{1/θ}/w_j.
        # T_m ≡ 1; T_rs, T_ns recovered in closed form at t=1 to match the 1950
        # scientist-share targets.
        G_roy = _gammafn(1.0 - 1.0 / θfr)
        # Only the closed-form T recovery below needs a common φ. The allocation
        # itself handles sector-specific φ_j through the scalar fixed point.
        φ_roy = (alloc !== :stonegeary && t == 1 &&
                 T_rs_fixed === nothing && T_ns_fixed === nothing) ?
                _roy_require_common_phi(φm, φrs, φns) : φm

        Xm  = max(V_m,  0.0) * ζmt[t]  * s.Λmt[t]^γ
        Xrs = max(V_rs, 0.0) * ζrst[t] * s.Λrst[t]^γ
        Xns = max(V_ns, 0.0) * ζnst[t] * s.Λnst[t]^γ

        # Allocation-only wedge. Ṅ_j below still uses ζ_j and the resulting H_j, so
        # τ^s shifts who researches where without changing research productivity.
        Xm_a  = Xm
        Xrs_a = Xrs * max(tau_s_rs, 1e-12)
        Xns_a = Xns * max(tau_s_ns, 1e-12)

        if t == 1 && T_rs_fixed === nothing && T_ns_fixed === nothing
            # T_j = [ (π_j*/π_m*)^{D} (X_m/X_j)^{θ} ]^{1/φ}, D = 1+(θ-1)(1-φ).
            rs_sh = clamp(sci_sh_rs_target[1], 1e-6, 0.98)
            ns_sh = clamp(sci_sh_ns_target[1], 1e-6, 0.98)
            m_sh  = max(1e-6, 1.0 - min(rs_sh + ns_sh, 0.98))
            if Xm > 0.0
                if Xrs > 0.0
                    T_rs_const[] = roy_recover_T(rs_sh, m_sh, Xrs, Xm, φ_roy, θfr)
                end
                if Xns > 0.0
                    T_ns_const[] = roy_recover_T(ns_sh, m_sh, Xns, Xm, φ_roy, θfr)
                end
            end
        elseif t == 1 && (T_rs_fixed === nothing || T_ns_fixed === nothing)
            error("T_rs_fixed and T_ns_fixed must be supplied together")
        end

        local πm, πrs, πns, H_m, H_rs, H_ns, w_m, w_rs, w_ns
        if alloc === :data
            # DIAGNOSTIC MODE. The researcher allocation is IMPOSED at its observed
            # path rather than solved. Free entry then no longer holds with a common
            # wage: each sector's zero-profit wage is w_j = X_j H_j^{phi_j-1}, and
            # the DISPERSION of those wages is exactly the distortion the data
            # require. w_j/w_m is the implied wedge, read off rather than assumed.
            rs_sh = clamp(sci_sh_rs_target[min(t, length(sci_sh_rs_target))], 1e-8, 0.98)
            ns_sh = clamp(sci_sh_ns_target[min(t, length(sci_sh_ns_target))], 1e-8, 0.98)
            m_sh  = max(1.0 - rs_sh - ns_sh, 1e-8)
            tot   = m_sh + rs_sh + ns_sh
            πm, πrs, πns = m_sh/tot, rs_sh/tot, ns_sh/tot
            H_m, H_rs, H_ns = πm*Ls_tot, πrs*Ls_tot, πns*Ls_tot
            w_m  = (Xm_a  > 0.0 && H_m  > 0.0) ? Xm_a  * H_m^(φm  - 1.0) : 0.0
            w_rs = (Xrs_a > 0.0 && H_rs > 0.0) ? Xrs_a * H_rs^(φrs - 1.0) : 0.0
            w_ns = (Xns_a > 0.0 && H_ns > 0.0) ? Xns_a * H_ns^(φns - 1.0) : 0.0
        elseif alloc === :stonegeary
            # Homogeneous researchers: headcount IS efficiency units and there is a
            # single research wage, so π_j is just the headcount share.
            L_m, L_rs, L_ns, w_sg =
                sg_solve(Xm_a, Xrs_a, Xns_a, φm, φrs, φns,
                         Lbar_m, Lbar_rs, Lbar_ns, Ls_tot)
            Ltot = L_m + L_rs + L_ns
            πm  = Ltot > 0.0 ? L_m  / Ltot : 0.0
            πrs = Ltot > 0.0 ? L_rs / Ltot : 0.0
            πns = Ltot > 0.0 ? L_ns / Ltot : 0.0
            H_m, H_rs, H_ns = L_m, L_rs, L_ns
            w_m = w_rs = w_ns = w_sg
        else
            πm, πrs, πns, H_m, H_rs, H_ns, w_m, w_rs, w_ns =
                roy_solve(T_m, T_rs_const[], T_ns_const[], Xm_a, Xrs_a, Xns_a,
                          φm, φrs, φns, θfr, G_roy, Ls_tot)
        end

        # Headcount allocation and equilibrium research wages per efficiency unit.
        Ls_m  = πm  * Ls_tot
        Ls_rs = πrs * Ls_tot
        Ls_ns = πns * Ls_tot
        s.wsmt[t]  = w_m
        s.wsrst[t] = w_rs
        s.wsnst[t] = w_ns

        s.Lsmt[t]  = Ls_m
        s.Lsrst[t] = Ls_rs
        s.Lsnst[t] = Ls_ns
        s.Hmt[t]   = H_m
        s.Hrst[t]  = H_rs
        s.Hnst[t]  = H_ns
        s.Lst[t]   = Ls_m + Ls_rs + Ls_ns

        # Variety creation flows Ṅ_j = ζ_j Λ_j^γ H_j^{φ_j} (efficiency units).
        s.Nmdot[t]  = clamp(Ndot_sector(ζmt[t],  s.Λmt[t],  γ, H_m,  φm),  0.0, FLOW_CAP)
        s.Nrsdot[t] = clamp(Ndot_sector(ζrst[t], s.Λrst[t], γ, H_rs, φrs), 0.0, FLOW_CAP)
        s.Nnsdot[t] = clamp(Ndot_sector(ζnst[t], s.Λnst[t], γ, H_ns, φns), 0.0, FLOW_CAP)

        bad = !isfinite(s.Yt[t]) || !isfinite(s.pmt[t]) ||
              !isfinite(s.Lst[t]) || s.Lst[t] < 0
        return !bad
    end

    # ----------------------------------------------------------------------
    # Outer fixed-point on V_j (PDV of profits)
    # ----------------------------------------------------------------------
    iter_v = 1
    diff_V = 1e5
    diff_V_hist = Float64[]
    Vmt    = (V_init_m  !== nothing && length(V_init_m)  == T) ? copy(V_init_m)  :
             (length(V_SEED_M)  == T ? copy(V_SEED_M)  : fill(0.5, T))
    Vrst   = (V_init_rs !== nothing && length(V_init_rs) == T) ? copy(V_init_rs) :
             (length(V_SEED_RS) == T ? copy(V_SEED_RS) : fill(0.5, T))
    Vnst   = (V_init_ns !== nothing && length(V_init_ns) == T) ? copy(V_init_ns) :
             (length(V_SEED_NS) == T ? copy(V_SEED_NS) : fill(0.5, T))
    θm     = something(theta_m_fixed, θm_in)
    θrs    = something(theta_rs_fixed, θrs_in)
    θns    = 1.0 - θm - θrs
    diff_θm  = 0.0
    diff_θrs = 0.0

    # Helper to return failure dict with diagnostics
    @inline function _populate_diag_and_return(fail_t::Int, fail_code::Int)
        d = state_to_dict(s)
        d["Vmt"]  = copy(Vmt);  d["Vrst"] = copy(Vrst); d["Vnst"] = copy(Vnst)
        d["converged"] = [0.0]
        d["fail_t"]    = [Float64(fail_t)]
        d["fail_code"] = [Float64(fail_code)]
        avg_pm = nfev_pm[] > 0 ? iter_pm_total[] / nfev_pm[] : 0.0
        avg_y  = nfev_y[]  > 0 ? iter_y_total[]  / nfev_y[]  : 0.0
        d["diag_n_v_iters"]   = [Float64(n_v_iters[])]
        d["diag_nfev_pm"]     = [Float64(nfev_pm[])]
        d["diag_nfev_y"]      = [Float64(nfev_y[])]
        d["diag_nfev_lp"]     = [0.0]
        d["diag_avg_iter_pm"] = [avg_pm]
        d["diag_avg_iter_y"]  = [avg_y]
        d["diag_avg_iter_lp"] = [0.0]
        d["diag_nfail_lp"]    = [0.0]
        d["diag_diff_V"]      = [diff_V]
        d["diff_θm"]          = [diff_θm]
        d["diff_θrs"]         = [diff_θrs]
        # Populate on the failure path too, so downstream printing never KeyErrors.
        d["θm"]               = [θm]
        d["θrs"]              = [θrs]
        d["θns"]              = [1.0 - θm - θrs]
        d["τpmt"]             = ones(T)
        d["τprst"]            = τprst
        d["τpnst"]            = τpnst
        d["T_m"]              = [T_m]
        d["T_rs"]             = [T_rs_const[]]
        d["T_ns"]             = [T_ns_const[]]
        d["θfr"]              = [θfr]
        return d
    end

    # Without fixed θ, recover the weights in closed form at t=1 to match the
    # data VA-share targets. With fixed θ, period-one shares remain endogenous.
    if !fixed_theta && !(isfinite(smt_1_target) && isfinite(srst_1_target) &&
         smt_1_target > 0.0 && srst_1_target > 0.0 &&
         (smt_1_target + srst_1_target) < 1.0)
        return _populate_diag_and_return(1, 3)
    end

    while (iter_v <= iter_max_v) && (diff_V > tol_v)
        n_v_iters[] += 1

        # ---- t=1: treat first period exactly as all other periods ----
        s.Lt[1]   = Lt[1]
        s.Nmt[1]  = ini_Nm
        s.Nrst[1] = ini_Nrs
        s.Nnst[1] = ini_Nns
        refresh_lambda!(1)

        θm1, θrs1, θns1 = path_theta ? (theta_m_path[1], theta_rs_path[1], theta_ns_path[1]) : (θm, θrs, θns)
        ok = solve_period_inner!(1, Vmt[1], Vrst[1], Vnst[1], θm1, θrs1, θns1, 5.0;
                 match_t1_shares=!fixed_theta && !path_theta)
        if !ok
            return _populate_diag_and_return(1, 1)
        end

        # ---- Recover θ in closed form from the t=1 equilibrium when requested ----
        # With imposed shares s_j(1)=s_j* and raw (pre-numeraire) prices p_j(1),
        #   θ_j = s_j* / (p_j(1)^(1-σ) · Y(1)^((η_j-1)(1-σ))),
        # which sums to one by construction (see solve_pm_given_Y!). These θ are
        # then used for the forward shoot t=2..T.
        if !fixed_theta && !path_theta
            θm_prev  = θm
            θrs_prev = θrs
            let εθ = 1.0 - σ, Y1 = s.Yt[1]
                if εθ == 0.0
                    θm  = clamp(smt_1_target,  1e-12, 1.0)
                    θrs = clamp(srst_1_target, 1e-12, 1.0)
                    θns = max(1.0 - θm - θrs, 1e-12)
                else
                    sm_star  = clamp(smt_1_target,  1e-8, 1.0 - 2e-8)
                    srs_star = clamp(srst_1_target, 1e-8, 1.0 - 2e-8)
                    sns_star = max(1.0 - sm_star - srs_star, 1e-8)
                    wm  = s.pmt[1]^εθ  * Y1^((ηm  - 1.0)*εθ)
                    wrs = s.prst[1]^εθ * Y1^((ηrs - 1.0)*εθ)
                    wns = s.pnst[1]^εθ * Y1^((ηns - 1.0)*εθ)
                    θm  = sm_star  / max(wm,  1e-300)
                    θrs = srs_star / max(wrs, 1e-300)
                    θns = sns_star / max(wns, 1e-300)
                    Σθ  = θm + θrs + θns
                    θm /= Σθ; θrs /= Σθ; θns /= Σθ
                end
            end
            diff_θm  = abs(θm  - θm_prev)
            diff_θrs = abs(θrs - θrs_prev)
        else
            diff_θm  = 0.0
            diff_θrs = 0.0
        end

        # Update varieties for t=2
        s.Nmt[2]  = clamp(s.Nmt[1]  + s.Nmdot[1],  1e-12, N_CAP)
        s.Nrst[2] = clamp(s.Nrst[1] + s.Nrsdot[1], 1e-12, N_CAP)
        s.Nnst[2] = clamp(s.Nnst[1] + s.Nnsdot[1], 1e-12, N_CAP)
        refresh_lambda!(2)

        # ---- forward shoot t=2..T ----
        for t in 2:T
            s.Lt[t] = Lt[t]
            if lock_N && isfinite(N_m_lock[t]) && isfinite(N_rs_lock[t]) && isfinite(N_ns_lock[t])
                s.Nmt[t] = N_m_lock[t]; s.Nrst[t] = N_rs_lock[t]; s.Nnst[t] = N_ns_lock[t]
                all(>(0.0), (s.Nmt[t], s.Nrst[t], s.Nnst[t])) || error("locked knowledge stocks must be positive")
                refresh_lambda!(t)
            end
            Y_ref = s.Yt[t-1]
            θmt, θrst, θnst = path_theta ? (theta_m_path[t], theta_rs_path[t], theta_ns_path[t]) : (θm, θrs, θns)
            ok = solve_period_inner!(t, Vmt[t], Vrst[t], Vnst[t], θmt, θrst, θnst, Y_ref)
            if !ok
                return _populate_diag_and_return(t, 1)
            end

            if t < T
                s.Nmt[t+1]  = clamp(s.Nmt[t]  + s.Nmdot[t],  1e-12, N_CAP)
                s.Nrst[t+1] = clamp(s.Nrst[t] + s.Nrsdot[t], 1e-12, N_CAP)
                s.Nnst[t+1] = clamp(s.Nnst[t] + s.Nnsdot[t], 1e-12, N_CAP)
                refresh_lambda!(t+1)
            end
        end

        # ---- Vectorized post-processing of equilibrium objects ----
        s.Ymt  .= Y_demand_.(s.smt,  s.Yt, s.pmt)
        s.Yrst .= Y_demand_.(s.srst, s.Yt, s.prst)
        s.Ynst .= Y_demand_.(s.snst, s.Yt, s.pnst)

        # Here s_rs and s_ns are already expenditure shares of total Y.
        s.Lpmt  .= Lp_conditional_.(s.smt,  s.Yt, s.pmt,  s.Nmt,  s.e_Y_mt,  s.Amt,  α)
        s.Lprst .= Lp_conditional_.(s.srst, s.Yt, s.prst, s.Nrst, s.e_Y_rst, s.Arst, α)
        s.Lpnst .= Lp_conditional_.(s.snst, s.Yt, s.pnst, s.Nnst, s.e_Y_nst, s.Anst, α)

        # Per-variety machine quantities and profits
        s.xmt  .= x_sym_.(μymt,  μxmt,  zmt,  s.pmt,  s.Ymt,  s.Nmt,  α)
        s.xrst .= x_sym_.(μyrst, μxrst, zrst, s.prst, s.Yrst, s.Nrst, α)
        s.xnst .= x_sym_.(μynst, μxnst, znst, s.pnst, s.Ynst, s.Nnst, α)
        s.πxmt .= π_x_.(μymt,  μxmt,  s.pmt,  s.Ymt,  s.Nmt,  α)
        s.πxrst.= π_x_.(μyrst, μxrst, s.prst, s.Yrst, s.Nrst, α)
        s.πxnst.= π_x_.(μynst, μxnst, s.pnst, s.Ynst, s.Nnst, α)

        # Aggregate consumption from final-good market clearing
        s.Ct .= C_.(s.Yt, s.μyμx_avt, α)

        # Real interest rate from Euler equation (per-capita c growth)
        @inbounds for t in 2:T
            ct   = s.Ct[t]   / s.Lt[t]
            ctm1 = s.Ct[t-1] / s.Lt[t-1]
            s.rt[t-1] = r_euler_(ρ, χ, ct, ctm1)
        end

        # Normalize all sector prices to one at t=1 (numeraire choice)
        pm0  = s.pmt[1]
        prs0 = s.prst[1]
        pns0 = s.pnst[1]
        if isfinite(pm0) && isfinite(prs0) && isfinite(pns0) && pm0 > 0 && prs0 > 0 && pns0 > 0
            s.pmt  ./= pm0
            s.prst ./= prs0
            s.pnst ./= pns0
        else
            return _populate_diag_and_return(1, 2)
        end

        # θ (θ_m, θ_rs, θ_ns) was recovered in closed form at t=1 above; its
        # V-iteration change is tracked in diff_θm / diff_θrs.

        # ---- update value functions: V_j as PDV of per-variety profits ----
        H = 500 - T   # forward-extension horizon for BGP tail of PDV
        pad       = H
        πxmt_ext  = vcat(s.πxmt,  fill(s.πxmt[end],  pad))
        πxrst_ext = vcat(s.πxrst, fill(s.πxrst[end], pad))
        πxnst_ext = vcat(s.πxnst, fill(s.πxnst[end], pad))
        r_last    = s.rt[T - 1]
        rt_ext    = vcat(s.rt[1:T-1], fill(r_last, pad + 1))

        Vmt_new  = similar(Vmt)
        Vrst_new = similar(Vrst)
        Vnst_new = similar(Vnst)

        cumsum_rt = vcat(0.0, cumsum(rt_ext))
        disc_buf  = zeros(H)

        @inbounds for t in 1:T
            for k in 1:H
                disc_buf[k] = (k == 1) ? 1.0 : exp(-(cumsum_rt[t + k - 1] - cumsum_rt[t]))
            end
            vm_acc  = 0.0
            vrs_acc = 0.0
            vns_acc = 0.0
            for k in 1:H
                vm_acc  += πxmt_ext[t + k - 1]  * disc_buf[k]
                vrs_acc += πxrst_ext[t + k - 1] * disc_buf[k]
                vns_acc += πxnst_ext[t + k - 1] * disc_buf[k]
            end
            Vmt_new[t]  = vm_acc
            Vrst_new[t] = vrs_acc
            Vnst_new[t] = vns_acc
        end

        diff_V = maximum(abs.(vcat(Vmt_new - Vmt, Vrst_new - Vrst, Vnst_new - Vnst)))
        push!(diff_V_hist, diff_V)

        # Damped V update
        Vmt  = ω_v .* Vmt_new  .+ (1.0 - ω_v) .* Vmt
        Vrst = ω_v .* Vrst_new .+ (1.0 - ω_v) .* Vrst
        Vnst = ω_v .* Vnst_new .+ (1.0 - ω_v) .* Vnst

        s.Vmt  .= Vmt
        s.Vrst .= Vrst
        s.Vnst .= Vnst

        iter_v += 1
    end

    # ---- Pack final state into Dict with metadata ----
    d = state_to_dict(s)
    d["θm"]        = [θm]
    d["θrs"]       = [θrs]
    d["θns"]       = [1.0 - θm - θrs]
    d["diff_θm"]   = [diff_θm]
    d["diff_θrs"]  = [diff_θrs]
    d["converged"] = [(diff_V <= tol_v) ? 1.0 : 0.0]
    d["diag_diff_V_hist"] = diff_V_hist
    d["fail_t"]    = [0.0]
    d["fail_code"] = [0.0]

    avg_pm = nfev_pm[] > 0 ? iter_pm_total[] / nfev_pm[] : 0.0
    avg_y  = nfev_y[]  > 0 ? iter_y_total[]  / nfev_y[]  : 0.0
    d["diag_n_v_iters"]   = [Float64(n_v_iters[])]
    d["diag_nfev_pm"]     = [Float64(nfev_pm[])]
    d["diag_nfev_y"]      = [Float64(nfev_y[])]
    d["diag_nfev_lp"]     = [0.0]
    d["diag_avg_iter_pm"] = [avg_pm]
    d["diag_avg_iter_y"]  = [avg_y]
    d["diag_avg_iter_lp"] = [0.0]
    d["diag_nfail_lp"]    = [0.0]
    d["diag_diff_V"]      = [diff_V]
    d["diff_θm"]          = [diff_θm]
    d["diff_θrs"]         = [diff_θrs]
    d["τpmt"]             = ones(T)
    d["τprst"]            = τprst
    d["τpnst"]            = τpnst
    d["T_m"]              = [T_m]
    d["T_rs"]             = [T_rs_const[]]
    d["T_ns"]             = [T_ns_const[]]
    d["θfr"]              = [θfr]

    return d
end

#------------------------------------------------------------------------------
# 6. Driver (only runs when this file is executed directly)
#------------------------------------------------------------------------------

if abspath(PROGRAM_FILE) == @__FILE__

    T = 100

    # --- Household preferences ---
    χ  = 2.000
    ρ  = 0.020

    # --- Knowledge spillover matrix (rows: receiving; cols: source; rows sum to 1) ---
    δmm = 0.920;  δmr = 0.050;  δmn = 0.030
    δrm = 0.570;  δrr = 0.270;  δrn = 0.160
    δnm = 0.720;  δnr = 0.100;  δnn = 0.180

    # --- Structural parameters ---
    α   = 0.70;   κ   = (m=0.866, rs=0.866, ns=0.866)
    σ   = 0.010
    ηm  = 0.505;  ηrs = 1.573;  ηns = 1.000
    γ   = 0.65;   φ   = 0.222
    θm  = 0.4144; θrs = 0.0667

    # --- Innovation productivity ---
    gζm = 0.0;  gζrs = 0.0;  gζns = 0.0
    ζm0 = 0.2;  ζrs0 = 0.2;  ζns0 = 0.2

    ini_Nm = 1.0;  ini_Nrs = 1.0;  ini_Nns = 1.0

    # --- Population / labor ---
    gl = fill(0.013, T - 1)
    Lt = cumprod([1.0; 1.0 .+ gl])

    # --- Exogenous paths (constant by default) ---
    μxmt   = fill(0.769, T);  μxrst = fill(0.677, T);  μxnst = fill(0.677, T)
    μymt   = fill(0.90, T);   μyrst = fill(0.90, T);   μynst = fill(0.90, T)
    zmt    = fill(1.00, T);   zrst  = fill(1.00, T);   znst  = fill(1.00, T)
    hbarmt = fill(1.00, T);   hbarrst = fill(1.00, T); hbarnst = fill(1.00, T)

    # --- Exogenous research share ℓˢ(t) and calibration targets ---
    ls_path          = fill(0.05, T)
    sci_sh_rs_target = fill(0.33, T)
    sci_sh_ns_target = fill(0.33, T)
    p_ratio_rs       = fill(1.0, T)
    p_ratio_ns       = fill(1.0, T)

    @time eq = solve_GE(
        T, χ, ρ,
        δmm, δmr, δmn, δnm, δnr, δnn, δrm, δrr, δrn,
        α, κ,
        σ, ηm, ηrs, ηns, γ, φ, θm, θrs,
        gζm, gζrs, gζns,
        ζm0, ζrs0, ζns0,
        ls_path,
        sci_sh_rs_target, sci_sh_ns_target,
        p_ratio_rs, p_ratio_ns,
        0.42, 0.07,
        ini_Nm, ini_Nrs, ini_Nns, gl,
        μxmt, μxrst, μxnst, μymt, μyrst, μynst,
        zmt, zrst, znst, hbarmt, hbarrst, hbarnst, Lt;
        mode=:calibrate,
        lambda_N0=(m=1/3, rs=1/3, ns=1/3),
        θfr=2.0, T_m=1.0,
    )

    println()
    println("converged = ", eq["converged"][1])
    if haskey(eq, "fail_t") && eq["fail_t"][1] > 0
        println("failure diagnostic: t=", Int(round(eq["fail_t"][1])),
                " code=", Int(round(eq["fail_code"][1])))
    end
    println()
    println("─── Solver diagnostics ────────────────────────────────────────────")
    @printf "  V iterations           : %d\n"   Int(eq["diag_n_v_iters"][1])
    @printf "  p_m calls (closed-form): %d\n"    Int(eq["diag_nfev_pm"][1])
    @printf "  Y  root-finder calls   : %d   (avg %.1f iters/call)\n" Int(eq["diag_nfev_y"][1]) eq["diag_avg_iter_y"][1]
    println()

    println("─── Equilibrium path (selected periods) ───────────────────────────")
    @printf "%-4s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s\n" "t" "s_m" "s_rs" "Y" "p_m" "Ls/L" "w^p_m" "w^s_rs" "π_rs"
    println("-"^88)
    for t in [1; 10:10:T]
        Ls_share = eq["Lst"][t] / max(1e-12, eq["Lt"][t])
        π_rs = eq["Lsrst"][t] / max(1e-12, eq["Lst"][t])
        @printf "%-4d  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f\n" t eq["smt"][t] eq["srst"][t] eq["Yt"][t] eq["pmt"][t] Ls_share eq["wpmt"][t] eq["wsrst"][t] π_rs
    end
    println()
    @printf "  Fréchet scales: θ_F=%.3f  T_m=%.4f  T_rs=%.4f  T_ns=%.4f\n" eq["θfr"][1] eq["T_m"][1] eq["T_rs"][1] eq["T_ns"][1]
end
