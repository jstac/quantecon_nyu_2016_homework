##############################################################################
##
## HW set 11
##
##############################################################################

# Model
# y_{t+1} = ξ_{t+1}f( y_t - c_t )

## Algorithm

# 1. Pick a guess σ¹
# 2. Evaluate v_σ¹
# 3. Choose σ² from optimal poliy from Bellman

# ..................................................................

using Distributions
using QuantEcon: qnwlogn, qnwlege, do_quad, compute_fixed_point, LAE, lae_est
using Interpolations
using Optim: optimize, Brent, GoldenSection
# using NLopt

function logistic(x, a =1, b=2, c = 20, d = 1)
    a + (b-a) * 1./( 1 + exp(-c * (x-d) ) )
end

"""
Planner Problem type

 max_{c} E ∑ u(c_t),    y_{t+1} = ξ_{t+1}f( y_t - c_t )

#### Fields
- `β`               : discount factor
- `u`               : utility function
- `f`               : growth function
- `grid`            : grid of y

- `phi`             : distribution of ξ
- `quad_nodes`      :
- `quad_weights`    :
- `ξ_sample`        : sample for MonteCarlo comp of Exp

- `v`               : value function storage
- `σ`               : policy function storage

"""
type PlannerProb

    β::Float64
    u::Function
    f::Function

    grid::AbstractVector

    phi::Distribution
    quad_nodes::Vector
    quad_weights::Vector
    ξ_samp::Vector

    v::Vector{Float64}
    σ::Vector{Float64}

end

function PlannerProblem(β; γ = 0.9, θ = 0.5, α = 0.5, mu =-0.1, sig = sqrt(0.2),
                        grid_max = 8, grid_size = 120)

# γ = 0.9; θ = 0.5; α = 0.5; mu =-0.1; sig = 0.2
# grid_max = 30; grid_size = 200
    #== Functions ==#
    u(c) = 1 - exp(-θ * c.^γ )
    f(s) = (s.^α) .* (logistic(s))

    #== Grid ==#
    # grid1 = linspace(1e-8, 3.0, 3*div(grid_size,4))
    # grid2 = linspace(3.0, grid_max, div(grid_size,4)+1)
    # grid  = [collect(grid1); collect(grid2[2:end])]

    # grid = powerspacegrid(1e-5, grid_max, grid_size, 0.5)
    grid = linspace(1e-8, grid_max, grid_size)

    #== ξ distribution and nodes for ∫ ==#
    phi     = LogNormal(mu, sig)
    sig_phi = std(phi)
    _int_min, _int_max = 0.5, 3.5 #exp(mu + 3 * sig)
    n, w = qnwlege(30, _int_min, _int_max)

    ξ_samp = rand(phi,1000)

    #== Value function & policy ==#
    v = zeros(grid)
    σ = copy(grid)

    PlannerProb(β, u, f, grid, phi, n, w, ξ_samp, v, σ)
end


function integrate(pp::PlannerProb, g::Function, use_quad::Bool,
                   qn::Array=pp.quad_nodes, qw::Array=pp.quad_weights)


    if use_quad
        #== Compute Integral with quad ==#
        int_func(x) = g(x) .* pdf(pp.phi, x)
        return do_quad(int_func, qn, qw)
    else
        #== Compute Integral with Mc ==#
        ξ = pp.ξ_samp
        return mean(g(ξ))
        # == Deterministic case ==#
        # return g(1.0)

    end
end

function compare_int!(pp, out)

    v, σ = pp.σ, pp.v
    f    = pp.f
    grid = pp.grid

    #== Method 1 ==#
    itp  = interpolate((collect(grid),), v, Gridded(Linear()) )
    # itp2 = extrapolate(itp, Flat())
    Lv(y)    = itp[y]

    for (i,x) in enumerate(grid)
        to_integrate(ξ) = Lv( ξ * f(x) )
        out[i,1] = integrate(pp, to_integrate)
    end

    #== Method 2 ==#
    ξ_samp = pp.ξ_samp
    N = length(ξ_samp)

    for (i,x) in enumerate(grid)
        sum = 0.0
        for ξ in ξ_samp
            sum +=  Lv( ξ * f(x) )
        end
        out[i,2] = sum/N
    end

    #== Method 3 ==#
    # mu   = mean(pp.phi)
    # sig2 = var(pp.phi)
    # nodes, weights = qnwlogn(21, -0.1, 0.2)
    # for (i,x) in enumerate(grid)
    #     out[i,3] = dot( Lv( nodes * f(x) ), weights )
    # end
    return Void
end


function bellman_operator(pp::PlannerProb, v::Vector; ret_policy=false, use_quad = false)

    ret_policy ? (out = zeros(pp.σ)) : (out = zeros(v))
    bellman_operator!(pp, v, out, ret_policy=ret_policy, use_quad = use_quad)
    return out
end

function bellman_operator!(pp, v, out; ret_policy = false, use_quad = false)

    u, f, β = pp.u, pp.f, pp.β
    grid = pp.grid

    c_min = 0.0

    # Lv(y) = interpolate((collect(grid),), v, Gridded(Linear()) )[y]
    itp  = interpolate((collect(grid),), v, Gridded(Linear()) )
    itp2 = extrapolate(itp, Flat())
    Lv(y)    = itp2[y]

    for (i,x) in enumerate(grid)
        function obj(c)
            to_integrate(ξ) = Lv( ξ * f(x - c) )
            return -u(c)  - β * integrate(pp, to_integrate, use_quad)
        end
        rgl = optimize(obj, c_min, x, method = GoldenSection())
        rbr  = optimize(obj, c_min, x, method = Brent())

        (rbr.f_minimum < rgl.f_minimum) ? (res = rbr) : (res = rgl)
        # res = rgl

        if ret_policy
            out[i] =   res.minimum
        else
            out[i] = - res.f_minimum
        end
    end
end

function rhs(pp, x, c, use_quad = false)

    u, f, β = pp.u, pp.f, pp.β
    grid = collect(pp.grid)

    # Lv(y) = interpolate((grid,), pp.v, Gridded(Linear()) )[y]
    itp  = interpolate((collect(grid),), pp.v, Gridded(Linear()) )
    itp2 = extrapolate(itp, Flat())
    Lv(y)    = itp2[y]

    to_integrate(ξ) = Lv( ξ * f(x - c) )
    return u(c) + β * integrate(pp, to_integrate, use_quad)
    # return β * integrate(pp, to_integrate, doquad)

end
"""
Iteration on the policy function operator


"""
function policy_operator!(pp::PlannerProb, v::AbstractVector, Tv::AbstractVector; use_quad::Bool=false)

    u, f, β  = pp.u, pp.f, pp.β
    grid, σ  = pp.grid, pp.σ

    # grid_int = collect(grid)
    itp  = interpolate((collect(grid),), v, Gridded(Linear()) )
    itp2 = extrapolate(itp, Flat())
    Lv(y)    = itp2[y]

    for (i, x) in enumerate(grid)
        to_integrate(ξ) = Lv( ξ * f( x - σ[i] ) )
        Tv[i] = u( σ[i] ) + β * integrate(pp, to_integrate, use_quad)
    end

end


function policy_operator(pp::PlannerProb, v::AbstractVector; use_quad::Bool = false)

    Tv = zeros(v)
    policy_operator!(pp, v, Tv; use_quad = use_quad)
    return Tv
end

"""
Compute the policy function for the PlannerProblem `pp`. Results are
updated on the instance.

##### Arguments

- `pp` : instance of PlannerProbl


"""
function compute_σ_howard(pp::PlannerProb; use_quad::Bool = false, pol_tol = 1e-6, max_it = 100)


    #== Use quadrature to compute Exp ==#

    #== Bellman oper ==#
    bellman_oper0(v, pol) = bellman_operator(pp, v; ret_policy=pol, use_quad = true)
    bellman_oper(v, pol)  = bellman_operator(pp, v; ret_policy=pol, use_quad = use_quad) # set ret policy true

    #== Values ==#
    v_init  = copy(pp.v)  # Initial condition
    v_new   = copy(pp.v)

    #== Policies ==#
    # σ_old  = copy(pp.σ)

    it = 1
    dist = 1

    #== Iterate a little on Bellman ==#
    # println("Initial iterations on Bellman")
    for j =1:100
        v_init = bellman_oper0(v_init,false)
    end

    copy!(pp.σ, bellman_oper0(v_init, true))

    while dist>pol_tol && it<=max_it

        # @printf("ITERATION %02d \n", it)
        # println("--------------------------------------------------------------")
        #== Policy operator ==#
        pol_oper(v_σ) = policy_operator(pp, v_σ, use_quad = use_quad)  # defined inside for bc pp.σ updates

        #== Policy value function ==#
        vσ = compute_fixed_point(pol_oper, v_init; err_tol=1e-8, max_iter=1000, verbose=false, print_skip=50)

        σ_new = bellman_oper(vσ, true)
        v_new = bellman_oper(vσ, false)

        #== comment ==#
        # dist = maximum(abs(pp.σ - σ_new))
        dist = maximum(abs(v_init - v_new))

        #== Update on pp ==#
        copy!(pp.σ, σ_new)
        copy!(v_init,vσ)
        @printf("Iteration %02d with error : %.6f\n", it, dist)
        # println("..............................................................\n")
        it += 1
    end

    #== Update value function on pp ==#
    copy!(pp.v, v_new)

    return Void
end

function compute_σ_vfi(pp::PlannerProb)

    bellman_oper(v) = bellman_operator(pp, v) # set ret policy false

    v_init = zeros(pp.v)
    v_star = compute_fixed_point(bellman_oper, v_init; err_tol=1e-6, max_iter=500, verbose=false, print_skip=10)
    σ_star = bellman_operator(pp,v_star; ret_policy = true)

    copy!(pp.v, v_star)
    copy!(pp.σ, σ_star)

    return Void
end

"""
1/n ∑ p(x,y) → ∫ p(x,y)ψ(x) = ψ₁(y)
"""
function asymp_liny_ψ(pp::PlannerProb,ygrid::AbstractVector)

    grid, σ = pp.grid, pp.σ
    phi     = pp.phi

    #== Interpolation ==#
    σpol(y)  = interpolate((collect(grid),), σ, Gridded(Linear()) )[y]
    function p(x, y)
        #=
        Stochastic kernel for model.
        Both x and y must be strictly positive.
        =#
        d = pp.f( x - σpol(x))

        # scipy silently evaluates the pdf of the lognormal dist at a negative
        # value as zero. It should be undefined and Julia recognizes this.
        pdf_arg = clamp(y ./ d, eps(), Inf)
        return pdf(phi, pdf_arg)  ./ d
    end

    n = 1000    # Number of observations at each date t
    T = 51      # Compute density of k_t at 1,...,T+1

    # Generate matrix s.t. t-th column is n observations of k_t
    y = Array(Float64, n, T)
    ξ = rand!(phi, Array(Float64, n, T))

    # Draw first column from initial distribution
    # match scale=0.5 and loc=2*i in python version
    y[:, 1] = 1.0
    for t=1:T-1
        y[:, t+1] = ξ[:, t] .* pp.f( y[:,t] - σpol(y[:,t]) )
    end
    # Generate T instances of LAE using this data, one for each date t
    laes = LAE(p, y[:, T])
    # laes = [LAE(p, y[:, t]) for t=T:-5:1]

    # fig, ax = subplots(figsize = (8,8))
    # greys = [string(g) for g in linspace(0.0, 1.0, T)]
    # for (psi, g) in zip(laes, greys)
        # ax[:plot](ygrid, lae_est(psi, ygrid), color=g, lw=4, alpha=0.6)
    # end
    # ax[:set_xlabel]("y")
    # ax[:set_xlim]([0,2.0])
    return lae_est(laes, ygrid)
end

"""
1/n ∑ p(x,y) → ∫ p(x,y)ψ(x) = ψ₁(y)
"""
function asymp_logy_ψ(pp::PlannerProb,ygrid::AbstractVector)

    grid, σ = pp.grid, pp.σ
    phi     = pp.phi
    μ,sig   = params(phi)
    normal = Normal(μ, sig)
    #== Interpolation ==#
    σpol(y)  = interpolate((collect(grid),), σ, Gridded(Linear()) )[y]

    function p(x, y)
        #=
        Stochastic kernel for model.
        =#
        pdf_arg = y .- log( pp.f( exp(x) - σpol(exp(x))) )
        return pdf(normal, pdf_arg)
    end

    n = 1000    # Number of observations at each date t
    T = 51      # Compute density of k_t at 1,...,T+1

    # Generate matrix s.t. t-th column is n observations of k_t
    y = Array(Float64, n, T)
    ξ = rand!(phi, Array(Float64, n, T))

    # Draw first column from initial distribution
    # match scale=0.5 and loc=2*i in python version
    y[:, 1] = 1.0
    for t=1:T-1
        y[:, t+1] = ξ[:, t] .* pp.f( y[:,t] - σpol(y[:,t]) )
    end
    # Generate T instances of LAE using this data, one for each date t
    laes = LAE(p, log(y[:, T]) )
    # laes = [LAE(p, y[:, t]) for t=T:-5:1]

    # fig, ax = subplots(figsize = (8,8))
    # greys = [string(g) for g in linspace(0.0, 1.0, T)]
    # for (psi, g) in zip(laes, greys)
        # ax[:plot](ygrid, lae_est(psi, ygrid), color=g, lw=4, alpha=0.6)
    # end
    # ax[:set_xlabel]("y")
    # ax[:set_xlim]([0,2.0])
    return lae_est(laes, ygrid)
end


##############################################################################
##
## Helping function
##
##############################################################################

function powerspacegrid(init::Real, eend::Real, n::Int64, k::Real = 1)

    if n<=2
        println("n has to be larger than 2")
        return
    end

    x = collect(linspace(0,1,n))

    z = x.^(1.0/float(k))

    return init +(eend - init)*z
end
