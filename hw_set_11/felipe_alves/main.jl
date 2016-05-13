##############################################################################
##
## Testing functions
##
##############################################################################
using PyPlot

include("hw11.jl")

βvals = collect(0.94:0.01:0.98)
pps = PlannerProb[PlannerProblem(β;  sig = 0.2, mu =-0.1, grid_size = 100)
                for β in βvals]

xmax = 8.0
# ygrid = linspace(0.01, xmax, 500)
ygrid = linspace(log(0.01), log(xmax), 500)

fig, axes = subplots(3,2,figsize = (8,8))

fig2 = figure()
ax2 = fig2[:gca](projection="3d")

ax2[:set_xlim3d](βvals[1],βvals[end])
ax2[:set_xticks](βvals)
ax2[:set_xlabel](L"\beta" ,fontsize = 20)
# ax2[:set_yticks](ygrid[1],ygrid[2])

# ax2[:set_ylim3d](0.20, exp(1.5))
ax2[:set_ylim3d](-1.5, 1.5)
ax2[:set_ylabel](L"y" ,fontsize = 20)
ax2[:set_zlim3d](0, 2)
# ax2[:set_zticks]((0.2, 0.4))

ind_lw  = searchsortedfirst(ygrid,-1.5)-1
ind_hi = searchsortedfirst(ygrid,1.5)

for i = 1:length(βvals)
    println("β = $(βvals[i]) ")
    pp = pps[i]

    #== Compute policy ==#
    compute_σ_howard(pp; use_quad = false)
    #== Asymptotic dist ==#
    # ψ = asymp_liny_ψ(pp, ygrid)
    ψ = asymp_logy_ψ(pp, ygrid)

    #== Policies figures ==#
    if i==1 || i ==length(βvals)
        fig0, ax0 = subplots(figsize = (10,10))
        grid = collect(pp.grid)
        ypri = mean(pp.phi)*pp.f(grid- pp.σ)
        ax0[:plot](grid, pp.σ, lw = 3, label = L"\sigma(y)")
        ax0[:plot](grid, ypri, lw = 3, label = L"\bar{\xi}f(y-\sigma(y))")
        ax0[:plot](grid, grid, "k", lw = 1)

        ax0[:legend](ncol = 2, fontsize = 15)
        ax0[:set_xlim]([0, 4.0])
        ax0[:set_ylim]([0, 3.0])
    end

    ax = axes[i]
    # ax[:set_xlim](0.0,exp(1.5))
    ax[:set_xlim](-1.5,1.5)
    ax[:set_title](latexstring("\\beta = $(βvals[i])"))
    ax[:set_xlabel]("y")
    ax[:plot](ygrid, ψ, lw=5, alpha=0.6)

    ax2[:plot3D](ygrid[ind_lw:ind_hi], ψ[ind_lw:ind_hi], βvals[i], zdir="x", color = "blue", alpha = 0.6, lw = 5)

    println("..................................................")
end
# ax[:legend](ncol = 2)

pp = pps[5]
#== Compute policy ==#
compute_σ_howard(pp; use_quad = false)
#== Asymptotic dist ==#
ψ = asymp_ψ(pp, ygrid)
