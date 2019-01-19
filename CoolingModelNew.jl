1+1
## ---

using StatGeochem, Statistics
include("UtilitiesNew.jl")
using Plots; gr(); default(fmt=:svg);
using ProgressMeter: @showprogress

## --- Define physical parameters and set up model

# Time step (Myr)
dt = 0.1
timevec = collect(0:dt:4000)

# Calculate past mantle heat production
p = Parameters()
H = bse_heat_production(p,mcdbse,timevec) * p.Mm
Qrm = H .* (p.Qrm_now / H[1])
## --- Run Newton's method

nNewton = 20 # Number of Newton's method iterations
dQrm_now = p.Qrm_now*1E-9 # Heat production step size
p.Rc = 650E3 # Plate bending radius (km)

Qrm_now_vec = Array{Float64}(undef, nNewton)
Qrm_now_vec[1] = p.Qrm_now # Initial guess for modern heat production
for i = 1:(nNewton-1)
    ll = coolingmodel_LL(p, Qrm * (Qrm_now_vec[i]/Qrm[1]), dt, timevec, TrelObs_time*1000, TrelObs, TrelObs_sigma)
    llp = coolingmodel_LL(p, Qrm * ((Qrm_now_vec[i]+dQrm_now)/Qrm[1]), dt, timevec, TrelObs_time*1000, TrelObs, TrelObs_sigma)
    dll_dQrm = (llp-ll) / dQrm_now
    Qrm_now_vec[i+1] = Qrm_now_vec[i] - ll/dll_dQrm/2
end

(Tm, Tc, Qm, Qc, dp) = coolingmodel(p, Qrm * (Qrm_now_vec[nNewton]/Qrm[1]), dt, timevec)
h = plot(TrelObs_time, TrelObs, yerror=2*TrelObs_sigma, seriestype=:scatter, markersize=2, color=:darkblue, markerstrokecolor=:auto, label="Data")
plot!(h, timevec/1000, Tm .- Tm[1], label="Model", xlabel="Age (Ma)", ylabel="\\Delta T (C)", framestyle=:box, legend=:topleft)
xlims!(h,(0,4))
display(h)
print("Urey Ratio convergence:\n")
display(Qrm_now_vec)

## --- Explore different Ur at high d (8.0)

# Get fresh parameters and set difficulty
p = Parameters()
p.Rc = (p.eta_L / (10.0^8))^(1/3)

nSteps = 10
UrMin = 0.1
UrMax = 0.5
Ur = reverse(linspace(UrMin,UrMax,nSteps))

cmap = cgrad((plasma[40:212]))
h = plot(TrelObs_time,TrelObs,yerror=2*TrelObs_sigma,seriestype=:scatter, markersize=2, color=:darkblue, markerstrokecolor=:auto, label="Data")
plottime = downsample(timevec/1000,50)
for i=1:nSteps
    (Tm, Tc, Qm, Qc, dp) = coolingmodel(p, Ur[i], Hmvec, dt, timevec)
    plotTrel = downsample(Tm .- Tm[1], 50)
    if i==1
        plot!(h,plottime,plotTrel,line_z=Ur[i],line=(cmap),label="Model") #color=cmap[i]
    else
        plot!(h,plottime,plotTrel,line_z=Ur[i],line=(cmap),label="") #color=cmap[i]
    end
end
plot!(h,legend=:topleft,xlabel="Age (Ma)",ylabel="\\Delta T (C)",framestyle=:box,grid=:off,fg_color_legend=:white)
ylims!(h,(-20,250))
xlims!(h,(0,4))
display(h)

savefig(h,"d=8variableUr.pdf")

## --- Explore different Ur at low d (5.5)

# Get fresh parameters and set difficulty
p = Parameters()
p.Rc = (p.eta_L / (10.0^5.5))^(1/3)

nSteps = 10
UrMin = 0.7555
UrMax = 0.7558

Ur = reverse(linspace(UrMin,UrMax,nSteps))
cmap = cgrad((plasma[40:212]))

h = plot(TrelObs_time,TrelObs,yerror=2*TrelObs_sigma,seriestype=:scatter, markersize=2, color=:darkblue, markerstrokecolor=:auto, label="Data")
plottime = downsample(timevec/1000,50)
for i=1:nSteps
    (Tm, Tc, Qm, Qc, dp) = coolingmodel(p, Ur[i], Hmvec, dt, timevec)
    plotTrel = downsample(Tm .- Tm[1], 50)
    if i==1
        plot!(h,plottime,plotTrel,line_z=Ur[i],line=(cmap),label="Model") #color=cmap[i]
    else
        plot!(h,plottime,plotTrel,line_z=Ur[i],line=(cmap),label="") #color=cmap[i]
    end
end
plot!(h,legend=:topleft,xlabel="Age (Ma)",ylabel="\\Delta T (C)",framestyle=:box,fg_color_legend=:white)
xlims!(h,(0,4))
ylims!(h,(-20,250))
display(h)
savefig(h,"d=5.5variableUr.pdf")


## --- Explore different Ur at low d (5.6)

# Get fresh parameters and set difficulty
p = Parameters()
p.Rc = (p.eta_L / (10.0^5.6))^(1/3)

nSteps = 12
UrMin = 0.7528
UrMax = 0.754

Ur = reverse(linspace(UrMin,UrMax,nSteps))
cmap = cgrad((plasma[40:212]))

h = plot(TrelObs_time,TrelObs,yerror=2*TrelObs_sigma,seriestype=:scatter, color=:darkblue, markersize=2, markerstrokecolor=:auto, label="Data")
plottime = downsample(timevec/1000,50)
for i=1:nSteps
    # Calculate the mantle temperature curve
    (Tm, Tc, Qm, Qc, dp) = coolingmodel(p, Ur[i], Hmvec, dt, timevec)

    # Plot the resutls
    plotTrel = downsample(Tm .- Tm[1], 50)
    if i==1
        plot!(h,plottime,plotTrel,line_z=Ur[i],line=(cmap),label="Model")
    else
        plot!(h,plottime,plotTrel,line_z=Ur[i],line=(cmap),label="")
    end
end
plot!(h,legend=:topleft,xlabel="Age (Ma)",ylabel="\\Delta T (C)",framestyle=:box,fg_color_legend=:white)
xlims!(h,(0,4))
ylims!(h,(-20,250))
display(h)
savefig(h,"d=5.6variableUr.pdf")

## --- Explore nl/Rc3 space - small range

nNewton = 20
dUr = 1E-9
n_difficulties = 100
minDiff = 4.5
maxDiff = 6.5
difficulies = linspace(minDiff,maxDiff,n_difficulties)
Ur_offset = 1E-3
p = Parameters() # Default parameters
descrip = "default"

# p.Q_addtl = 15E12 # Add 15 TW additional heat flux
# descrip = "Qaddtl15"

# p.Q_addtl = 3E12 # Add 10 TW additional heat flux
# p.Qc_now = 12E12 # Turn up present core heat flux from 6 TW
# p.Tc_now = 6400+273.15 # Turn up present core temperature from 4400 K
# p.Qm_now = 25E12 # Turn down present-day mantle heat flux from 35 TW.
# descrip = "multiple"

# p.Qc_now = 20E12 # Turn up present core heat flux from 6 TW
# p.Tc_now = 6400+273.15 # Turn up present core temperature from 4400 K
# descrip = "Qc20"

p.Qm_now = 13E12 # Turn down present-day mantle heat flux from 35 TW.
descrip = "Qm13"

Rc_t = (p.eta_L ./ (10.0.^difficulies)).^(1/3)
Ur_best = Array{Float64}(undef, size(Rc_t))
residual_best = Array{Float64}(undef, size(Rc_t))
residual_offset_mille = Array{Float64}(undef, size(Rc_t))
residual_offset_cent = Array{Float64}(undef, size(Rc_t))

cmap = cgrad(reverse(viridis[1:212]))
h1 = plot(TrelObs_time,TrelObs,yerror=2*TrelObs_sigma,seriestype=:scatter, color=:darkblue, markersize=2, markerstrokecolor=:auto, label="Data")
h2 = plot(xlabel="Age (Ma)",ylabel="Core temperature (K)",framestyle=:box)
h3 = plot(xlabel="Age (Ma)",ylabel="Core heat flux (W)",framestyle=:box)
h4 = plot(xlabel="Age (Ma)",ylabel="Mantle heat flux (W)",framestyle=:box)
plottime = downsample(timevec/1000,50)
@time for n=1:n_difficulties
    # Find best-fit Ur for this Rc_t[n]
    p.Rc = Rc_t[n]
    Ur_best[n] = coolingmodel_Newton_Ur(p,nNewton,Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma)

    # Keep track of results
    residual_best[n] = coolingmodel_resid(p,Ur_best[n],Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma)
    residual_offset_mille[n] = coolingmodel_resid(p,Ur_best[n]+0.001,Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma) / 2
    residual_offset_mille[n] += coolingmodel_resid(p,Ur_best[n]-0.001,Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma) / 2
    residual_offset_cent[n] = coolingmodel_resid(p,Ur_best[n]+0.01,Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma) / 2
    residual_offset_cent[n] += coolingmodel_resid(p,Ur_best[n]-0.01,Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma) / 2
    (Tm, Tc, Qm, Qc, dp) = coolingmodel(p,Ur_best[n],Hmvec,dt,timevec)

    # Plot results
    plotTrel = downsample(Tm .- Tm[1], 50)
    if n == n_difficulties
        plot!(h1,plottime,plotTrel,line_z=difficulies[n],line=(cmap),label="Model")
    else
        plot!(h1,plottime,plotTrel,line_z=difficulies[n],line=(cmap),label="")
    end
    plot!(h2,plottime,downsample(Tc,50),line_z=difficulies[n],line=(cmap),label="")
    plot!(h3,plottime,downsample(Qc,50),line_z=difficulies[n],line=(cmap),label="")
    plot!(h4,plottime,downsample(Qm,50),line_z=difficulies[n],line=(cmap),label="")
end
plot!(h1,legend=:topleft,xlabel="Age (Ma)",ylabel="\\Delta T (C)",framestyle=:box)
xlims!(h1,(0,4))
ylims!(h1,(-20,250))
display(h1)
savefig(h1,"$(descrip)_CoolingCurveVsDifficulty.pdf")
display(h2)
savefig(h2,"$(descrip)_CoolingCurveVsDifficulty_coretemp.pdf")
display(h3)
savefig(h3,"$(descrip)_CoolingCurveVsDifficulty_coreflux.pdf")
display(h4)
savefig(h4,"$(descrip)_CoolingCurveVsDifficulty_mantleflux.pdf")


# Two-axis plot for best-fit d and Ur
h = plot(xlabel="log_{10}(\\eta_L/Rc^3)",ylabel="Sum squared residual",xlims=(minDiff,maxDiff),legend=:topright,fg_color_legend=:white)
plot!(log10.(p.eta_L./Rc_t.^3),residual_best,label="Residual at best-fit Ur",grid=:off)
plot!(log10.(p.eta_L./Rc_t.^3),residual_offset_mille,label="with Ur offset 1 per mille",ylims=(0,1.5*maximum(residual_best)),color=:blue)
plot!(log10.(p.eta_L./Rc_t.^3),residual_offset_cent,label="with Ur offset 1 percent",ylims=(0,1.5*maximum(residual_best)),color=:darkblue)
plot!(twinx(),log10.(p.eta_L./Rc_t.^3),Ur_best,label="",ylims=(0,1),grid=:off,color=:darkred,ylabel="Best-fit Urey ratio",xlims=(minDiff,maxDiff))
display(h)
savefig(h,"$(descrip)_Ur=f(Difficulty).pdf")


## --- Explore nl/Rc3 space -- Full range

nNewton = 40
n_difficulties = 100
minDiff = 4
maxDiff = 9
difficulies = linspace(minDiff,maxDiff,n_difficulties)
Ur_offset = 0.001
p = Parameters()
descrip = "default"

p.Qc_now = 20E12 # Turn up present core heat flux from 6 TW
p.Tc_now = 6400+273.15 # Turn up present core temperature from 4400 K
descrip = "Qc20"

p.Qm_now = 13E12 # Turn down present-day mantle heat flux from 35 TW.
descrip = "Qm13"

Rc_t = (p.eta_L ./ (10 .^ difficulies)).^(1/3)
Ur_best = Array{Float64}(undef, size(Rc_t))
residual_best = Array{Float64}(undef, size(Rc_t))
residual_offset_mille = Array{Float64}(undef, size(Rc_t))
residual_offset_cent = Array{Float64}(undef, size(Rc_t))

cmap = cgrad(reverse(viridis[1:212]))
h = plot(TrelObs_time,TrelObs,yerror=2*TrelObs_sigma,seriestype=:scatter, color=:darkblue, markersize=2, markerstrokecolor=:auto, label="Data")
plottime = downsample(timevec/1000,50)
@time for n=1:n_difficulties
    # Find best-fit Ur for this Rc_t[n]
    p.Rc = Rc_t[n]
    Ur_best[n] = coolingmodel_Newton_Ur(p,nNewton,Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma)

    # Keep track of results
    residual_best[n] = coolingmodel_resid(p,Ur_best[n],Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma)
    residual_offset_mille[n] = coolingmodel_resid(p,Ur_best[n]+0.001,Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma) / 2
    residual_offset_mille[n] += coolingmodel_resid(p,Ur_best[n]-0.001,Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma) / 2
    residual_offset_cent[n] = coolingmodel_resid(p,Ur_best[n]+0.01,Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma) / 2
    residual_offset_cent[n] += coolingmodel_resid(p,Ur_best[n]-0.01,Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma) / 2
    (Tm, Tc, Qm, Qc, dp) = coolingmodel(p,Ur_best[n],Hmvec,dt,timevec)

    # Plot results
    plotTrel = downsample(Tm .- Tm[1], 50)
    if n == n_difficulties
        plot!(h,plottime,plotTrel,line_z=difficulies[n],line=(cmap),label="Model")
    elseif mod(n,2)==0
        plot!(h,plottime,plotTrel,line_z=difficulies[n],line=(cmap),label="")
    end
end
plot!(h,legend=:topleft,xlabel="Age (Ma)",ylabel="\\Delta T (C)",framestyle=:box)
xlims!(h,(0,4))
ylims!(h,(-20,250))
display(h)
savefig(h,"$(descrip)_CoolingCurveVsDifficultyHighRes.pdf")

# Two-axis plot for best-fit d and Ur
h = plot(xlabel="log_{10}(\\eta_L/Rc^3)",ylabel="Sum squared residual",xlims=(minDiff,maxDiff),legend=:topleft,fg_color_legend=:white)
plot!(log10.(p.eta_L./Rc_t.^3),residual_best,label="Residual at best-fit Ur",grid=:off)
plot!(log10.(p.eta_L./Rc_t.^3),residual_offset_mille,label="with Ur offset 1 per mille",ylims=(0,1.5*maximum(residual_best)),color=:blue)
plot!(log10.(p.eta_L./Rc_t.^3),residual_offset_cent,label="with Ur offset 1 percent",ylims=(0,1.5*maximum(residual_best)),color=:darkblue)
plot!(twinx(),log10.(p.eta_L./Rc_t.^3),Ur_best,label="",ylims=(0,1),grid=:off,color=:darkred,ylabel="Best-fit Urey ratio",xlims=(minDiff,maxDiff))
display(h)
savefig(h,"$(descrip)_Ur=f(Difficulty)HighRes.pdf")


## --- MCMC-Newton hybrid inversion

nSteps = 10^3
burnin = 5*10^2
d = 6 # Initial guess for difficulty
d_step_sigma = 0.5 # Standard deviation of Gaussian proposal distribution by which we adjust difficulty

p = Parameters() # Start with fresh parameters
p.Qm_now = 15E12 # Turn down present-day mantle heat flux from 35 TW.

# Define a function to conduct a hybrid MCMC-Newton inversion,
# adjusting d randomly, and finding the best-fit Ur and Tm path
function MCMC_Newton_SMC(nSteps, d, d_step_sigma, p)

    # Allocate arrays to record the stationary distributions
    d_dist = Array{Float64}(undef, nSteps)
    Ur_dist = Array{Float64}(undef, nSteps)
    Trel_dist = Array{Float64}(undef, nSteps,length(TrelObs_time))
    ll_dist = Array{Float64}(undef, nSteps)
    acceptancedist = fill(false,nSteps)

    # First markov chain step
    p.Rc = (p.eta_L/(10.0^d))^(1/3)
    Ur = coolingmodel_Newton_Ur(p,nNewton,Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma)
    Tm = coolingmodel(p, Ur, Hmvec, dt, timevec)[1]
    Trel = linterp1(timevec,Tm .- Tm[1], TrelObs_time*1000)
    ll = sum( -(Trel - TrelObs).^2 ./ (2*TrelObs_sigma.^2) - log.(sqrt.(2*pi*TrelObs_sigma)))

    d_dist[1] = d
    Ur_dist[1] = Ur
    ll_dist[1] = ll

    @showprogress for n=1:nSteps
       d_prop = d + randn()*d_step_sigma
       p.Rc = (p.eta_L / (10.0^d_prop))^(1/3)
       Ur_prop = coolingmodel_Newton_Ur(p,nNewton,Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma)
       Tm_prop = coolingmodel(p, Ur_prop, Hmvec, dt, timevec)[1]
       Trel_prop = linterp1(timevec, Tm_prop .- Tm_prop[1], TrelObs_time*1000)
       ll_prop = sum( -(Trel_prop - TrelObs).^2 ./ (2*TrelObs_sigma.^2) - log.(sqrt.(2*pi*TrelObs_sigma)))

       # Accept or reject proposal based on likelihood
       if rand() < exp(ll_prop-ll)
           d_step_sigma = 2.9*abs(d_prop-d)
           Trel = Trel_prop
           d = d_prop
           Ur = Ur_prop
           ll = ll_prop
           acceptancedist[n] = true
       end

       # Record results
       Trel_dist[n,:]=Trel
       d_dist[n] = d
       Ur_dist[n] = Ur
       ll_dist[n] = ll
    end
    return (d_dist, Ur_dist, Trel_dist, ll_dist, acceptancedist)
end

# Run the hybrid inversion
(d_dist, Ur_dist, Trel_dist, ll_dist, acceptancedist) = MCMC_Newton_SMC(nSteps, d, d_step_sigma, p)

Ur_mu = nanmean(Ur_dist[burnin:end])
Ur_sigma = nanstd(Ur_dist[burnin:end])

d_mu = nanmean(d_dist[burnin:end])
d_sigma = nanstd(d_dist[burnin:end])

p.Rc = (p.eta_L/(10.0^d_mu))^(1/3)
(Tm_mu, Tc_mu, Qm_mu, Qc_mu, dp_mu) = coolingmodel(p,Ur_mu,Hmvec,dt,timevec)

h = plot(TrelObs_time, TrelObs, yerror=2*TrelObs_sigma, seriestype=:scatter, color=:darkblue, markersize=2, markerstrokecolor=:auto, label="Data")
plot!(h,timevec/1000, Tm_mu .- Tm_mu[1], label="Model: Ur=$(round(Ur_mu,sigdigits=4)),d=$(round(d_mu,sigdigits=4))", legend=:topleft)
plot!(h,xlabel="Age (Ga)",ylabel="\\Delta T (C)")
display(h)
savefig(h,"MCMC-Newton Best Fit.pdf")

print(" Ur: $Ur_mu +/- $Ur_sigma\n d: $d_mu +/- $d_sigma\n\n")


## --- Full unconstrained MCMC inversion

function LL(x_prop,x,sigma)
    return sum( -(x_prop .- x).^2 ./ (2*sigma.^2) .- log.(sqrt.(2*pi*sigma)))
end

nSteps = 10^6
burnin = 5*10^4
jumpingsigmafactor = 2.718  # This can be adjusted to optimize acceptance likelihood

# Start with fresh parameters
p = Parameters()

# Initial guesses
# d = 6
# p.Rc = (p.eta_L/(10.0^d))^(1/3)
# p.Qm_now = 25E12
# p.Qc_now = 13E12
# p.Tc_now = 6320
# p.Q_addtl = 1E12

# Better initial guesses
d = 8
p.Rc = (p.eta_L/(10.0^d))^(1/3)
p.Qm_now = 35E12
p.Qc_now = 10E12
p.Tc_now = 4400+273
p.Q_addtl = 0E12
descrip = "Qm35-10_Qc10-10"

# Define a function to conduct a full MCMC inversion
function MCMC_SMC(nSteps, d, jumpingsigmafactor, p)

    # Inital constraints
    Qm_now = p.Qm_now
    Qc_now = p.Qc_now
    Tc_now = p.Tc_now
    Q_addtl = p.Q_addtl

    # Allocate arrays to record the stationary distributions
    acceptancedist = fill(false,nSteps)
    ll_dist = Array{Float64}(undef, nSteps)
    d_dist = Array{Float64}(undef, nSteps)
    Ur_dist = Array{Float64}(undef, nSteps)
    Qm_now_dist = Array{Float64}(undef, nSteps)
    Qc_now_dist = Array{Float64}(undef, nSteps)
    Tc_now_dist = Array{Float64}(undef, nSteps)
    Q_addtl_dist = Array{Float64}(undef, nSteps)
    Trel_dist = Array{Float64}(undef, nSteps,length(TrelObs_time))

    # Proposal distributions
    d_step_sigma = 0.1
    Ur_step_sigma = 0.005
    Qm_step_sigma = 0.1E12
    Qc_step_sigma = 0.1E12
    Tc_step_sigma = 100
    Q_addtl_step_sigma = 0.1E12

    # First markov chain step
    Ur = coolingmodel_Newton_Ur(p,20,Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma)
    # Ur = 0.3289
    Tm = coolingmodel(p, Ur, Hmvec, dt, timevec)[1]
    Trel = linterp1(timevec,Tm .- Tm[1], TrelObs_time*1000)
    h = plot(TrelObs_time,TrelObs,yerror=2*TrelObs_sigma,seriestype=:scatter, color=:darkblue, markersize=2, markerstrokecolor=:auto, label="Data")
    plot!(h,timevec/1000,Tm .- Tm[1], label="Model: Ur=$(round(Ur,sigdigits=3)),d=$(round(d,sigdigits=3))")
    plot!(h,xlabel="Age (Ga)",ylabel="\\Delta T (C)",grid=:off,legend=:topleft,fg_color_legend=:white)
    summarystring = " d:\t\t$(round(d,sigdigits=2))\n" *
        " Ur: $(round(Ur,sigdigits=3))\n\n" *
        " Qc: $(round(Qc_now/1E12,sigdigits=1))\n" *
        " Qm: $(round(Qm_now/1E12,sigdigits=1))\n" *
        " Tc: $(round(Int,Tc_now))\n"
    plot!(h,1,1,seriestype=:scatter, ann=(0.16,130,text(summarystring,8,"Helvetica",:left)),label="")
    display(h)
    savefig(h,"MCMC_firstproposal.pdf")

    ll = LL(Trel,TrelObs,TrelObs_sigma*2) + LL(Ur,0.3289,0.01) + LL(d,5,5)
        LL(Qm_now/1E12,35,10) + LL(Qc_now/1E12,10,10) + LL(Tc_now,4400+273.15,500)

    # Accept and record first proposal
    ll_dist[1] = ll
    d_dist[1] = d
    Ur_dist[1] = Ur
    Qm_now_dist[1] = Qm_now
    Qc_now_dist[1] = Qc_now
    Tc_now_dist[1] = Tc_now
    Q_addtl_dist[1] = Q_addtl
    Trel_dist[1,:] = Trel

    # Run the Markov chain / Metropolis walker
    @showprogress for n=1:nSteps

        # Randomly choose which parameter to adjust. The main benefit of his
        # form of "Metropolis-in-Gibbs" / "Alternate conditional sampling"
        # is that it allows us to optimize the proposal (jumping) distributions
        # for each parameter individually
        case = ceil(rand()*5)

        # Propose new parameters
        Ur_prop = case==1 ? abs(Ur + randn()*Ur_step_sigma) : Ur
        d_prop  = case==2 ? abs(d + randn()*d_step_sigma) : d
        p.Rc = abs(p.eta_L / (10.0^d_prop))^(1/3)
        p.Qm_now = case==3 ? abs(Qm_now + randn()*Qm_step_sigma) : Qm_now
        p.Qc_now = case==4 ? abs(Qc_now + randn()*Qc_step_sigma) : Qc_now
        p.Tc_now = case==5 ? abs(Tc_now + randn()*Tc_step_sigma) : Tc_now
        p.Q_addtl = case==6 ? abs(Q_addtl + randn()*Q_addtl_step_sigma) : Q_addtl

        # If proposal is so bad that we hit a math error, let ll be NaN to
        # auto-reject the proposal
        Trel_prop = NaN
        ll_prop = NaN
        try
            # Calculate proposed solution
            Tm_prop = coolingmodel(p, Ur_prop, Hmvec, dt, timevec)[1]
            Trel_prop = linterp1(timevec, Tm_prop .- Tm_prop[1], TrelObs_time*1000)

            # Calculate log likelihood of new proposal
            ll_prop = LL(Trel_prop,TrelObs,TrelObs_sigma*2) + LL(Ur_prop,0.3289,0.01) + LL(d_prop,5,5) +
               LL(p.Qm_now/1E12,35,10) + LL(p.Qc_now/1E12,10,10) + LL(p.Tc_now,4400+273.15,500)
        catch
            @info("Model failed to converge, trying again")
        end

        # Accept or reject proposal based on likelihood
        if rand() < exp(ll_prop-ll)
            # Update proposal (jumping) distributions
            Ur_step_sigma = case==1 ? jumpingsigmafactor*abs(Ur_prop - Ur) : Ur_step_sigma
            d_step_sigma  = case==2 ? jumpingsigmafactor*abs(d_prop-d) : d_step_sigma
            Qm_step_sigma = case==3 ? jumpingsigmafactor*abs(p.Qm_now - Qm_now) : Qm_step_sigma
            Qc_step_sigma = case==4 ? jumpingsigmafactor*abs(p.Qc_now - Qc_now) : Qc_step_sigma
            Tc_step_sigma = case==5 ? jumpingsigmafactor*abs(p.Tc_now - Tc_now) : Tc_step_sigma
            Q_addtl_step_sigma = case==6 ? jumpingsigmafactor*abs(p.Q_addtl - Q_addtl) : Q_addtl_step_sigma

            # Accept results
            ll = ll_prop
            d = d_prop
            Ur = Ur_prop
            Qm_now  = p.Qm_now
            Qc_now  = p.Qc_now
            Tc_now  = p.Tc_now
            Q_addtl = p.Q_addtl
            Trel = Trel_prop

            # Record accptance
            acceptancedist[n] = true
        end

        # Record results
        ll_dist[n] = ll
        d_dist[n] = d
        Ur_dist[n] = Ur
        Qm_now_dist[n] = Qm_now
        Qc_now_dist[n] = Qc_now
        Tc_now_dist[n] = Tc_now
        Q_addtl_dist[n] = Q_addtl
        Trel_dist[n,:] = Trel
    end

    return (acceptancedist, ll_dist, d_dist, Ur_dist, Qm_now_dist, Qc_now_dist, Tc_now_dist, Q_addtl_dist, Trel_dist)
end

# Run the full MCMC inversion
(acceptancedist, ll_dist, d_dist, Ur_dist, Qm_now_dist, Qc_now_dist, Tc_now_dist, Q_addtl_dist, Trel_dist) =  MCMC_SMC(nSteps, d, jumpingsigmafactor, p)

# Summary statistics
d_mu = nanmean(d_dist[burnin:end])
d_sigma = nanstd(d_dist[burnin:end])
Ur_mu = nanmean(Ur_dist[burnin:end])
Ur_sigma = nanstd(Ur_dist[burnin:end])
Qm_now_mu = mean(Qm_now_dist[burnin:end])
Qm_now_sigma = std(Qm_now_dist[burnin:end])
Qc_now_mu = mean(Qc_now_dist[burnin:end])
Qc_now_sigma = std(Qc_now_dist[burnin:end])
Tc_now_mu = mean(Tc_now_dist[burnin:end])
Tc_now_sigma = std(Tc_now_dist[burnin:end])
Q_addtl_mu = mean(Q_addtl_dist[burnin:end])
Q_addtl_sigma = std(Q_addtl_dist[burnin:end])

# Run model with average burnt-in parameters
p.Rc = (p.eta_L/(10.0^d_mu))^(1/3)
p.Qm_now = Qm_now_mu
p.Qc_now = Qc_now_mu
p.Tc_now = Tc_now_mu
p.Q_addtl = Q_addtl_mu
(Tm_mu, Tc_mu, Qm_mu, Qc_mu, dp_mu) = coolingmodel(p,Ur_mu,Hmvec,dt,timevec)

# Plot results
lines_to_plot = 500
cmap = cgrad(reverse(viridis[128:end]))
dsfactor = round(Int, (nSteps-burnin)/lines_to_plot)
zdata = downsample(d_dist[burnin:end], dsfactor)'
Trel_dist_plot = downsample(Trel_dist[burnin:end,:], dsfactor,1)'
h = plot(TrelObs_time, Trel_dist_plot, line_z=zdata, line=(cmap), linealpha=0.3, label="")
plot!(h, TrelObs_time, TrelObs, yerror=2*TrelObs_sigma, seriestype=:scatter, color=:darkblue, markersize=2, markerstrokecolor=:auto, label="Data")
plot!(h, timevec/1000, Tm_mu .- Tm_mu[1], color=:black, label="Model")
plot!(h, xlabel="Age (Ma)", ylabel="\\Delta T (C)" ,xlims=(0,4), framestyle=:box, grid=:off, legend=:topleft, fg_color_legend=:white)
summarystring = " d:\t\t$(round(d_mu,sigdigits=2))\n" *
    " Ur: $(round(Ur_mu,sigdigits=3))\n\n" *
    " Qc: $(round(Qc_now_mu/1E12,sigdigits=1))\n" *
    " Qm: $(round(Qm_now_mu/1E12,sigdigits=1))\n" *
    " Tc: $(round(Int,Tc_now_mu))\n"
plot!(h,1,1,seriestype=:scatter, ann=(0.16,130,text(summarystring,8,"Helvetica",:left)),label="")
display(h)
savefig(h,"MCMC_bestfit_$descrip.pdf")

# Print results
resultstring = " Ur: $(round(Ur_mu,sigdigits=6)) +/- $(round(Ur_sigma,sigdigits=6))\n" *
 " d: $(round(d_mu,sigdigits=4)) +/- $(round(d_sigma,sigdigits=4))\n" *
 " Qm_now: $(round(Qm_now_mu,sigdigits=4)) +/- $(round(Qm_now_sigma,sigdigits=4))\n" *
 " Qc_now: $(round(Qc_now_mu,sigdigits=4)) +/- $(round(Qc_now_sigma,sigdigits=4))\n" *
 " Tc_now: $(round(Tc_now_mu,sigdigits=4)) +/- $(round(Tc_now_sigma,sigdigits=4))\n" *
 " Q_addtl: $(round(Q_addtl_mu,sigdigits=4)) +/- $(round(Q_addtl_sigma,sigdigits=4))\n"
print(resultstring * " acceptance probability: $(round(mean(acceptancedist[burnin:end]),sigdigits=2))\n")

## --- End of File