    1+1
## ---

    using StatGeochem
    include("Utilities.jl")
    using Plots; gr()
    using ProgressMeter: @showprogress

## --- Define physical parameters and set up model

    # Time step (Myr)
    dt = 0.1
    timevec = collect(0:dt:4000)

    # Calculate past mantle heat production
    p = Parameters()
    Hmvec = mantle_heat_production(p,mcdbse,timevec)

## --- Run Newton's method

    nNewton = 20 # Number of Newton's method iterations
    dUr = 1E-9 # Urey ratio step size
    p.Rc = 750E3 # Plate bending radius (km)

    Ur = Array{Float64}(nNewton)
    Ur[1] = 0.5 # Initial guess for Urey ratio
    for i = 1:nNewton-1
        ll = coolingmodel_LL(p,Ur[i],Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma)
        llp = coolingmodel_LL(p,Ur[i]+dUr,Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma)
        dll_dUr = (llp-ll)/dUr
        Ur[i+1] = Ur[i] - ll/dll_dUr/2
    end

    (Tm, Tc, Qm, Qc, dp) = coolingmodel(p, Ur[nNewton], Hmvec, dt, timevec)
    h = plot(TrelObs_time, TrelObs, yerror=2*TrelObs_sigma, seriestype=:scatter, markersize=2, color=:darkblue, markerstrokecolor=:auto, label="Data")
    plot!(h, timevec/1000,Tm-Tm[1], label="Model", xlabel="Age (Ma)", ylabel="\\Delta T (C)", legend=:topleft, framestyle=:box)
    xlims!(h,(0,4))
    display(h)
    print("Urey Ratio convergence:\n")
    display(Ur)

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
        (Tm, Tc, Qm, Qc, dp) = coolingmodel(p,Ur[i],Hmvec,dt,timevec)
        plotTrel = downsample(Tm-Tm[1],50)
        if i==1
            plot!(h,plottime,plotTrel,line_z=Ur[i],line=(cmap),label="Model") #color=cmap[i]
        else
            plot!(h,plottime,plotTrel,line_z=Ur[i],line=(cmap),label="") #color=cmap[i]
        end
    end
    plot!(h,legend=:topleft,xlabel="Age (Ma)",ylabel="\\Delta T (C)",framestyle=:box)
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
        (Tm, Tc, Qm, Qc, dp) = coolingmodel(p,Ur[i],Hmvec,dt,timevec)
        plotTrel = downsample(Tm-Tm[1],50)
        if i==1
            plot!(h,plottime,plotTrel,line_z=Ur[i],line=(cmap),label="Model") #color=cmap[i]
        else
            plot!(h,plottime,plotTrel,line_z=Ur[i],line=(cmap),label="") #color=cmap[i]
        end
    end
    plot!(h,legend=:topleft,xlabel="Age (Ma)",ylabel="\\Delta T (C)",framestyle=:box)
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
        (Tm, Tc, Qm, Qc, dp) = coolingmodel(p,Ur[i],Hmvec,dt,timevec)

        # Plot the resutls
        plotTrel = downsample(Tm-Tm[1],50)
        if i==1
            plot!(h,plottime,plotTrel,line_z=Ur[i],line=(cmap),label="Model")
        else
            plot!(h,plottime,plotTrel,line_z=Ur[i],line=(cmap),label="")
        end
    end
    plot!(h,legend=:topleft,xlabel="Age (Ma)",ylabel="\\Delta T (C)",framestyle=:box)
    xlims!(h,(0,4))
    ylims!(h,(-20,250))
    display(h)
    savefig(h,"d=5.6variableUr.pdf")

## --- Explore nl/Rc3 space - small range

    nNewton = 20
    dUr = 1E-9
    n_difficulties = 20
    minDiff = 4.5
    maxDiff = 6.5
    difficulies = linspace(minDiff,maxDiff,n_difficulties)
    Ur_offset = 1E-3
    p = Parameters()
    # p.Q_addtl = 15E12 # Add 15 TW additional heat flux

    # p.Q_addtl = 3E12 # Add 10 TW additional heat flux
    # p.Qc_now = 12E12 # Turn up present core heat flux from 6 TW
    # p.Tc_now = 6400+273.15 # Turn up present core temperature from 4400 K
    # p.Qm_now = 25E12 # Turn down present-day mantle heat flux from 35 TW.

    # p.Qc_now = 20E12 # Turn up present core heat flux from 6 TW
    # p.Tc_now = 6400+273.15 # Turn up present core temperature from 4400 K

    # p.Qm_now = 13E12 # Turn down present-day mantle heat flux from 35 TW.


    Rc_t = (p.eta_L ./ (10.0.^difficulies)).^(1/3)
    Ur_best = Array{Float64}(size(Rc_t))
    residual_best = Array{Float64}(size(Rc_t))
    residual_offset_mille = Array{Float64}(size(Rc_t))
    residual_offset_cent = Array{Float64}(size(Rc_t))

    cmap = cgrad(flipdim(viridis[1:212],1))
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
        plotTrel = downsample(Tm-Tm[1],50)
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
    savefig(h1,"CoolingCurveVsDifficulty.pdf")

    display(h2)
    display(h3)
    display(h4)

    # Two-axis plot for best-fit d and Ur
    h = plot(xlabel="log_{10}(\\eta_L/Rc^3)",ylabel="Sum squared residual",xlims=(minDiff,maxDiff),legend=:topleft,fg_color_legend=:white)
    plot!(log10.(p.eta_L./Rc_t.^3),residual_best,label="Residual at best-fit Ur",grid=:off)
    plot!(log10.(p.eta_L./Rc_t.^3),residual_offset_mille,label="with Ur offset 1 per mille",ylims=(0,1.5*maximum(residual_best)),grid=:off,color=:darkblue)
    plot!(log10.(p.eta_L./Rc_t.^3),residual_offset_cent,label="with Ur offset 1 percent",ylims=(0,1.5*maximum(residual_best)),grid=:off,color=:black)
    plot!(twinx(),log10.(p.eta_L./Rc_t.^3),Ur_best,label="",ylims=(0,1),grid=:off,color=:darkred,ylabel="Best-fit Urey ratio",xlims=(minDiff,maxDiff))
    display(h)


## --- Explore nl/Rc3 space -- Full range

    p.eta_L = 1E23 # Plate viscsity (Pa s)

    nNewton = 20
    n_difficulties = 100
    minDiff = 4
    maxDiff = 9
    difficulies = linspace(minDiff,maxDiff,n_difficulties)
    Ur_offset = 0.001
    p = Parameters()

    Rc_t = (p.eta_L ./ (10 .^ difficulies)).^(1/3)
    Ur_best = Array{Float64}(size(Rc_t))
    residual_best = Array{Float64}(size(Rc_t))
    residual_offset_mille = Array{Float64}(size(Rc_t))
    residual_offset_cent = Array{Float64}(size(Rc_t))

    cmap = cgrad(flipdim(viridis[1:212],1))
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
        plotTrel = downsample(Tm-Tm[1],50)
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
    savefig(h,"CoolingCurveVsDifficultyHighRes.pdf")

    # Two-axis plot for best-fit d and Ur
    h = plot(xlabel="log_{10}(\\eta_L/Rc^3)",ylabel="Sum squared residual",xlims=(minDiff,maxDiff),legend=:topleft,fg_color_legend=:white)
    plot!(log10.(p.eta_L./Rc_t.^3),residual_best,label="Residual at best-fit Ur",grid=:off)
    plot!(log10.(p.eta_L./Rc_t.^3),residual_offset_mille,label="with Ur offset 1 per mille",ylims=(0,1.5*maximum(residual_best)),grid=:off,color=:darkblue)
    plot!(log10.(p.eta_L./Rc_t.^3),residual_offset_cent,label="with Ur offset 1 percent",ylims=(0,1.5*maximum(residual_best)),grid=:off,color=:black)
    plot!(twinx(),log10.(p.eta_L./Rc_t.^3),Ur_best,label="",ylims=(0,1),grid=:off,color=:darkred,ylabel="Best-fit Urey ratio",xlims=(minDiff,maxDiff))
    display(h)
    savefig(h,"Ur=f(Difficulty).pdf")


## --- MCMC-Newton hybrid inversion

    nSteps = 10^3
    burnin = 5*10^2

    d_dist = Array{Float64}(nSteps)
    Ur_dist = Array{Float64}(nSteps)
    Trel_dist = Array{Float64}(nSteps,length(TrelObs_time))
    ll_dist = Array{Float64}(nSteps)
    acceptancedist = fill(false,nSteps)

    p = Parameters() # Start with fresh parameters
    d = 6 # Initial guess for difficulty
    d_step_sigma = 0.5

    p.Qm_now = 15E12 # Turn down present-day mantle heat flux from 35 TW.

    # First markov chain step
    p.Rc = (p.eta_L/(10.0^d))^(1/3)
    Ur = coolingmodel_Newton_Ur(p,nNewton,Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma)
    Tm = coolingmodel(p, Ur, Hmvec, dt, timevec)[1]
    Trel = linterp1(timevec,Tm-Tm[1],TrelObs_time*1000)
    ll = sum( -(Trel - TrelObs).^2 ./ (2*TrelObs_sigma.^2) - log.(sqrt.(2*pi*TrelObs_sigma)))

    d_dist[1] = d
    Ur_dist[1] = Ur
    ll_dist[1] = ll

    @showprogress for n=1:nSteps
       d_prop = d + randn()*d_step_sigma
       p.Rc = (p.eta_L / (10.0^d_prop))^(1/3)
       Ur_prop = coolingmodel_Newton_Ur(p,nNewton,Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma)
       Tm_prop = coolingmodel(p, Ur_prop, Hmvec, dt, timevec)[1]
       Trel_prop = linterp1(timevec,Tm_prop-Tm_prop[1],TrelObs_time*1000)
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

    Ur_mu = nanmean(Ur_dist[burnin:end])
    Ur_sigma = nanstd(Ur_dist[burnin:end])

    d_mu = nanmean(d_dist[burnin:end])
    d_sigma = nanstd(d_dist[burnin:end])

    p.Rc = (p.eta_L/(10.0^d_mu))^(1/3)
    (Tm_mu, Tc_mu, Qm_mu, Qc_mu, dp_mu) = coolingmodel(p,Ur_mu,Hmvec,dt,timevec)

    h = plot(TrelObs_time,TrelObs,yerror=2*TrelObs_sigma,seriestype=:scatter, color=:darkblue, markersize=2, markerstrokecolor=:auto, label="Data")
    plot!(h,timevec/1000,Tm_mu-Tm_mu[1],label="Model: Ur=$(round(Ur_mu,4)),d=$(round(d_mu,4))",legend=:topleft)
    plot!(h,xlabel="Age (Ga)",ylabel="\\Delta T (C)")
    display(h)
    savefig(h,"MCMC-Newton Best Fit.pdf")

    print(" Ur: $Ur_mu +/- $Ur_sigma\n d: $d_mu +/- $d_sigma\n\n")


## --- Full unconstrained MCMC inversion

    function LL(x_prop,x,sigma)
        return sum( -(x_prop - x).^2 ./ (2*sigma.^2) - log.(sqrt.(2*pi*sigma)))
    end

    nSteps = 10^4
    burnin = 5*10^3
    jumpingsigmafactor = 2.718

    acceptancedist = fill(false,nSteps)
    ll_dist = Array{Float64}(nSteps)
    d_dist = Array{Float64}(nSteps)
    Ur_dist = Array{Float64}(nSteps)
    Qm_now_dist = Array{Float64}(nSteps)
    Qc_now_dist = Array{Float64}(nSteps)
    Tc_now_dist = Array{Float64}(nSteps)
    Q_addtl_dist = Array{Float64}(nSteps)
    Trel_dist = Array{Float64}(nSteps,length(TrelObs_time))

    # Initial guesses
    # d = 6
    # Qm_now = 25E12
    # Qc_now = 13E12
    # Tc_now = 6320
    # Q_addtl = 1E12
    d = 8
    Qm_now = 35E12
    Qc_now = 6E12
    Tc_now = 4400+273
    Q_addtl = 0E12

    # Start with fresh parameters
    p = Parameters()
    p.Rc = (p.eta_L/(10.0^d))^(1/3)
    p.Qm_now = Qm_now
    p.Qc_now = Qc_now
    p.Tc_now = Tc_now
    p.Q_addtl = Q_addtl

    # Proposal distributions
    d_step_sigma = 0.1
    Ur_step_sigma = 0.005
    Qm_step_sigma = 0.1E12
    Qc_step_sigma = 0.1E12
    Tc_step_sigma = 100
    Q_addtl_step_sigma = 0.1E12

    # First markov chain step
    Ur = coolingmodel_Newton_Ur(p,20,Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma)
    # Ur = 0.31
    Tm = coolingmodel(p, Ur, Hmvec, dt, timevec)[1]
    Trel = linterp1(timevec,Tm-Tm[1],TrelObs_time*1000)
    h = plot(TrelObs_time,TrelObs,yerror=2*TrelObs_sigma,seriestype=:scatter, color=:darkblue, markersize=2, markerstrokecolor=:auto, label="Data")
    plot!(h,timevec/1000,Tm-Tm[1],label="Model: Ur=$(round(Ur,4)),d=$(round(d,4))",legend=:topleft)
    plot!(h,xlabel="Age (Ga)",ylabel="\\Delta T (C)")
    display(h)

    ll = LL(Trel,TrelObs,TrelObs_sigma) + LL(Ur,0.31,0.01) + LL(d,5,5)
        LL(Qm/1E12,35,2) + LL(Qc/1E12,10,5) + LL(Tc,4400+273.15,500)

    Trel_dist[1,:]=Trel
    d_dist[1] = d
    Ur_dist[1] = Ur
    ll_dist[1] = ll

    @showprogress for n=1:nSteps
        # Propose new parameters
        case = ceil(rand()*5) # Which parameter to adjust

        Ur_prop = case==1 ? abs(Ur + randn()*Ur_step_sigma) : Ur
        d_prop  = case==2 ? abs(d + randn()*d_step_sigma) : d
        p.Rc = abs(p.eta_L / (10.0^d_prop))^(1/3)
        p.Qm_now = case==3 ? abs(Qm_now + randn()*Qm_step_sigma) : Qm_now
        p.Qc_now = case==4 ? abs(Qc_now + randn()*Qc_step_sigma) : Qc_now
        p.Tc_now = case==5 ? abs(Tc_now + randn()*Tc_step_sigma) : Tc_now
        p.Q_addtl = case==6 ? abs(Q_addtl + randn()*Q_addtl_step_sigma) : Q_addtl

        # If proposal is so bad we hit a math error, auto-reject the proposal
        Trel_prop = NaN
        ll_prop = NaN
        try
            # Calculate proposed solution
            Tm_prop = coolingmodel(p, Ur_prop, Hmvec, dt, timevec)[1]
            Trel_prop = linterp1(timevec,Tm_prop-Tm_prop[1],TrelObs_time*1000)

            # Calculate log likelihood of new proposal
            ll_prop = LL(Trel_prop,TrelObs,TrelObs_sigma) + LL(Ur_prop,0.31,0.01) + LL(d_prop,5,5) +
               LL(p.Qm_now/1E12,35,2) + LL(p.Qc_now/1E12,10,5) + LL(p.Tc_now,4400+273.15,500)
        catch
            warn("Model failed to converge!")
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

    Ur_mu = nanmean(Ur_dist[burnin:end])
    Ur_sigma = nanstd(Ur_dist[burnin:end])

    d_mu = nanmean(d_dist[burnin:end])
    d_sigma = nanstd(d_dist[burnin:end])

    p.Rc = (p.eta_L/(10.0^d_mu))^(1/3)
    (Tm_mu, Tc_mu, Qm_mu, Qc_mu, dp_mu) = coolingmodel(p,Ur_mu,Hmvec,dt,timevec)

    h = plot(TrelObs_time,TrelObs,yerror=2*TrelObs_sigma,seriestype=:scatter, color=:darkblue, markersize=2, markerstrokecolor=:auto, label="Data")
    plot!(h,timevec/1000,Tm_mu-Tm_mu[1],label="Model: Ur=$(round(Ur_mu,4)),d=$(round(d_mu,4))",legend=:topleft, fg_color_legend=:white, framestyle=:box)
    plot!(h,xlabel="Age (Ga)",ylabel="\\Delta T (C)")
    display(h)
    savefig(h,"MCMC-Newton Best Fit.pdf")

    print(" Ur: $(signif(Ur_mu,6)) +/- $(signif(Ur_sigma,6))\n" *
     " d: $(signif(d_mu,4)) +/- $(signif(d_sigma,4))\n" *
     " Qm_now: $(signif(mean(Qm_now_dist[burnin:end]),4)) +/- $(signif(std(Qm_now_dist[burnin:end]),4))\n" *
     " Qc_now: $(signif(mean(Qc_now_dist[burnin:end]),4)) +/- $(signif(std(Qc_now_dist[burnin:end]),4))\n" *
     " Tc_now: $(signif(mean(Tc_now_dist[burnin:end]),4)) +/- $(signif(std(Tc_now_dist[burnin:end]),4))\n" *
     " Q_addtl: $(signif(mean(Q_addtl_dist[burnin:end]),4)) +/- $(signif(std(Q_addtl_dist[burnin:end]),4))\n" *
     " acceptance probability: $(signif(mean(acceptancedist[burnin:end]),2))\n\n")

## ---

    lines_to_plot = 500
    plot(TrelObs_time,Trel_dist[end-lines_to_plot:end,:]',line_z=d_dist[end-lines_to_plot:end]',line=(:viridis_r),linealpha=0.5,label="")
    plot!(TrelObs_time,TrelObs,yerror=2*TrelObs_sigma,seriestype=:scatter, color=:darkblue, markersize=2, markerstrokecolor=:auto, label="")
    plot!(xlabel="Age (Ma)", ylabel="\\Delta T (C)")

## --- End of File
