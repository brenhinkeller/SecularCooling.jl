    1+1
## ---

    using StatGeochem
    include("Utilities.jl")
    using Plots; gr()
    using LaTeXStrings

## --- Define physical parameters and set up model

    p = Parameters()
    p.eta_L = 1E23; # Plate viscsity (Pa s)

    dt = 0.1; # Time step (Myr)
    timevec = collect(0:dt:4000);

    # Calculate past mantle heat production
    Hmvec = mantle_heat_production(p,mcdbse,timevec)

## --- Run Newton's method

    Rc = 750E3 # Plate bending radius (km)
    dUr = 1E-9 # Urey ratio step size

    nSteps = 20
    Ur = Array{Float64}(nSteps)
    Ur[1] = 0.5
    for i = 1:nSteps-1
        ll = coolingmodel_LL(p,Rc,Ur[i],Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma);
        llp = coolingmodel_LL(p,Rc,Ur[i]+dUr,Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma);
        dll_dUr = (llp-ll)/dUr;
        Ur[i+1] = Ur[i] - ll/dll_dUr/2;
    end

    Tm = coolingmodel(p, Rc, Ur[nSteps], Hmvec, dt, timevec)
    h = plot(TrelObs_time, TrelObs, yerror=2*TrelObs_sigma, seriestype=:scatter, markersize=2, color=:darkblue, markerstrokecolor=:auto, label="Data")
    plot!(h, timevec/1000,Tm-Tm[1], label="Model", xlabel="Age (Ma)", ylabel="\\Delta T (C)", legend=:topleft)
    xlims!(h,(0,4))
    display(h)
    print("Urey Ratio convergence:\n")
    display(Ur)

## --- Explore different Ur at high d (8.0)

    Rc = (p.eta_L / (10.^8))^(1/3)

    nSteps = 10
    UrMin = 0.1
    UrMax = 0.5
    Ur = reverse(linspace(UrMin,UrMax,nSteps))

    cmap = cgrad((plasma[40:212]))
    h = plot(TrelObs_time,TrelObs,yerror=2*TrelObs_sigma,seriestype=:scatter, markersize=2, color=:darkblue, markerstrokecolor=:auto, label="Data")
    plottime = downsample(timevec/1000,50)
    for i=1:nSteps
        Tm = coolingmodel(p,Rc,Ur[i],Hmvec,dt,timevec);
        plotTrel = downsample(Tm-Tm[1],50)
        if i==1
            plot!(h,plottime,plotTrel,line_z=Ur[i],line=(cmap),label="Model") #color=cmap[i]
        else
            plot!(h,plottime,plotTrel,line_z=Ur[i],line=(cmap),label="") #color=cmap[i]
        end
    end
    plot!(h,legend=:topleft,xlabel="Age (Ma)",ylabel="\\Delta T (C)")
    ylims!(h,(-20, 300))
    xlims!(h,(0,4))
    display(h)

    savefig(h,"d=8variableUr.pdf")

## --- Explore different Ur at low d (5.5)

    nSteps = 10;
    Rc = (p.eta_L / (10^5.5))^(1/3);
    UrMin = 0.7555
    UrMax = 0.7558

    Ur = reverse(linspace(UrMin,UrMax,nSteps))
    cmap = cgrad((plasma[40:212]))

    h = plot(TrelObs_time,TrelObs,yerror=2*TrelObs_sigma,seriestype=:scatter, markersize=2, color=:darkblue, markerstrokecolor=:auto, label="Data")
    plottime = downsample(timevec/1000,50)
    for i=1:nSteps
        Tm = coolingmodel(p,Rc,Ur[i],Hmvec,dt,timevec);
        plotTrel = downsample(Tm-Tm[1],50)
        if i==1
            plot!(h,plottime,plotTrel,line_z=Ur[i],line=(cmap),label="Model") #color=cmap[i]
        else
            plot!(h,plottime,plotTrel,line_z=Ur[i],line=(cmap),label="") #color=cmap[i]
        end
    end
    plot!(h,legend=:topleft,xlabel="Age (Ma)",ylabel="\\Delta T (C)")
    xlims!(h,(0,4))
    ylims!(h,(-20,250))
    display(h)
    savefig(h,"d=5.5variableUr.pdf")


## --- Explore different Ur at low d (5.6)

    nSteps = 12;
    Rc = (p.eta_L / (10^5.6))^(1/3)
    UrMin = 0.7528
    UrMax = 0.754

    Ur = reverse(linspace(UrMin,UrMax,nSteps))
    cmap = cgrad((plasma[40:212]))

    h = plot(TrelObs_time,TrelObs,yerror=2*TrelObs_sigma,seriestype=:scatter, color=:darkblue, markersize=2, markerstrokecolor=:auto, label="Data")
    plottime = downsample(timevec/1000,50)
    for i=1:nSteps
        # Calculate the mantle temperature curve
        Tm = coolingmodel(p,Rc,Ur[i],Hmvec,dt,timevec)

        # Plot the resutls
        plotTrel = downsample(Tm-Tm[1],50)
        if i==1
            plot!(h,plottime,plotTrel,line_z=Ur[i],line=(cmap),label="Model")
        else
            plot!(h,plottime,plotTrel,line_z=Ur[i],line=(cmap),label="")
        end
    end
    plot!(h,legend=:topleft,xlabel="Age (Ma)",ylabel="\\Delta T (C)")
    xlims!(h,(0,4))
    ylims!(h,(-20,250))
    display(h)
    savefig(h,"d=5.6variableUr.pdf")

## --- Explore nl/Rc3 space

    p.eta_L = 1E23 # Plate viscsity (Pa s)
    nSteps = 40

    n_difficulties = 50
    minDiff = 4
    maxDiff = 9
    difficulies = linspace(minDiff,maxDiff,n_difficulties)

    Rc_t = (p.eta_L ./ (10 .^ difficulies)).^(1/3)
    Ur_b = Array{Float64}(size(Rc_t))
    R_b = Array{Float64}(size(Rc_t))

    cmap = cgrad(flipdim(viridis[1:212],1))
    h = plot(TrelObs_time,TrelObs,yerror=2*TrelObs_sigma,seriestype=:scatter, color=:darkblue, markersize=2, markerstrokecolor=:auto, label="Data")
    plottime = downsample(timevec/1000,50)
    for n=1:n_difficulties
        # Find best-fit Ur for this Rc_t[n]
        Ur_b[n] = coolingmodel_Newton_Ur(p,nSteps,Rc_t[n],Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma)
        Tm = coolingmodel(p,Rc_t[n],Ur_b[n],Hmvec,dt,timevec)

        # Plot results
        plotTrel = downsample(Tm-Tm[1],50)
        if n == n_difficulties
            plot!(h,plottime,plotTrel,line_z=difficulies[n],line=(cmap),label="Model")
        else
            plot!(h,plottime,plotTrel,line_z=difficulies[n],line=(cmap),label="")
        end
        # Keep track of results
        Trel = linterp1(timevec/1000,Tm-Tm[1],TrelObs_time)
        R_b[n] = sum((Trel-TrelObs).^2./(2.*TrelObs_sigma.^2))
    end
    plot!(h,legend=:topleft,xlabel="Age (Ma)",ylabel="\\Delta T (C)")
    xlims!(h,(0,4))
    ylims!(h,(-20,250))
    display(h)
    savefig(h,"CoolingCurveVsDifficultyHighRes.pdf")

    # Two-axis plot for best-fit d and Ur
    # figure
    # xlabel('log_{10}(\eta_L/Rc^3)')
    # xlim([minDiff,maxDiff])
    # yyaxis left
    # plot(log10(p.eta_L./Rc_t.^3),Ur_b)
    # ylabel('Best-fit Ur')
    # yyaxis right
    # plot(log10(p.eta_L./Rc_t.^3),R_b)
    # ylabel('Sum squared residual')
    # formatfigure
    # fig = gcf
    # fig.PaperSize = [fig.PaperPosition(3) fig.PaperPosition(4)]
    # saveas(gcf,'Ur=f(Difficulty).pdf')

## --- Explore nl/Rc3 space - smaller range

    nSteps = 20
    dUr = 1E-9
    n_difficulties = 5
    minDiff = 4.5
    maxDiff = 6.5
    difficulies = linspace(minDiff,maxDiff,n_difficulties)

    Rc_t = (p.eta_L ./ (10 .^ difficulies)).^(1/3)
    Ur_b = Array{Float64}(size(Rc_t))
    R_b = Array{Float64}(size(Rc_t))

    cmap = cgrad(flipdim(viridis[1:212],1))
    h = plot(TrelObs_time,TrelObs,yerror=2*TrelObs_sigma,seriestype=:scatter, color=:darkblue, markersize=2, markerstrokecolor=:auto, label="Data")
    plottime = downsample(timevec/1000,50)
    for n=1:n_difficulties
        # Find best-fit Ur for this Rc_t[n]
        Ur_b[n] = coolingmodel_Newton_Ur(p,nSteps,Rc_t[n],Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma)
        Tm = coolingmodel(p,Rc_t[n],Ur_b[n],Hmvec,dt,timevec)

        # Plot results
        plotTrel = downsample(Tm-Tm[1],50)
        if n == n_difficulties
            plot!(h,plottime,plotTrel,line_z=difficulies[n],line=(cmap),label="Model")
        else
            plot!(h,plottime,plotTrel,line_z=difficulies[n],line=(cmap),label="")
        end
        # Keep track of results
        Trel = linterp1(timevec/1000,Tm-Tm[1],TrelObs_time)
        R_b[n] = sum((Trel-TrelObs).^2./(2.*TrelObs_sigma.^2))
    end
    plot!(h,legend=:topleft,xlabel="Age (Ma)",ylabel="\\Delta T (C)")
    xlims!(h,(0,4))
    ylims!(h,(-20,250))
    display(h)
    savefig(h,"CoolingCurveVsDifficulty.pdf")

## --- End of File
