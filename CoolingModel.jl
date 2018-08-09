    1+1
## ---

    using StatGeochem
    include("Utilities.jl")
    using Plots; gr()

## --- Define physical parameters

    p = Parameters()
    p.eta_L = 1E23; # Plate viscsity (Pa s)

## --- Load observed cooling against which to fit

    using MAT
    mcdbse = matread("mcdbse.mat")["mcdbse"]
    obs = matread("CoolingAverage.mat")
    TrelObs_time = obs["c"][:]
    TrelObs = obs["TrelAve"][:]
    TrelObs_sigma = obs["TrelAve_sigma"][:]
    TrelObs_sigma[1] = 0.5;

## --- Set up model

    dt = 0.1; # Time step (Myr)
    timevec = collect(0:dt:4000);

    # Heat production
    # Composition: Bulk Silicate Earth
    K40i  = 0.00012 .* mcdbse["K"]  .* exp.(p.lambda40K  .* timevec*10^6);
    Rb87i = 0.27835 .* mcdbse["Rb"] .* exp.(p.lambda87Rb .* timevec*10^6);
    U235i = 0.00720 .* mcdbse["U"]  .* exp.(p.lambda235U .* timevec*10^6);
    U238i = 0.99274 .* mcdbse["U"]  .* exp.(p.lambda238U .* timevec*10^6);
    Th232i = mcdbse["Th"] .* exp.(p.lambda232Th .* timevec*10^6)

    # Mantle heat production in W/Kg, from Rybach, 1985
    Hmvec = 1000 .* K40i   ./10^6 .* p.lambda40K/p.s_yr  .* (0.6*0.8928 + 1.460*0.1072)*p.J_MeV*p.mol/40 +
        1000 .* Rb87i  ./10^6 .* p.lambda87Rb/p.s_yr .* 1/3*0.283*p.J_MeV*p.mol/87 +
        1000 .* U235i  ./10^6 .* p.lambda235U/p.s_yr .* 45.26*p.J_MeV*p.mol/235 +
        1000 .* U238i  ./10^6 .* p.lambda238U/p.s_yr .* 46.34*p.J_MeV*p.mol/238 +
        1000 .* Th232i ./10^6 .* p.lambda232Th/p.s_yr .* 38.99*p.J_MeV*p.mol/232

## --- Run Newton's method

    Rc = 750E3 # Plate bending radius (km)
    dUr = 1E-9 # Urey ratio step size

    nSteps = 20
    Ur = Array{Float64}(nSteps)
    Ur[1] = 0.5
    for i = 1:nSteps-1
        ll = CoolingModelLL(p,Rc,Ur[i],Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma);
        llp = CoolingModelLL(p,Rc,Ur[i]+dUr,Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma);
        dll_dUr = (llp-ll)/dUr;
        Ur[i+1] = Ur[i] - ll/dll_dUr/2;
    end

    Tm = CoolingModel(p, Rc, Ur[nSteps], Hmvec, dt, timevec)
    h = plot(TrelObs_time, TrelObs, yerror=2*TrelObs_sigma, seriestype=:scatter, markersize=2, markerstrokecolor=:auto, label="data")
    plot!(h, timevec/1000,Tm-Tm[1], label="model", xlabel="Age (Ma)", ylabel="Delta T (C)", legend=:topleft)
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
    # cmap = resize_colormap(flipdim(viridis[1:212],1), nSteps)
    cmap = cgrad(flipdim(viridis[1:212],1))
    h = plot(TrelObs_time,TrelObs,yerror=2*TrelObs_sigma,seriestype=:scatter, markersize=2, color=:darkred, markerstrokecolor=:auto, label="Data")
    plottime = downsample(timevec/1000,50)
    for i=1:nSteps
        Tm = CoolingModel(p,Rc,Ur[i],Hmvec,dt,timevec);
        plotTrel = downsample(Tm-Tm[1],50)
        if i==1
            plot!(h,plottime,plotTrel,line_z=Ur[i],line=(cmap),label="Model") #color=cmap[i]
        else
            plot!(h,plottime,plotTrel,line_z=Ur[i],line=(cmap),label="") #color=cmap[i]
        end
    end
    plot!(h,legend=:topleft,xlabel="Age (Ma)",ylabel="Delta T (C)")
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
    cmap = cgrad(flipdim(viridis[1:212],1))

    h = plot(TrelObs_time,TrelObs,yerror=2*TrelObs_sigma,seriestype=:scatter, markersize=2, color=:darkred, markerstrokecolor=:auto, label="Data")
    plottime = downsample(timevec/1000,50)
    for i=1:nSteps
        Tm = CoolingModel(p,Rc,Ur[i],Hmvec,dt,timevec);
        plotTrel = downsample(Tm-Tm[1],50)
        if i==1
            plot!(h,plottime,plotTrel,line_z=Ur[i],line=(cmap),label="Model") #color=cmap[i]
        else
            plot!(h,plottime,plotTrel,line_z=Ur[i],line=(cmap),label="") #color=cmap[i]
        end
    end
    plot!(h,legend=:topleft,xlabel="Age (Ma)",ylabel="Delta T (C)")
    xlims!(h,(0,4))
    ylims!(h,(-20,250))
    display(h)
    savefig(h,"d=5.5variableUr.pdf")


## --- Explore different Ur at low d (5.6)

    nSteps = 10;
    Rc = (p.eta_L / (10^5.6))^(1/3)
    UrMin = 0.753
    UrMax = 0.754

    Ur = reverse(linspace(UrMin,UrMax,nSteps))
    cmap = cgrad(flipdim(viridis[1:212],1))

    h = plot(TrelObs_time,TrelObs,yerror=2*TrelObs_sigma,seriestype=:scatter, color=:darkred, markersize=2, markerstrokecolor=:auto, label="Data")
    plottime = downsample(timevec/1000,50)
    for i=1:nSteps
        # Calculate the mantle temperature curve
        Tm = CoolingModel(p,Rc,Ur[i],Hmvec,dt,timevec)

        # Plot the resutls
        plotTrel = downsample(Tm-Tm[1],50)
        if i==1
            plot!(h,plottime,plotTrel,line_z=Ur[i],line=(cmap),label="Model")
        else
            plot!(h,plottime,plotTrel,line_z=Ur[i],line=(cmap),label="")
        end
    end
    plot!(h,legend=:topleft,xlabel="Age (Ma)",ylabel="Delta T (C)")
    xlims!(h,(0,4))
    ylims!(h,(-20,250))
    display(h)
    savefig(h,"d=5.6variableUr.pdf")

## --- Explore nl/Rc3 space

    p.eta_L = 1E23 # Plate viscsity (Pa s)
    nSteps = 40
    dUr = 1E-9

    n_difficulties = 50
    minDiff = 4
    maxDiff = 9
    difficulies = linspace(minDiff,maxDiff,n_difficulties)

    Rc_t = (p.eta_L ./ (10 .^ difficulies)).^(1/3)
    Ur_b = Array{Float64}(size(Rc_t))
    R_b = Array{Float64}(size(Rc_t))

    cmap = cgrad(flipdim(viridis[1:212],1))
    h = plot(TrelObs_time,TrelObs,yerror=2*TrelObs_sigma,seriestype=:scatter, color=:darkred, markersize=2, markerstrokecolor=:auto, label="Data")
    plottime = downsample(timevec/1000,50)
    for n=1:n_difficulties
        Rc = Rc_t[n]
        Ur = 0.5
        for i=1:nSteps-1
            ll = CoolingModelLL(p, Rc,Ur,Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma)
            llp = CoolingModelLL(p, Rc,Ur+dUr,Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma)
            dll_dUr = (llp-ll)/dUr
            Ur = Ur - ll/dll_dUr/2
        end

        Ur_b[n] = Ur
        Tm = CoolingModel(p,Rc,Ur,Hmvec,dt,timevec)

        plotTrel = downsample(Tm-Tm[1],50)
        if n == n_difficulties
            plot!(h,plottime,plotTrel,line_z=difficulies[n],line=(cmap),label="Model")
        else
            plot!(h,plottime,plotTrel,line_z=difficulies[n],line=(cmap),label="")
        end
        Trel = linterp1(timevec/1000,Tm-Tm[1],TrelObs_time)
        R_b[n] = sum((Trel-TrelObs).^2./(2.*TrelObs_sigma.^2))
    end
    plot!(h,legend=:topleft,xlabel="Age (Ma)",ylabel="Delta T (C)")
    xlims!(h,(0,4))
    ylims!(h,(-20,250))
    display(h)
    savefig(h,"CoolingCurveVsDifficultyHighRes.pdf")


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
    h = plot(TrelObs_time,TrelObs,yerror=2*TrelObs_sigma,seriestype=:scatter, color=:darkred, markersize=2, markerstrokecolor=:auto, label="Data")
    plottime = downsample(timevec/1000,50)
    for n=1:n_difficulties
        Rc = Rc_t[n]
        Ur = 0.5
        for i=1:nSteps-1
            ll = CoolingModelLL(p, Rc,Ur,Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma)
            llp = CoolingModelLL(p, Rc,Ur+dUr,Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma)
            dll_dUr = (llp-ll)/dUr
            Ur = Ur - ll/dll_dUr/2
        end

        Ur_b[n] = Ur
        Tm = CoolingModel(p,Rc,Ur,Hmvec,dt,timevec)

        plotTrel = downsample(Tm-Tm[1],50)
        if n == n_difficulties
            plot!(h,plottime,plotTrel,line_z=difficulies[n],line=(cmap),label="Model")
        else
            plot!(h,plottime,plotTrel,line_z=difficulies[n],line=(cmap),label="")
        end
        Trel = linterp1(timevec/1000,Tm-Tm[1],TrelObs_time)
        R_b[n] = sum((Trel-TrelObs).^2./(2.*TrelObs_sigma.^2))
    end
    plot!(h,legend=:topleft,xlabel="Age (Ma)",ylabel="Delta T (C)")
    xlims!(h,(0,4))
    ylims!(h,(-20,250))
    display(h)
    savefig(h,"CoolingCurveVsDifficulty.pdf")
