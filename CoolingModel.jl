    using StatGeochem
    include("Utilities.jl")
    using Plots
    gr();

## --- Define physical parameters

    p = Parameters(0) # Create struct
    p.R = 8.3145; # Universal gas constant (J/Mol/K)
    p.g = 9.8; # Gravitational acceleration (m/s^2)
    p.s_yr = 31556926; # Seconds per year
    p.J_MeV=1.602176565E-13; # Joules per megaelectron-volt
    p.mol=6.022E23; # Avogadro's number (Atoms per mole)

    p.Mm = 4E24; # Mass of mantle (kg)
    p.Mc = 2E24; # Mass of core (kg)
    p.Re = 6371E3; # Radius of Earth (m)
    p.Rp = 1000E3; # Radius of Davies' "plume feeding area"
    p.C_m = 1000;  # Specific heat of mantle (J/kg/C)
    p.C_c = 800;  # Specific heat of core (J/kg/C)
    p.rho_m = 3500; # Density of mantle (kg/m^3)
    p.rho_c = 12300; # Density of core (kg/m^3)

    p.alpha = 2.5E-5; # Thermal expansion coefficient
    p.kappa = 0.86E-6; # Thermal diffusivity (m^2/s)

    p.Trl = 3500+273.15; # Reference temperature, lower mantle (C)
    p.mur_l = 1E19; # Reference viscosity, lower mantle (Pa s)
    p.nc = 1/5; # core exponent
    p.Kc = 30; # Condutivity (W/m/C)

    p.Tru = 1300+273.15; # Reference temperature, upper mantle (C)
    p.mur_u = 1E21; # Reference viscosity, upper mantle (Pa s)
    p.Ea = 400E3; # Activation energy (upper mantle) (J/mol)
    p.Ha = 800E3; # Activation enthalpy (lower mantle) (J/mol)
    p.Km = 3; # Conductivity (W/m/C)
    p.Xm = 1.2; # Mantle heat capacity correction for average temperature
    p.Xc = 1.11; # Core heat capacity correction for average temperature

    p.eta_L = 1E23; # Plate viscsity (Pa s)
    # Done defining parameters

## --- Load observed cooling data to fit against

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
    lambda40K = log(2) / 1.2483E9
    lambda87Rb = log(2) / 4.88E10
    lambda235U = log(2) / 7.04E8
    lambda238U = log(2) / 4.468E9
    lambda232Th = log(2) / 1.405E10

    # Composition: Bulk Silicate Earth
    K40i  = 0.00012 .* mcdbse["K"]  .* exp(lambda40K  .* timevec*10^6);
    Rb87i = 0.27835 .* mcdbse["Rb"] .* exp(lambda87Rb .* timevec*10^6);
    U235i = 0.00720 .* mcdbse["U"]  .* exp(lambda235U .* timevec*10^6);
    U238i = 0.99274 .* mcdbse["U"]  .* exp(lambda238U .* timevec*10^6);
    Th232i = mcdbse["Th"] .* exp(lambda232Th .* timevec*10^6)

    # Mantle heat production in W/Kg, from Rybach, 1985
    Hmvec = 1000 .* K40i   ./10^6 .* lambda40K/p.s_yr  .* (0.6*0.8928 + 1.460*0.1072)*p.J_MeV*p.mol/40 +
        1000 .* Rb87i  ./10^6 .* lambda87Rb/p.s_yr .* 1/3*0.283*p.J_MeV*p.mol/87 +
        1000 .* U235i  ./10^6 .* lambda235U/p.s_yr .* 45.26*p.J_MeV*p.mol/235 +
        1000 .* U238i  ./10^6 .* lambda238U/p.s_yr .* 46.34*p.J_MeV*p.mol/238 +
        1000 .* Th232i ./10^6 .* lambda232Th/p.s_yr .* 38.99*p.J_MeV*p.mol/232

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

## --- Explore different Ur at high Rc

    Rc = (p.eta_L / (10.^8))^(1/3)

    nSteps = 10
    UrMin = 0.1
    UrMax = 0.5
    Ur = linspace(UrMin,UrMax,nSteps)
    cmap = resize_colormap(flipdim(viridis[1:212],1), nSteps)

    h = plot(TrelObs_time,TrelObs,yerror=2*TrelObs_sigma,seriestype=:scatter, markersize=2, markerstrokecolor=:auto, label="Data")
    plot!(h,legend=:topleft,xlabel="Age (Ma)",ylabel="Delta T (C)",fg_color_legend=:white)
    for i=1:nSteps
        Tm = CoolingModel(p,Rc,Ur[i],Hmvec,dt,timevec);
        if i==1
            plot!(h,downsample(timevec/1000,10),downsample(Tm-Tm[1],10),color=cmap[i],label="Model")
        else
            plot!(h,downsample(timevec/1000,10),downsample(Tm-Tm[1],10),color=cmap[i],label="")
        end
    end
    xlims!(h,(0,4))
    ylims!(h, (-20, 300))
    display(h)

    plot!(h,colorbar=true,clims=(UrMin,UrMax),colorbar_title="",color_palette=flipdim(viridis[1:212],1))
    display(h)
    # cb.Label.String = 'Ur';
    savefig(h,"d=8variableUr.pdf")

## --- Explore different Ur at low Rc

    Rc =(p.eta_L./(10.^5.5)).^(1/3);

    nSteps = 10;
    UrMin = 0.7555;
    UrMax = 0.7558;
    Ur = flipud(linspace(UrMin,UrMax,nSteps));

    load viridis
    cmap = flipud(viridis(1:212,:));
    color = interpColor(cmap,nSteps);
    figure; errorbar(TrelObs_time,TrelObs,2*TrelObs_sigma,'.r')
    for i=1:nSteps
        Tm = CoolingModel(p,Rc,Ur[i],Hmvec,dt,timevec);
        hold on; plot(timevec/1000,Tm-Tm[1],'color',color(i,:))
    end

    ylim([-20 250])

    xlabel('Age (Ga)'); ylabel('\Delta T');
    colormap(cmap)
    cb = colorbar;
    cb.Label.String = 'Ur';
    caxis([UrMin,UrMax])
    formatfigure;
    fig = gcf;
    fig.PaperSize = [fig.PaperPosition(3) fig.PaperPosition(4)];
    saveas(gcf,'d=5.5variableUr.pdf')

## --- Explore different Ur at low Rc

    p.eta_L = 1E23; # Plate viscsity (Pa s)
    Rc =(p.eta_L./(10.^5.6)).^(1/3);

    nSteps = 10;
    UrMin = 0.753;
    UrMax = 0.754;
    Ur = flipud(linspace(UrMin,UrMax,nSteps));

    load viridis
    cmap = flipud(viridis(1:212,:));
    color = interpColor(cmap,nSteps);
    figure; errorbar(TrelObs_time,TrelObs,2*TrelObs_sigma,'.r')
    for i=1:nSteps
        Tm = CoolingModel(p,Rc,Ur[i],Hmvec,dt,timevec);
        hold on; plot(timevec/1000,Tm-Tm[1],'color',color(i,:))
    end

    ylim([-20 250])

    xlabel('Age (Ga)'); ylabel('\Delta T');
    colormap(cmap)
    cb = colorbar;
    cb.Label.String = 'Ur';
    caxis([UrMin,UrMax])
    formatfigure;
    fig = gcf;
    fig.PaperSize = [fig.PaperPosition(3) fig.PaperPosition(4)];
    saveas(gcf,'d=5.6variableUr.pdf')


## --- Explore Rc space
#
#     p.eta_L = 1E23; # Plate viscsity (Pa s)
#     nSteps = 20;
#     dUr = 1E-9;
#
#     nDifficulties = 10;
#     Rc_t = linspace(10E3,1000E3,nDifficulties);
#     Ur_b = Array{Float64}(size(Rc_t));
#
#     load viridis
#     color = interpColor(flipud(viridis(1:212,:)),nDifficulties);
#     figure; errorbar(TrelObs_time,TrelObs,2*TrelObs_sigma,'.r')
#     for n=1:nDifficulties
#         Rc = Rc_t(n);
#         Ur = 0.5;
#         for i=1:nSteps-1;
#             ll = CoolingModelLL(p, Rc,Ur,Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma);
#             llp = CoolingModelLL(p, Rc,Ur+dUr,Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma);
#             dll_dUr = (llp-ll)/dUr;
#             Ur = Ur - ll/dll_dUr/2;
#         end
#
#         Ur_b(n) = Ur;
#         Tm = CoolingModel(p,Rc,Ur,Hmvec,dt,timevec);
#         hold on; plot(timevec/1000,Tm-Tm[1],'color',color(n,:))
#     end
#
#     figure; plot(Rc_t,Ur_b)
#     xlabel('Rc'); ylabel('Best fit Ur')
#     figure; plot(log10(p.eta_L./Rc_t.^3),Ur_b)
#     xlabel('\eta_l/Rc'); ylabel('Best fit Ur')


## --- Explore nl/Rc3 space

    p.eta_L = 1E23; # Plate viscsity (Pa s)
    nSteps = 40;
    dUr = 1E-9;

#     nDifficulties = 10;
#     Rc_t = linspace(100E3,1000E3,nDifficulties);

    nDifficulties = 50;
#     nDifficulties = 100;
    minDiff = 4;
    maxDiff = 9;
    Rc_t = (p.eta_L./(10.^linspace(minDiff,maxDiff,nDifficulties))).^(1/3);
    Ur_b = Array{Float64}(size(Rc_t));
    R_b = Array{Float64}(size(Rc_t));

    load viridis
    cmap = flipud(viridis(1:212,:));
    color = interpColor(cmap,nDifficulties);
    figure; errorbar(TrelObs_time,TrelObs,2*TrelObs_sigma,'.r')
    for n=1:nDifficulties
        Rc = Rc_t(n);
        Ur = 0.5;
        for i=1:nSteps-1;
            ll = CoolingModelLL(p, Rc,Ur,Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma);
            llp = CoolingModelLL(p, Rc,Ur+dUr,Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma);
            dll_dUr = (llp-ll)/dUr;
            Ur = Ur - ll/dll_dUr/2;
        end

        Ur_b(n) = Ur;
        Tm = CoolingModel(p,Rc,Ur,Hmvec,dt,timevec);
        hold on; plot(timevec/1000,Tm-Tm[1],'color',color(n,:))
        Trel = interp1(timevec/1000,Tm-Tm[1],TrelObs_time);
        R_b(n) = sum((Trel-TrelObs).^2./(2.*TrelObs_sigma.^2));
    end
    xlabel('Age (Ga)'); ylabel('\Delta T');
    colormap(cmap)
    cb = colorbar;
    cb.Label.String = 'log_{10}(\eta_L/Rc^3)';
    caxis([minDiff,maxDiff])
    formatfigure;
    fig = gcf;
    fig.PaperSize = [fig.PaperPosition(3) fig.PaperPosition(4)];
    set(fig,'renderer','painters')
    saveas(fig,'CoolingCurveVsDifficultyHighRes.pdf')

    figure;
    xlabel('log_{10}(\eta_L/Rc^3)');
    xlim([minDiff,maxDiff])
    yyaxis left
    plot(log10(p.eta_L./Rc_t.^3),Ur_b)
    ylabel('Best-fit Ur');
    yyaxis right
    plot(log10(p.eta_L./Rc_t.^3),R_b)
    ylabel('Sum squared residual');
    formatfigure;
    fig = gcf;
    fig.PaperSize = [fig.PaperPosition(3) fig.PaperPosition(4)];
    saveas(gcf,'Ur=f(Difficulty).pdf')


    ## --- Explore nl/Rc3 space - smaller range

    p.eta_L = 1E23; # Plate viscsity (Pa s)
    nSteps = 20;
    dUr = 1E-9;

    nDifficulties = 5;
    minDiff = 4.5;
    maxDiff = 6.5;
    Rc_t = (p.eta_L./(10.^linspace(minDiff,maxDiff,nDifficulties))).^(1/3);
    Ur_b = Array{Float64}(size(Rc_t));

    load viridis
    cmap = flipud(viridis(1:212,:));
    color = interpColor(cmap,nDifficulties);
    figure; errorbar(TrelObs_time,TrelObs,2*TrelObs_sigma,'.r')
    for n=1:nDifficulties
        Rc = Rc_t(n);
        Ur = 0.5;
        for i=1:nSteps-1;
            ll = CoolingModelLL(p, Rc,Ur,Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma);
            llp = CoolingModelLL(p, Rc,Ur+dUr,Hmvec,dt,timevec,TrelObs_time*1000,TrelObs,TrelObs_sigma);
            dll_dUr = (llp-ll)/dUr;
            Ur = Ur - ll/dll_dUr/2;
        end

        Ur_b[n] = Ur;
        Tm = CoolingModel(p,Rc,Ur,Hmvec,dt,timevec);
        hold on; plot(timevec/1000,Tm-Tm[1],'color',color(n,:))
    end
    xlabel('Age (Ga)'); ylabel('\Delta T (C)');
    colormap(cmap)
    cb = colorbar;
    cb.Label.String = 'log_{10}(\eta_L/Rc^3)';
    caxis([minDiff,maxDiff])
    formatfigure;
    fig = gcf;
    fig.PaperSize = [fig.PaperPosition(3) fig.PaperPosition(4)];
    saveas(gcf,'CoolingCurveVsDifficultyScale.pdf')
    ylim([-20 320])
    delete(cb)
    saveas(gcf,'CoolingCurveVsDifficulty.pdf')
