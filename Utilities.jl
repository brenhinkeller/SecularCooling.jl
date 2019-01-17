## --- Load observed cooling against which to fit

    using MAT
    mcdbse = matread("mcdbse.mat")["mcdbse"]
    obs = matread("CoolingAverage.mat")
    TrelObs_time = obs["c"][:]
    TrelObs = obs["TrelAve"][:]
    TrelObs_sigma = obs["TrelAve_sigma"][:]
    TrelObs_sigma[1] = 0.5;

## -- Constant parameters

    mutable struct Parameters
        R::Float64 # Universal gas constant (J/Mol/K)
        g::Float64 # Gravitational acceleration (m/s^2)
        s_yr::Float64 # Seconds per year
        J_MeV::Float64 # Joules per megaelectron-volt
        mol::Float64 # Avogadro's number (Atoms per mole)
        lambda40K::Float64   # K-40 decay constant (1/years)
        lambda87Rb::Float64  # Rb-87 decay constant (1/years)
        lambda235U::Float64  # U-235 decay constant (1/years)
        lambda238U::Float64  # U-338 decay constant (1/years)
        lambda232Th::Float64 # Th-232 decay constant (1/years)

        Mm::Float64 # Mass of mantle (kg)
        Mc::Float64 # Mass of core (kg)
        Re::Float64 # Radius of Earth (m)
        Rp::Float64 # Radius of "plume feeding area"
        Cp_m::Float64  # Specific heat of mantle (J/kg/C)
        Cp_c::Float64  # Specific heat of mantle (J/kg/C)
        rho_m::Float64 # Density of mantle (kg/m^3)
        rho_c::Float64 # Density of mantle (kg/m^3)
        Qm_now::Float64 # Present mantle heat flux (W)
        Qc_now::Float64 # Present core heat flux (W)
        Tm_now::Float64 # Present mantle temperature (K)
        Tc_now::Float64 # Present core temperature (K)

        alpha::Float64 # Thermal expansion coefficient
        kappa::Float64 # Thermal diffusivity (m^2/s)

        Trl::Float64 # Reference temperature, CMB (C)
        mur_l::Float64 # Reference viscosity, CMB (Pa s)
        nc::Float64 # Core exponent
        Kc::Float64 # Condutivity of the core (W/m/C)

        Tru::Float64 # Reference temperature, upper mantle (C)
        mur_u::Float64 # Reference viscosity, upper mantle (Pa s)
        Ea::Float64 # Activation energy, upper mantle (J/mol)
        Ha::Float64 # Activation enthalpy, CMB (J/mol)
        Km::Float64 # Conductivity of the mantle (W/m/C)
        Xm::Float64 # Mantle heat capacity correction for average temperature
        Xc::Float64 # Core heat capacity correction for average temperature

        eta_L::Float64 # Lithospheric viscosity (Pa s)
        Rc::Float64 # Lithspheric plate bending radius of curvature (m)
        Q_addtl::Float64 # Constant additional heat flux to base of mantle (for sensitivity testing)
    end

    Parameters() = Parameters(8.3145, 9.8, 3.1556926e7, 1.602176565e-13, 6.022e23, 5.552729156131902e-10, 1.4203835667211993e-11, 9.845840632953767e-10, 1.551358953804712e-10, 4.933431890106372e-11, 4.0e24, 2.0e24, 6.371e6, 1.0e6, 1000.0, 800.0, 3500.0, 12300.0, 3.5e13, 6.0e12, 1573.15, 4673.15, 2.5e-5, 8.6e-7, 3773.15, 1.0e19, 0.2, 30.0, 1573.15, 1.0e21, 400000.0, 800000.0, 3.0, 1.2, 1.11, 1.0e23, 450000.0, 0.0)

    # # Define parameters
    # p = Parameters() # Create struct
    # p.R = 8.3145; # Universal gas constant (J/Mol/K)
    # p.g = 9.8; # Gravitational acceleration (m/s^2)
    # p.s_yr = 31556926; # Seconds per year
    # p.J_MeV=1.602176565E-13; # Joules per megaelectron-volt
    # p.mol=6.022E23; # Avogadro's number (Atoms per mole)
    # p.lambda40K  = log(2) / 1.2483E9 # K-40 decay constant (1/years)
    # p.lambda87Rb = log(2) / 4.88E10  # Rb-87 decay constant (1/years)
    # p.lambda235U = log(2) / 7.04E8   # U-235 decay constant (1/years)
    # p.lambda238U = log(2) / 4.468E9  # U-238 decay constant (1/years)
    # p.lambda232Th = log(2) / 1.405E10 # Th-232 decay constant (1/years)
    #
    # p.Mm = 4E24; # Mass of mantle (kg)
    # p.Mc = 2E24; # Mass of core (kg)
    # p.Re = 6371E3; # Radius of Earth (m)
    # p.Rp = 1000E3; # Radius of Davies' "plume feeding area"
    # p.Cp_m = 1000;  # Specific heat of mantle (J/kg/C)
    # p.Cp_c = 800;  # Specific heat of core (J/kg/C)
    # p.rho_m = 3500; # Density of mantle (kg/m^3)
    # p.rho_c = 12300; # Density of core (kg/m^3)
    # p.Qm_now = 35*10^12 # Present mantle heat flux (W)
    # p.Qc_now = 6*10^12 # Present core heat flux (W)
    # p.Tm_now = 1300+273.15 # Present mantle temperature (K)
    # p.Tc_now = 4400+273.15 # Present core temperature (K). I/O boundary 6230 +/- 500 according to 10.1126/science.1233514
    #
    # p.alpha = 2.5E-5; # Thermal expansion coefficient
    # p.kappa = 0.86E-6; # Thermal diffusivity (m^2/s)
    #
    # p.Trl = 3500+273.15; # Reference temperature, CMB (C)
    # p.mur_l = 1E19; # Reference viscosity, CMB (Pa s)
    # p.nc = 1/5; # Core exponent
    # p.Kc = 30; # Condutivity of the core (W/m/C)
    #
    # p.Tru = 1300+273.15; # Reference temperature, upper mantle (C)
    # p.mur_u = 1E21; # Reference viscosity, upper mantle (Pa s)
    # p.Ea = 400E3; # Activation energy, upper mantle (J/mol)
    # p.Ha = 800E3; # Activation enthalpy CMB (J/mol)
    # p.Km = 3; # Conductivity of the mantle (W/m/C)
    # p.Xm = 1.2; # Mantle heat capacity correction for average temperature
    # p.Xc = 1.11; # Core heat capacity correction for average temperature
    #
    # p.eta_L = 1E23; # Plate viscsity (Pa s)
    # p.Rc = 450000; # Lithspheric plate bending radius of curvature (m)
    # p.Q_addtl = 0 # Additional (constant) heat flux
    # # Done defining parameters

## --- Calculate past mantle heat production model for a given time vector

    function mantle_heat_production(p::Parameters,composition::Dict,timevec)
        # Composition: Bulk Silicate Earth
        K40i  = 0.00012 .* composition["K"]  .* exp.(p.lambda40K  .* timevec*10^6);
        Rb87i = 0.27835 .* composition["Rb"] .* exp.(p.lambda87Rb .* timevec*10^6);
        U235i = 0.00720 .* composition["U"]  .* exp.(p.lambda235U .* timevec*10^6);
        U238i = 0.99274 .* composition["U"]  .* exp.(p.lambda238U .* timevec*10^6);
        Th232i = composition["Th"] .* exp.(p.lambda232Th .* timevec*10^6)

        # Mantle heat production in W/Kg, from Rybach, 1985
        Hmvec = 1000 .*K40i./10^6 .* p.lambda40K/p.s_yr  .* (0.6*0.8928 + 1.460*0.1072)*p.J_MeV*p.mol/40 +
            1000 .* Rb87i  ./10^6 .* p.lambda87Rb/p.s_yr .* 1/3*0.283*p.J_MeV*p.mol/87 +
            1000 .* U235i  ./10^6 .* p.lambda235U/p.s_yr .* 45.26*p.J_MeV*p.mol/235 +
            1000 .* U238i  ./10^6 .* p.lambda238U/p.s_yr .* 46.34*p.J_MeV*p.mol/238 +
            1000 .* Th232i ./10^6 .* p.lambda232Th/p.s_yr .* 38.99*p.J_MeV*p.mol/232

        return Hmvec
    end

## --- Model functions

    # Finite-difference cooling model
    function coolingmodel(p::Parameters, Ur, Hm, dt, timevec)

        # Allocate arrays
        Qc = Array{Float64}(undef, size(timevec))
        Qm = Array{Float64}(undef, size(timevec))
        dTm_dt = Array{Float64}(undef, size(timevec))
        dTc_dt = Array{Float64}(undef, size(timevec))
        mu_u = Array{Float64}(undef, size(timevec))
        mu_l = Array{Float64}(undef, size(timevec))
        Tm = Array{Float64}(undef, size(timevec))
        Tc = Array{Float64}(undef, size(timevec))
        dp = Array{Float64}(undef, size(timevec))

        # Consolidated core constants
        a_c = 4pi * p.Rc^2 * p.Kc * ((2 * p.g * p.rho_c * p.alpha) / (3 * p.kappa * p.Rp^2))^p.nc

        # ---  Variable parameters to adjust --- #
        Qm[1] = p.Qm_now # Present mantle heat flux (W)
        Qc[1] = p.Qc_now # Present core heat flux (W)
        Tm[1] = p.Tm_now # Present mantle temperature (K)
        Tc[1] = p.Tc_now # Present core temperature (K)

        mu_u[1] = p.mur_u * exp(p.Ea/p.R*(1/Tm[1]-1/p.Tru)) # Inital mantle viscosity
        mu_l[1] = p.mur_l * exp(p.Ha/p.R*(1/Tc[1]-1/p.Trl)) # Initial lower mantle / CMB viscosity
        epsilon = Tm[1]^2 / (p.Ha/p.R*(Tc[1]-Tm[1])) # Initial CMB boundary layer thickness ratio

        # Intial lithospheric plate thickness (piecewise fit to Korenaga)
        dp[1] = 115E3 + 275*(Tm[1]-(1400+273.15)) +
            (Tm[1]<(1410+273.15)) *
            (0.003649*(Tm[1]-273.15).^2 - 10.4*(Tm[1]-273.15) + 7409.5)

        # Adjustment to fit modern (imposed) core and mantle heat fluxes
        Beta_c = Qc[1] /
            (a_c * epsilon^(4*p.nc) * ((Tc[1]-Tm[1])^(1+p.nc) / mu_l[1]^p.nc))
        Beta_d = Qm[1] /
            (4*pi*p.Re^2*p.Km*Tm[1]*
                (p.g*p.rho_m*p.alpha*Tm[1]*dp[1] /
                    (4*mu_u[1]*p.kappa + (4/3*p.eta_L*p.kappa)*(dp[1]/p.Rc)^3)
                )^(1/2)
            )

        # Scale mantle heat production to satisfy Urey ratio
        Hm = Hm .* Ur / (Hm[1]*p.Mm/Qm[1])

        # Initial core and mantle cooling rates
        dTc_dt[1] = (p.Mc.*0 - Qc[1]) /
            (p.Xc * p.Mc * p.Cp_c) * p.s_yr*1E6
        dTm_dt[1] = (p.Mm.*Hm[1] + Qc[1] - Qm[1]) /
            (p.Xm * p.Mm * p.Cp_m) * p.s_yr * 1E6

        for i=2:length(timevec)
            Tm[i] = Tm[i-1] - dTm_dt[i-1]*dt # New mantle temperature
            Tc[i] = Tc[i-1] - dTc_dt[i-1]*dt # New core temperature

            mu_u[i] = p.mur_u * exp(p.Ea/p.R*(1/Tm[i]-1/p.Tru)) # New mantle viscosity
            mu_l[i] = p.mur_l * exp(p.Ha/p.R*(1/Tc[i]-1/p.Trl)) # New lower mantle / CMB viscosity
            epsilon = Tm[i]^2 / (p.Ha/p.R*(Tc[i]-Tm[i])) # New CMB boundary layer thickness ratio

            # Lithospheric plate thickness (piecewise fit to Korenaga)
            dp[i] = 115E3 + 275*(Tm[i]-(1400+273.15)) +
                (Tm[i]<(1410+273.15)) *
                (0.003649*(Tm[i]-273.15).^2 - 10.4*(Tm[i]-273.15) + 7409.5)

            # Core heat flux
            Qc[i] = Beta_c * a_c * abs(epsilon)^(4*p.nc) *
                abs(Tc[i]-Tm[i])^(1+p.nc) / abs(mu_l[i])^p.nc

            # Mantle heat flux
            Qm[i] = Beta_d*4*pi*p.Re^2*p.Km*Tm[i] *
                sqrt(abs(p.g*p.rho_m*p.alpha*Tm[i]*dp[i]/(4*mu_u[i]*p.kappa +
                    (4/3*p.eta_L*p.kappa)*(dp[i]/p.Rc)^3)
                ))

            # Core and mantle cooling rates by energy balance
            dTc_dt[i] = (p.Mc.*0 - Qc[i]) /
                (p.Xc * p.Mc * p.Cp_c) * p.s_yr * 1E6
            dTm_dt[i] = (p.Mm.*Hm[i] + Qc[i] + p.Q_addtl - Qm[i]) /
                (p.Xm * p.Mm * p.Cp_m) * p.s_yr * 1E6
        end

        return (Tm, Tc, Qm, Qc, dp)
    end

    # Return the signed log-likelihood for a given temperature curve
    function signed_LL(Tm, timevec, TrelAve_time, TrelAve, TrelAve_sigma)
        Trel = linterp1(timevec, Tm .- Tm[1], TrelAve_time)
        signed_LL = sum(sign.(Trel - TrelAve) .*
                (Trel - TrelAve).^2 ./ (2 .* TrelAve_sigma.^2)
                .- log.(sqrt.(2*pi*TrelAve_sigma))
            )
        return signed_LL
    end

    # Run model and calculate its signed log likelihood
    function coolingmodel_LL(p::Parameters, Ur, Hm, dt, timevec, TrelAve_time, TrelAve, TrelAve_sigma)
        (Tm, Tc, Qm, Qc, dp) = coolingmodel(p, Ur, Hmvec, dt, timevec)
        return signed_LL(Tm, timevec, TrelAve_time, TrelAve, TrelAve_sigma)
    end

    # Return the residual for a given temperature curve
    function resid(Tm, timevec, TrelAve_time, TrelAve, TrelAve_sigma)
        Trel = linterp1(timevec, Tm .- Tm[1], TrelAve_time)
        return sum((Trel-TrelObs).^2 ./ (2 .* TrelObs_sigma.^2))
    end

    # Run model and calculate its residual
    function coolingmodel_resid(p::Parameters, Ur, Hm, dt, timevec, TrelAve_time, TrelAve, TrelAve_sigma)
        (Tm, Tc, Qm, Qc, dp) = coolingmodel(p, Ur, Hmvec, dt, timevec)
        return resid(Tm, timevec, TrelAve_time, TrelAve, TrelAve_sigma)
    end

    function coolingmodel_Newton_Ur(p::Parameters,nSteps,Hmvec,dt,timevec,TrelObs_time,TrelObs,TrelObs_sigma)
        Ur = 0.5
        dUr = 1E-9
        for i=1:nSteps-1
            ll = coolingmodel_LL(p,Ur,Hmvec,dt,timevec,TrelObs_time,TrelObs,TrelObs_sigma)
            llp = coolingmodel_LL(p,Ur+dUr,Hmvec,dt,timevec,TrelObs_time,TrelObs,TrelObs_sigma)
            dll_dUr = (llp-ll)/dUr
            Ur = Ur - ll/dll_dUr/2
        end
        return Ur
    end
