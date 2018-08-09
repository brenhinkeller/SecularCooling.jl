mutable struct Parameters
    R::Float64 # Universal gas constant (J/Mol/K)
    g::Float64 # Gravitational acceleration (m/s^2)
    s_yr::Float64 # Seconds per year
    J_MeV::Float64 # Joules per megaelectron-volt
    mol::Float64 # Avogadro's number (Atoms per mole)

    Mm::Float64 # Mass of mantle (kg)
    Mc::Float64 # Mass of core (kg)
    Re::Float64 # Radius of Earth (m)
    Rp::Float64 # Radius of "plume feeding area"
    C_m::Float64  # Specific heat of mantle (J/kg/C)
    C_c::Float64  # Specific heat of mantle (J/kg/C)
    rho_m::Float64 # Density of mantle (kg/m^3)
    rho_c::Float64 # Density of mantle (kg/m^3)

    alpha::Float64 # Thermal expansion coefficient
    kappa::Float64 # Thermal diffusivity (m^2/s)

    Trl::Float64 # Reference temperature, lower mantle (C)
    mur_l::Float64 # Reference viscosity, lower mantle (Pa s)
    nc::Float64 # core exponent
    Kc::Float64 # Condutivity (W/m/C)

    Tru::Float64 # Reference temperature, upper mantle (C)
    mur_u::Float64 # Reference viscosity, upper mantle (Pa s)
    Ea::Float64 # Activation energy (upper mantle) (J/mol)
    Ha::Float64 # Activation enthalpy (lower mantle) (J/mol)
    Km::Float64 # Conductivity (W/m/C)
    Xm::Float64 # Mantle heat capacity correction for average temperature
    Xc::Float64 # Core heat capacity correction for average temperature

    eta_L::Float64 # Lithospheric viscosity (Pa s)

    lambda40K::Float64   # K-40 decay constant (1/years)
    lambda87Rb::Float64  # Rb-87 decay constant (1/years)
    lambda235U::Float64  # U-235 decay constant (1/years)
    lambda238U::Float64  # U-338 decay constant (1/years)
    lambda232Th::Float64 # Th-232 decay constant (1/years)
end

Parameters(x) = Parameters(x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x)
Parameters() = Parameters(8.3145, 9.8, 3.1556926e7, 1.602176565e-13, 6.022e23, 4.0e24, 2.0e24, 6.371e6, 1.0e6, 1000.0, 800.0, 3500.0, 12300.0, 2.5e-5, 8.6e-7, 3773.15, 1.0e19, 0.2, 30.0, 1573.15, 1.0e21, 400000.0, 800000.0, 3.0, 1.2, 1.11, 1.0e23, 5.552729156131902e-10, 1.4203835667211993e-11, 9.845840632953767e-10, 1.551358953804712e-10, 4.933431890106372e-11)

    # Define parameters
    # p = Parameters(0) # Create struct
    # p.R = 8.3145; # Universal gas constant (J/Mol/K)
    # p.g = 9.8; # Gravitational acceleration (m/s^2)
    # p.s_yr = 31556926; # Seconds per year
    # p.J_MeV=1.602176565E-13; # Joules per megaelectron-volt
    # p.mol=6.022E23; # Avogadro's number (Atoms per mole)
    #
    # p.Mm = 4E24; # Mass of mantle (kg)
    # p.Mc = 2E24; # Mass of core (kg)
    # p.Re = 6371E3; # Radius of Earth (m)
    # p.Rp = 1000E3; # Radius of Davies' "plume feeding area"
    # p.C_m = 1000;  # Specific heat of mantle (J/kg/C)
    # p.C_c = 800;  # Specific heat of core (J/kg/C)
    # p.rho_m = 3500; # Density of mantle (kg/m^3)
    # p.rho_c = 12300; # Density of core (kg/m^3)
    #
    # p.alpha = 2.5E-5; # Thermal expansion coefficient
    # p.kappa = 0.86E-6; # Thermal diffusivity (m^2/s)
    #
    # p.Trl = 3500+273.15; # Reference temperature, lower mantle (C)
    # p.mur_l = 1E19; # Reference viscosity, lower mantle (Pa s)
    # p.nc = 1/5; # core exponent
    # p.Kc = 30; # Condutivity (W/m/C)
    #
    # p.Tru = 1300+273.15; # Reference temperature, upper mantle (C)
    # p.mur_u = 1E21; # Reference viscosity, upper mantle (Pa s)
    # p.Ea = 400E3; # Activation energy (upper mantle) (J/mol)
    # p.Ha = 800E3; # Activation enthalpy (lower mantle) (J/mol)
    # p.Km = 3; # Conductivity (W/m/C)
    # p.Xm = 1.2; # Mantle heat capacity correction for average temperature
    # p.Xc = 1.11; # Core heat capacity correction for average temperature
    #
    # p.eta_L = 1E23; # Plate viscsity (Pa s)
    #
    # p.lambda40K  = log(2) / 1.2483E9 # K-40 decay constant (1/years)
    # p.lambda87Rb = log(2) / 4.88E10  # Rb-87 decay constant (1/years)
    # p.lambda235U = log(2) / 7.04E8   # U-235 decay constant (1/years)
    # p.lambda238U = log(2) / 4.468E9  # U-238 decay constant (1/years)
    # p.lambda232Th = log(2) / 1.405E10 # Th-232 decay constant (1/years)
    # # Done defining parameters


function CoolingModel(p::Parameters, Rc, Ur, Hm, dt, timevec)
    a_c = 4pi * Rc^2 * p.Kc * ((2 * p.g * p.rho_c * p.alpha) / (3 * p.kappa * p.Rp^2))^p.nc

    # Allocate arrays
    Qc = Array{Float64}(size(timevec))
    Qm = Array{Float64}(size(timevec))
    dTm_dt = Array{Float64}(size(timevec))
    dTc_dt = Array{Float64}(size(timevec))
    mu_u = Array{Float64}(size(timevec))
    mu_l = Array{Float64}(size(timevec))
    Tm = Array{Float64}(size(timevec))
    Tc = Array{Float64}(size(timevec))
    dp = Array{Float64}(size(timevec))

    # Variable parameters to adjust

    Qm[1] = 35*10^12 # Present mantle heat flux (W)
    Qc[1] = 6*10^12 # Present core heat flux (W)
    Tm[1] = 1300+273.15 # Present mantle temperature (K)
    Tc[1] = 4400+273.15 # Present core temperature (K)
    # dp[1] = 120E3 + 275*(Tm[1]-(1400+273.15))  # Old equation from Davies
    dp[1] = 115E3 + 275*(Tm[1]-(1400+273.15)) + (Tm[1]<(1410+273.15)).*(0.003649*(Tm[1]-273.15).^2 -10.4*(Tm[1]-273.15) + 7409.5) # Piecewise fit to Korenaga
    mu_u[1] = p.mur_u * exp(p.Ea/p.R*(1/Tm[1]-1/p.Tru))
    mu_l[1] = p.mur_l * exp(p.Ha/p.R*(1/Tc[1]-1/p.Trl))
    dTc_dt[1] = (p.Mc.*0 - Qc[1]) / (p.Xc * p.Mc * p.C_c) * p.s_yr*1E6
    epsilon = Tm[1]^2 / (p.Ha/p.R*(Tc[1]-Tm[1]))
    Beta_c = Qc[1] /
        (a_c * epsilon^(4*p.nc) * ((Tc[1]-Tm[1])^(1+p.nc) / mu_l[1]^p.nc))
    Beta_d = Qm[1] /
        (4*pi*p.Re^2*p.Km*Tm[1]*
            (p.g*p.rho_m*p.alpha*Tm[1]*dp[1] /
                (4*mu_u[1]*p.kappa + (4/3*p.eta_L*p.kappa)*(dp[1]/Rc)^3)
            )^(1/2)
        )
    Hm = Hm .* Ur / (Hm[1]*p.Mm/Qm[1])
    dTm_dt[1] = (p.Mm.*Hm[1] + Qc[1] - Qm[1]) /
        (p.Xm * p.Mm * p.C_m) * p.s_yr * 1E6

    for i=2:length(timevec)
        Tm[i] = Tm[i-1] - dTm_dt[i-1]*dt # New mantle temperature
        Tc[i] = Tc[i-1] - dTc_dt[i-1]*dt # New core temperature
        mu_l[i] = p.mur_l * exp(p.Ha/p.R*(1/Tc[i]-1/p.Trl))
        epsilon = Tm[i]^2 / (p.Ha/p.R*(Tc[i]-Tm[i]))
        Qc[i] = Beta_c * a_c * abs(epsilon)^(4*p.nc) *
            abs(Tc[i]-Tm[i])^(1+p.nc) / abs(mu_l[i])^p.nc
        mu_u[i] = p.mur_u*exp(p.Ea/p.R*(1/Tm[i]-1/p.Tru))
        # dp[i] = 120E3 + 275*(Tm[i]-(1400+273.15)) # Old equation from Davies
        # Piecewise fit to Korenaga plate scaling
        dp[i] = 115E3 + 275*(Tm[i]-(1400+273.15)) +
            (Tm[i]<(1410+273.15)) .*
            (0.003649*(Tm[i]-273.15).^2 -10.4*(Tm[i]-273.15) + 7409.5)
        Qm[i] = Beta_d*4*pi*p.Re^2*p.Km*Tm[i] *
                sqrt(abs(p.g*p.rho_m*p.alpha*Tm[i]*dp[i]/(4*mu_u[i]*p.kappa +
                    (4/3*p.eta_L*p.kappa)*(dp[i]/Rc)^3)
                ))
        dTm_dt[i] = (p.Mm.*Hm[i] + Qc[i] - Qm[i]) /
            (p.Xm * p.Mm * p.C_m) * p.s_yr * 1E6
        dTc_dt[i] = (p.Mc.*0 - Qc[i]) /
            (p.Xc * p.Mc * p.C_c) * p.s_yr * 1E6
    end

    Tm
end

function sLL(Tm,timevec,TrelAve_time,TrelAve,TrelAve_sigma)
    Trel = linterp1(timevec,Tm-Tm[1],TrelAve_time)
    sll = sum(sign.(Trel - TrelAve) .*
            (Trel - TrelAve).^2./(2.*TrelAve_sigma.^2)
            -log.(sqrt.(2*pi*TrelAve_sigma))
        )
    return sll
end


function CoolingModelLL(p::Parameters, Rc, Ur, Hm, dt, timevec, TrelAve_time, TrelAve, TrelAve_sigma)
    Tm = CoolingModel(p, Rc, Ur, Hmvec, dt, timevec)
    return sLL(Tm, timevec, TrelAve_time, TrelAve, TrelAve_sigma)
end
