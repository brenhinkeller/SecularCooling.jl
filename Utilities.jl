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

    eta_L::Float64
end

Parameters(x) = Parameters(x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x)

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
