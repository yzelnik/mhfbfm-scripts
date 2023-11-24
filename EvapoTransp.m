function [ ET, PET_PT, PET, E_eq, E_aero, gs, ga, gas ] = EvapoTransp( Q_net, T, t_day, u, z, e_a, LAI, hc, gsmolww, psis)
% by Stefano Manzoni (11/11/2016)
% with additional commments and changes by Giulia Vico (06/02/2023)

% INPUTS (all inputs can be vectors):
% Q_net average daily radiation (in W/m2)
% T air temperature (deg C) 
% t_day day length (hours)
% u wind speed at height z (m/s) 
% z reference height at which u is measured (m)
% e_a partial vapor pressure (kPa)
% LAI leaf area index (m2/m2)
% hc canopy height (m) 
% gsmolww: leaf-level stomatal conductance (mol air/m2 leaf/s), considered to be well watered 
% psis soil water potential (MPa, with sign), to reduce gsmol; no reduction is applied if psis=0

%NOTES:
%gs:
%   reasonable values are of the order of 0.23 mol air/m2 leaf/s for needleleaf evergreen trees, 0.44 mol air/m2 leaf/s for arable cropland
%   for more values, see Appendix 2 of Schulze et al "Relationships among maximum stomatal conductance, 
%   ecosystem surface conductance, carbon assimuilation rate, and plant nitrogen nutrition, Annu Rev Ecol Syst 1994
%   any other value of gs there (in mm/s) needs to be multiplied by 10^(-3) and rho_mol (see below) to get mol air/m2 leaf/s
%   rho_mol decreases with temperature, approx 3% per 10 deg C change in temp
%z:
%   must be at least 3-4 times greater than hc; if z just above hc, ga becomes unrealistically high and Penman Monteith 
%   rates unrealistically high (in particular PET), whereas Priestly Taylor gives a lower than expected value
%the assumption of stable, adiabatic conditions is made, i.e., buoyancy from heating and cooling are neglected; for empirical correction for unstable conditions, see Viney 1991 
%   for possible solutions not requiring H, but requires knowing surface temperature

% --- %example of input:
% Q_net=300; 
% t_day=12;
% T=20;
% u=3; 
% z=100; 
% e_a=1;
% LAI=4;
% hc=15;
% gsmol=0.22;

% OUTPUTS: 
%ACTUAL ET:
% ET Penman Monteith equation for actual evapotranspiration including canopy and aerodynamic conductance (mm/d)
%POTENTIAL ET:
% PET_PT Priestly Taylor approximation of potential evapotranspiration (mm/d)
% PET Penman equation for potential evapotranspiration (canopy conductance not limiting) (mm/d)
% E_eq equilibrium evaporation from Penman Monteith equation (mm/d)
% E_aero aerodynamic component of Penman Monteith equation (mm/d)
%CONDUCTANCES:
% gs canopy conductance (m/s)
% ga aerodynamic conductance (m/s)
% gas combined aerodynamic and canopy conductance (m/s)


% --- physical constants
% constants for the partial pressure of water at saturation esa = a * exp((b.*T)./(T + c)) [kPa, deg C]
a = 0.611;              % [kPa]
b = 17.502;             % [unitless]
c = 240.97;             % [deg C]
%--- Air and water vapor
Rgas = 8.314;              % gas constant [J/mol/K]
m_air = 28.8;           % molecular weight of air [g/mol]
m_w = 18;               % molecular weight of water [g/mol]
epsilon = m_w/m_air;    % ratio of molecular weights=0.622 [-]
rho_w = 10^3;           % liquid water density [kg/m3]
Pair = 101;             % atmospheric pressure [kPa]
Cp = 1012;              % heat capacity of air [J/kg/degC]
%--- other constants
kvonKarman = 0.41;      % von Karman constant
alpha_PT = 1.26;        % Priestly-Taylor constant (Priestly Taylor 1972)


% --- Daily net radiation
Q_net = Q_net*3600*t_day; % J/m2/d

% --- Slope of the saturation water pressure curve (T in deg C)
esa = a .* exp((b.*T)./(T + c)); %kPa
s = a*b*c*exp(b*T./(c+T))./(c+T).^2; %kPa/degC
VPD = esa-e_a; %kPa

% --- latent heat of vaporization (Bonan, 2008)
Lv = (2502-2.308*T)*1000;      % [J/kg]

% --- Psychrometric constant
gamma = Cp*Pair./Lv./epsilon; % kPa/degC

% --- air density
rho_a = Pair./(Rgas/m_air.*(T+273.15)); % kg/m3

% --- aerodynamic conductance 
% Brutsaert, W. 1982. Evaporation in the atmosphere. Theory, history, and applications
% with parameters from Katul et al (2004) One- and two- equation models for canopy turbulence, Boundary layer met,  
% Table 2 values for spruce, Jack pine
d = 2/3*hc; % zero plane displacement height for dense canopies (m)
z0 = 0.1*hc; % momentum roughness height (m) - could be 0.08 for Pinus taeda and hardwood forest
z0v = 0.2*z0; % vapor roughness height (m) - the factor 0.2 is a common one (see Campbell Ch 7, eq 7.19)

ga = u*kvonKarman^2./log((z-d)./z0)./log((z-d)./z0v); % m/s - assuming zero stability corrections, i.e. stable conditions
ra = 1./ga; % s/m

% --- canopy-scale stomatal conductance
rho_mol = 1000*Pair./(Rgas.*(T+273.15)); % molar density of air at atmospheric pressure (mol/m3) - multiply by 1000 to convert kPa to Pa
if psis<0
    gsmol=gs_waterstressed(gsmolww, psis);
else
    gsmol=gsmolww; 
end
gs = gsmol*LAI./rho_mol; % m/s 
rs = 1./gs; % s/m

% --- combined stomatal and aerodynamic conductances
gas = ga.*gs./(ga+gs); % m/s

% --- drying power of air
Ea = epsilon*rho_a.*ga.*VPD./Pair.*3600*t_day; % kg/m2/d

% --- Penman-Monteith equation for daily evapotranspiration
E_mass = s./(s+gamma.*(1+ga./gs)).*(Q_net./Lv)+gamma./(s+gamma.*(1+ga./gs)).*Ea; % kg/m2/d
ET = E_mass./rho_w*1000; % mm/d
E_eq = s./(s+gamma.*(1+ga./gs)).*(Q_net./Lv)./rho_w*1000; % mm/d %evaporation due to the energy warming up the surface and hence the air water holding capacity
E_aero = gamma./(s+gamma.*(1+ga./gs)).*Ea./rho_w*1000; % mm/d %effect of wind, present even in the absence of warming up

% --- Penman equation for daily evapotranspiration (canopy conductance not limiting)
PET_mass = s./(s+gamma).*(Q_net./Lv)+gamma./(s+gamma).*Ea; % kg/m2/d
PET = PET_mass./rho_w*1000; % mm/d

% --- Priestly-Taylor equation for daily potential evapotranspiration
%for comparison with Penman
PET_PT = alpha_PT*E_eq; % mm/d

return


function gsmolws=gs_waterstressed(gsmolww, psi)

%Giulia Vico, Feb 2023
%reduction of stomatal conductance with soil water potential, assuming an exponential decline

%INPUT:
%gsmolww: stomatal conductance under well watered conditions
%psi: soil water potential at which gs ought to be determined (MPa, with sign)

%OUTPUT:
%gsmolws: stomatal conductance at psi


%--- parameter: 
psi50=-0.5; %MPa soil water potential at which gs is reduced by half; if -0.5, only 5% of gsww remains at -1.5 MPa

%--- stomatal conducnatce under water stress
gsmolws=gsmolww*exp(-(psi/(psi50/log(2))));


return

