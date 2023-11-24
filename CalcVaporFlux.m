function VF=CalcVaporFlux(spechumids,soilwcs,adr,TA)
% by Yuval Zelnik (24/11/2023)
% Calculate vapor flux (from soil to dead fuel), based on Zhao et al. 2022
% "Representing vapour and capillary rise from the soil improves a leaf litter moisture model"

% for connection between adr (aerodynamic resistance) and gamma parameter,
% see commented lines after end of this function

% make sure one is on a row and one in a column
spechumids=spechumids(:)';
soilwcs=soilwcs(:);

% specific humidity of fuel (in the surrounding air) and saturated sx`pecific humidity at the soil
q_sats = spechumids(1);
q_fuel = spechumids(2:end); 

% volumetric soil water content, and value at field capacity
theta_fc = soilwcs(1);
theta    = soilwcs(2:end);

% soil_h the relative humidity of the air within the soil pore space where evaporation occurs (fixed at 1),
soil_h = 1;

% get air density
Pair = 101;    % air pressure, kPa
Rgas = 8.314;  % gas constant [J/mol/K]
m_air= 28.8;   % molecular weight of air [g/mol]
C2K  = 273.15; % conversion from Celcius to Kelvin
rho_a= Pair./((Rgas/m_air).*(TA+C2K)); % air density [kg/m^3]

% calculate beta: scale the effect of soil resistance on the vapour flux depending on soil wetness
beta = ones(size(theta)); % start with ones
nonfc= theta<theta_fc;    % when theta<theta_fc, calculate:
beta(nonfc) = 0.25*(1-cos(pi*theta(nonfc)./theta_fc)).^2; % eq.2 in Zhao2022

% replicate vectors
beta=repmat(beta,1,size(q_fuel,2));
q_fuel=repmat(q_fuel,size(beta,1),1);

% calculte vapor flux
VF = -rho_a.*(beta./adr).*(q_fuel-soil_h.*q_sats);
VF = VF*3600; % give the vapor flux in units of hours rather than seconds

end

% connection between gamma and r_a
% fbd=0.05; % fuel bed depth [m] -- this is consistent with Schimmel1997 Fig.5
% r_a = fbd/(24.5*10^-6); % (r_a = 1/K_vc vapor conductance), see eq.3 in Matthews 2005
% gamma=20*fbd; % areal density of fuel, assuming density of 20kg per m^3
% this result is consistent with our value of gamma=1 [kg/m^2]


