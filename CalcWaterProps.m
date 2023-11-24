function [prmsoil,prmfuel,prmair]=CalcWaterProps(TA,RH,P0,mm)
% by Yuval Zelnik (24/11/2023)
% Calculate water-related properties, based on air temperature (TA),
% relative humidity (RA), and soil water potential (P0) 
% TA is given in degrees Celcius [C], RH is a unitless fraction, P0 is given in [MPa] with sign
% OUTPUTS:
% prmsoil=[theta theta_fc] contains two soil volumetric water contents [m3 water/m3 soil] 
% prmfuel=[emc frh] contains the fuel moisture [g water/g dry-biomass] and relative humidity of air around it 
% prmair=[e_sat e_a q_sat q_a] is a vector air parameters, of vapor pressure and specific humidity 
% where vapor pressure (e_sat,e_a) is given in [kPa] and specific humidity (q_sat,q_a) in [g_H2O/g_air]

% Get theta, based on Campbell (1985) relationship:
% soil volumetric water content (theta, m3 water/m3 soil) and soil water potential (P0, with sign, MPa) 
psis=P0;
% parameters for sandy loam (taken from table 2 in Clapp&Hornberger 1978)
psissat=-0.7/10^3; %soil water potential at saturation - WITH SIGN -, MPa
bsoil=4.90; %exponent of the soil water retention curve
nsp=0.43; %soil porosity

% soil water retention curve
theta=nsp*(psis/psissat).^(-1/bsoil);

% same thing, now with values for field capacity
psifc=-0.033; % from Rai et al. 2017, chapter 17 page 513.
theta_fc=nsp*(psifc/psissat).^(-1/bsoil);

prmsoil=[theta_fc; theta(:)]; % soil volumetric water content, itself an at field capacity


% Get equilibirum moisture content, based on a semi-empirical function (Nelson Jr. 1984)
% General constants 
Rc = 1.987; % universal gas constant [cal/mol*K]
Mw = 18.015;  % molar mass of water [g/mol]
C2K= 273.15;    % conversion from Celcius to Kelvin
% Shape parameters for log function, based on data from literature (Nelson Jr. 1984, Catchpole et al. 2001)
a = 0.31;
b = -0.062; 

emc = a + b*log((-Rc/Mw)*(TA+C2K).*log(RH)); % calculate equilibirum moisture content


% we calculate the humidity in the air within the fuel (i.e. in the region between soil and top of fuel) 
% shape parameters for exp function, based on values from Matthews (2006)
A = 5.2;
B = -19; 
% Using an inverted form of Nelson Jr. (1984) formula, following eq.16 in Matthews (2006)
% Note the moisture input (mm) is the actual fuel moisture, and not emc calculated above
frh = exp(-(Mw./(Rc*(TA+C2K))).*exp(A+mm*B));
prmfuel=[emc; frh(:)]; % equilibirum moisture content and fuel relative humidity

% Get vapor pressure and specific humidity
% relationship air temperature (deg C) and relative humidity (RH, fraction) and saturated vapor pressure adnd vapor pressure deficit

%inputs: RH (), relative humidity; TA (deg C), air temperature
Pair  = 101;            % air pressure, kPa
m_air = 28.8;           % molecular weight of air [g/mol]
m_w   = 18;             % molecular weight of water [g/mol]
epsilon = m_w/m_air;    % ratio of molecular weights=0.622 [-]
% constants for the partial pressure of water at saturation esa = a * exp((b.*T)./(T + c)) [kPa, deg C]
% Based on Campbell&Norman (2000), eq. 3.8, also see Buck (1981)
a = 0.611;              % [kPa]
b = 17.502;             % [unitless]
c = 240.97;             % [deg C]

% vapor pressure 
e_sat=a.*exp((b.*TA)./(TA+c)); % air vapor pressure at saturation [kPa, with temp in deg C]
e_a=e_sat*RH; % partial vapor pressure [kPa]

% specific humidity
q_sat=e_sat/Pair*epsilon; % saturated. [g_H2O/g_air]
q_a=e_a/Pair*epsilon; % [g_H2O/g_air]


prmair=[e_sat e_a q_sat q_a]; % make a vector of the vapor pressures and specific humidities


end
