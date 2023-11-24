function [actws,dfm,lfm,extras]=GetPreFireState(standstate,p0land,topws,tarh,extprms)
% by Yuval Zelnik (24/11/2023)
% get all the necessary (spatial) enviornmental parameters before a fire
% INPUTS:
% standstate is a 3-d matrix, with the third dimension for age and density
% p0land is similar matrix for the soil-water potential
% topws is the reference speed (at height windz=100 [m])
% tarh=[ta rh] with ta the temperature and rh the relative humidity of air
% extprms are additional parameters for dead and live fuel moisture
% OUTPUTS:
% actws [m] is the wind speed inside the canopy (at half canopy height)
% dfm [gH2O/gDW] is the dead fuel moisture
% lfm [gH2O/gDW] is the live fuel moisture
% extras is a cell structure of: extras={canopy-height,lead-area-index,evapotranspiration}

if nargin<5 extprms=[2 1 2000]; end;


% extprms = [tauh gamma adr windz qnet tday gsmww]. as default = [2 1 2000 100 300 12 0.22]
tauh  = extprms(1);  % [hour]   timescale of equilibiration of dead fuel 
gamma = extprms(2);  % [kg/m^2] (spatial ensity of dead fuel)
adr   = extprms(3);  % [s/m]    (aerodynamic resistance)
windz = extprms(4);  % [m]      (height above ground, assumed for measurements)
Q_net = extprms(5);  % [W/m^2]  (solar radiation)
t_day = extprms(6);  % [hour]   (length of day)
gsmww = extprms(7);  % [mol air/m^2 leaf/s]  (stomatal conductance, at well watered conditions)


relz  = 0.5; % relative height at which to assume fire spread (half of canopy height)
absz  = 0;   % additional height (not proportional to canopy height)

% stand age and density
sa=standstate(:,:,1);
sd=standstate(:,:,2);

% calculate tree height (hc, [m]) and leaf area (LAI [m^2/m^2])
[hc,lai]=CalcPlotProps(sa,sd);

TA = tarh(1); % temperature [C]
RH = tarh(2); % relative humidity []

% Numerically find  the dead-fuel moisture
mgrad=(1:500)/1e3;
% Find out what are the different values of P0 (this is more efficient this way)
p0vals = unique(p0land);
dfm = zeros(size(p0land));
% For each P0 value, numerically find (approximately) the value of VF that leads to equilibrium, thus giving us the value of dmc
for ii=1:length(p0vals)
  % Calculate water properties, including a gradiant of possible fuel moisture values
  [prmsoil,prmfuel,prmair]=CalcWaterProps(TA,RH,p0vals(ii),mgrad);
  
  % For this fuel moisture gradient, what are the values of vapor flux?
  VF=CalcVaporFlux(prmair(3).*[1;prmfuel(2:end)],prmsoil,adr,TA); 
  % Now equate the two parts of the equilibrium equation
  % and find the minimum --> numerically close to 0
  [~,minind]=min(abs((mgrad-prmfuel(1))/tauh-VF/gamma)); % find which value of VF will lead to equilibrium
  
  % translate this into the moisture value that gives equilibrium
  dmcvals(ii)=mgrad(minind);
  dfm(p0land==p0vals(ii))=dmcvals(ii);
end;
%VF=CalcVaporFlux(prmair(3).*[1;prmfuel(1+minind)],prmsoil,adr,TA); 

% Estimate evapotranspiration
ET=EvapoTransp(Q_net,TA,t_day,topws,windz,prmair(2),lai,hc,gsmww,p0land);
% With this, now estimate live fuel moisture (using tree water potential)
lfm=CalcLFM(ET,t_day,p0land,hc,lai);

% get the wind speed below canopy for fire spread
actws=CorrectWind(hc,windz,topws,hc*relz+absz);
% return more info [canopy height, leaf area index, evapotranspitation]
extras={hc,lai,ET};

end


