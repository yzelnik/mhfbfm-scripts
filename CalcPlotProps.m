function [hc,lai]=CalcPlotProps(sa,sd)
% by Yuval Zelnik (24/11/2023)
% Calculate plot properties: canopy height (hc) and leaf area index (LAI)
% This is based on empirical relationships with stand age (sa)
% (with LAI being normalized by stand density (sd))
% sa is given in years, sd in stems per hectare
% hc is given in meters, lai in [leaf m^2/ ground m^2]

% calculate canopy height from empirical relationship with stand age
hcfunc=@(x,b) (b(1)*(1-exp(-x/b(2))).^b(3));
hcprms = [22 25 2]; %function is: 22*(1-exp(-sa/25))^2
hc = hcfunc(sa,hcprms); % canopy height

% calculate LAI using some intermediaries 

% calculate basal area (per tree) from empirical relationship with stand age
bafunc=@(x,b) (b(1)*(1-exp(-x/b(2))).^b(3));
baprms = [0.05 30 3]; %function is: 0.05*(1-exp(-sa/30)).^3
bapt = bafunc(sa,baprms);  % basal area per tree

% calculate dbh per tree (assuming round trees)
dbhpt=100*sqrt(bapt/pi)*2;  % dbh given in cm, therefore factor of 100

% calculate leaf area per tree 
% based on amended relationships from the literature: Marklund (1988) & Majasalmiet al. (2013)
lapt=5.5*exp(7.8*(dbhpt./(dbhpt+10))-3);

% calculate total leaf area (using stem density) [leaf m^2 / ground m^2]
lai=1e-4*lapt.*sd; % density is given per hectare, therefore factor of 0.0001

end
