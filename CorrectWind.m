function u=CorrectWind(hc,zref,uref,z)
% Correct wind speed to values around or below canopy

%Giulia Vico February 2023
%based on Campbell and Norman, Ch 5

%input
%hc:   canopy height (m)
%zref: height of measurement (m)
%uref: wind speed at measurement (m/s)
%z:    vertical positions

%NOTE: if using also EvapoTransp, zref and uref here should match z and u there

%output: 
%u: wind speed (m/s)


%parameters
%--- above canopy parameters - %NOTE: if changing d and z0, change them also in EvapoTransp.m
d = 2/3*hc; % zero plane displacement height for dense canopies (m); the zero plane displacement is such that d+zm is the height at which the wind from the logaritmic profile extrapolates to zero
z0 = 0.1*hc; % momentum roughness height (m) - from Katul et al 2004: could be 0.08 for Pinus taeda and hardwood forest
%from Campbell and Norman, Table5.1: z0=110/100  for coniferous forest, 0.1/100 for grass, closely mowed, 9/100 for thick grass (>50 cm)
%--- within canopy, top 90%
a=1; %attenuation coefficient; 1 for larch, 1.1 for Xmas trees: can also be calculated given canopy properties
%--- within canopy, bottom 10%
dsurface=0;
z0surface=3/100;  %(m) momentum roughness height of the surface, assumed intermediate between closely mowed grass and thick grass and similar to alfalfa
vonkarmankappa=0.4;

extfactor = (1-2/3)/0.1; % just a short cut (see above) for: (hc-d)/z0

%setup:
isabove=find(z>=hc);
isbelow90=find(z<hc & z>0.1*hc);
isbelow10=find(z<=0.1*hc);

% lengths of 3 parts
lens=[length(isbelow10) length(isbelow90) length(isabove)];
u=zeros(size(hc)); % initilize


% above the canopy (eq 5.1 Campbell and Norman)
ustar=vonkarmankappa*uref./log((zref-d)./z0); %friction velocity, from inversion of eq 5.1 at zref (the friction velocity is constant at all heights)
if(lens(1)>0)
  u(isabove)=ustar(isabove)./vonkarmankappa.*log((z(isabove)-d(isabove))./z0(isabove));
end;

if(lens(2)>0)
  ucanopy = ustar(isbelow90)/vonkarmankappa.*log(extfactor);
  % upper 90% of the canopy (eq. 5.4 Campbell and Norman)
  % alternatively, it is possible to estimate that based on LAI
  u(isbelow90)=ucanopy.*exp(a.*(z(isbelow90)./hc(isbelow90)-1)); 
end;

if(lens(3)>0)
  ubetween = (ustar(isbelow10)/vonkarmankappa.*log(extfactor)) .*exp(a.*(0.1-1)); 
  % bottom 10% (Campbell and Norman page 69)
  ustar=vonkarmankappa*ubetween./log((z0(isbelow10)-dsurface)/z0surface); %friction velocity, from inversion of eq 5.1 at zref (the frixtion velocity is constant at all heights)
  u(isbelow10)=max(0,ustar./vonkarmankappa.*log((z(isbelow10)-dsurface)/z0surface));
end;

end
