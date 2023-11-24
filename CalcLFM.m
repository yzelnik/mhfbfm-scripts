function lfm=CalcLFM(ET,tday,p0,hc,lai)
% Calculate live fuel moisture
% by Yuval Zelnik (24/11/2023)
% largely based around the TreeWaterPotential function by Giulia Vico

% Paraeters used to estimate water potential
a  = 1.07;            % 1.07 in Table 3 Couvrer; 1.6: value for Pinus Ponderosa in Couvrer
bTs= [-5.79,-0.0192]; % Couvrer Table 3, does not match exactly the value for either Psudotsuga or Pinus in Table 2 though
bs = 0.08         ;   % MPa/m: average value for conifer in Couvrer
h0s= [4.3,-0.050];    % dependence of Huber value on L; Couvrer Table 3; return Huber in c2sap/m2leaf
k0s= [1.4,0.11];      % dependece of saturated hydraulic conductivity at the base of the tree (z = 0) on L; Couvrer, Table 3; for Pinus, [1.1,0.11]

maxval=1.3; % max moisture value [gH2O/gDW]

% get water potential
P=TreeWaterPotential(ET,tday,p0,hc,lai, k0s,h0s,a,bTs,bs);

afitInv = [-7.924 -0.568]; % Fitting parameters for Picea needle 

% get Relative-Water-Content from the water potential (using fitting parameters)
RWC=(afitInv(1)-P.*(1+afitInv(2)))./(afitInv(1)-afitInv(2).*P);
lfm = maxval*RWC; % multiply by maximum value of live fuel moisture

end