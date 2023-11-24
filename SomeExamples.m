% this file is meant to illustrate the way the simulations work and their output
% the different sections are meant to be run separately, 
% with at least the first section to be run first
%% general setup

funcind=2; % which landscpae function: 1 for aggredated, 2 for delineated, 3 for parceled
landsize  = [250 250]; % size of landscape in cells

tradeoff = [0 2.5*10^4 -0.75]; % trade-off function between stand density and stand age
landfuncs = {@SegmentAggregated,@SegmentDelineated,@SegmentParceled}; % functions for creating landscape (stand arrangement)

rhoa = @(x,prms) tanh(x/prms(1));               % function of stand age [yrs]
rhod = @(x,prms) tanh(x/prms(1));               % function of stand density [#/ha]
alpw = @(ws,wa,prms) exp((ws/prms(1))*cos(wa)); % function of wind speed [m/s] and angle [rad]
alps = @(x,prms) exp(tanh(x/prms(1)));          % function of slope angle (rad)
rhom = @(fmd,fml,prms) (prms(2)./(prms(2)+(fmd+fml*prms(1))./(1+prms(1)))).^prms(3);     % function of fuel moisture content, both live and dead 

standarrbp  = 1;    % base level of proportion for x-to-y stand size
staddarrwa  = 0.5;  % distribution size for stand area
standarrwp  = 0.5;  % distribution size for stand proportions
densitywdth = 1;    % width of "noise" around stand density
agemax      = 99;   % 100 year range of stang age [yr]
moistbase   = -0.65; % avg value of soil-water potential [MPa]
moistwdth   = 1.2;  % max-min of values of soil-water potential  [MPa]
altdiff     = 0;   % difference between low and high points in landscape [m]
landprms = [[standarrbp staddarrwa standarrwp] densitywdth agemax moistwdth moistbase altdiff];

baseprms   = [40 1000 4 0.5 0.2 1.2]; % paramertes for fire spread: [rho_a rho_d c_w c_s k z]
extdfmprms = [2 1 2000 100 300 12 0.22]; % parameters for fuel-moisture: [tauh gamma adr windz qnet tday gsmww]
sigma      = 0.05; % sigma = proportion between dead and live fuel moisture in determining fire spread

p01=[0.5 1.0]; % p0 and p1 (baseline values of fire sprad probability)

pix2area=(20/100).^2; % converting pixel size to area in hectares

%% a simple plot of fire spread

parcelsize = 40; % average length (1-dimension) of a square stand [pixels]
topws= 2; % wind speed above canopy [m/s]
tarh = [20 0.1]; % temperature and relative humidity; [C] and []

% make a landscape, including stands and soil-water-potential
rndkey=1;
matt=SetupLandscape(landfuncs{funcind},[landsize parcelsize],rndkey,landprms,tradeoff);

% get wind-speed (in canopy) and moisture values (lfm, dfm)
[actws,dfm,lfm]=GetPreFireState(matt(:,:,1:2),matt(:,:,3),topws,tarh,extdfmprms);
% put this and more into a structure
st=struct('fa',rhoa,'fd',rhod,'fw',alpw,'fs',alps,'fm',rhom,'prma',baseprms(1),'prmd',baseprms(2),'prmw',baseprms(3),'prms',baseprms(4),'prmm',[sigma baseprms(5:6)],'p0',p01(1),'p1',p01(2),'slp',0,'ws',actws,'wa',0,'fmd',dfm,'fml',lfm);
   
% run a single fire simulation
rng(1);
[outmat,burntime]=RunFire(matt,[100;1],st);
imagesc(outmat)

%% plot different landscape aspects and fire spread scenarios


parcelsize = 50;  % define the stand size
topws= 8;         % wind-speed above canopy
tarh = [20 0.2];  % temperature and relative-humidity (RH)

% some numbers for plotting
plotfactors=[100 -1.5];
plotcolorg =[0. 0 0.1;  0. 0.3 0.9];
plotcoldst =[0.3 1 0.1; 1 1 1];


rndkey=5; % choose randomization
% build the landscape
matt=SetupLandscape(landfuncs{funcind},[landsize parcelsize],rndkey,landprms,tradeoff);

% show the plot age, and the soil-water-potential
for ii=1:2
  subplot(2,2,ii)
  tmp=matt(:,:,ii*2-1)/plotfactors(ii);
  imagesc(repmat(reshape(plotcolorg(ii,:),[1 1 3]),landsize)+repmat(tmp,[1 1 3]).*repmat(reshape(plotcoldst(ii,:)-plotcolorg(ii,:),[1 1 3]),landsize))
end;

% simulate and show the results of two fires
for ii=1:2
  subplot(2,2,ii+2)
  % calculate conditions, and put it into structure
  [actws,dfm,lfm]=GetPreFireState(matt(:,:,1:2),matt(:,:,3),topws,tarh,extdfmprms);
  st=struct('fa',rhoa,'fd',rhod,'fw',alpw,'fs',alps,'fm',rhom,'prma',baseprms(1),'prmd',baseprms(2),'prmw',baseprms(3),'prms',baseprms(4),'prmm',[sigma baseprms(5:6)],'p0',p01(1),'p1',p01(2),'slp',0,'ws',actws,'wa',0,'fmd',dfm,'fml',lfm);
  % run fire (with 10 cells burnt as initial conditions)
  rng(2);
  [outmat,burntime]=RunFire(matt,[100:109;ones(1,10)],st);
  imagesc(outmat) % show fire result

  tarh=[20 0.3]; % update RH for the second scenario (RH=0.3)
end;


%% run multiple simulations to assess fire risk
% this may take quite a while to run
% plots are in next sections

% average stand size (in cell length)
parcelszs = [160 80 40 20 10 5 2];

landnum  = 50; % number of landscapes per parameter set 
simnum   = 100;% number of simulations per landscape
tarh     = [20 0.2]; % temperature and relative-humidity
topws    = 2;  % wind speed above canopy
startpnts= 10; % number of burning cells at the start of the fire

simstats=[]; firesizes=[]; % initilize
% go over landscape parameter -- stand size
for lndind=1:length(parcelszs)
  tic; 
  % go over randomization of landscapes, i.e. landscape index
  for rndkey=1:landnum
    rng(rndkey+lndind*1e3);
    % create landscape
    matt=SetupLandscape(landfuncs{funcind},[landsize parcelszs(lndind)],rndkey,landprms,tradeoff);

    % get wind-speed (in canopy) and moisture values (lfm, dfm); 
    [actws,dfm,lfm,extras]=GetPreFireState(matt(:,:,1:2),matt(:,:,3),topws,tarh,extdfmprms);
    % save these and more into a structure
    st=struct('fa',rhoa,'fd',rhod,'fw',alpw,'fs',alps,'fm',rhom,'prma',baseprms(1),'prmd',baseprms(2),'prmw',baseprms(3),'prms',baseprms(4),'prmm',[sigma baseprms(5:6)],'p0',p01(1),'p1',p01(2),'cellsz',20,'alt',matt(:,:,4),'ws',actws,'wa',0,'fmd',dfm,'fml',lfm);
      
    tmpstats=[];
    % go over fire simulations per landscape
    for ii=1:simnum
        rng(rndkey*1e3+lndind*1e6+ii) % randomize
        st.wa = (rand(1)-0.5)*pi/2; % choose random wind angle
        
        % choose random locations for start of fire on western edge of landscape
        [~,tmplocs]=sort(rand(1,landsize(1)*0.5,1)); 
        starty=round(landsize(1)*0.25+tmplocs(1:startpnts));
        
        rng(rndkey*1e3+lndind*1e6+ii);% randomize
        % run fire simulation
        outmat=RunFire(matt,[starty;ones(1,startpnts)],st);
        % save info on how many cells were burnt in the end
        tmpstats(ii)=nnz(outmat);
    end;
    % save the statistics of fire sizes
    firesizes(lndind,rndkey,:)=tmpstats;
    % save simulation statistics, of average and maximum
    simstats(lndind,rndkey,:)=[mean(tmpstats) max(tmpstats)];
  end;
 
  toc;
end;



%% plot out the results of fire risk statistics
cell2area=(20/100).^2; % 20[m]/100[m] or for area: 20^2/10000 [m^2] = cell-area/hectare 

% plot average fire risk
plot(log10(cell2area*parcelszs.^2),log10(cell2area*mean(simstats(:,:,1),2)),'k-o')
set(gca,'XTick',[-1:3],'XTickLabel',10.^[-1:3],'YTick',log10([10 20 50 100 200 500 1000]),'YTickLabel',[10 20 50 100 200 500 1000])

% add on to that the median and max
hold on;
plot(log10(cell2area*parcelszs.^2),log10(cell2area*median(simstats(:,:,1),2)),'b-*')
plot(log10(cell2area*parcelszs.^2),log10(cell2area*max(simstats(:,:,1),[],2)),'r-+')
hold off


%% plot out histograms of fire risk and fire size
cell2area=(20/100).^2; % 20[m]/100[m] or for area: 20^2/10000 [m^2] = cell-area/hectare 

% which stand size (index) do we want to look at?
landchs = 1;  % compare 1 and 6

subplot(1,2,1)
% histogram of fire risk
histogram(log10(cell2area*simstats(landchs,:,1)),10)
set(gca,'XTick',[-1:3],'XTickLabel',10.^[-1:3])
title('fire risk distribution','fontSize',20)


subplot(1,2,2)
% histogram of fire size
tmpfs=firesizes(landchs,:,:); % choose one of the stand size average
histogram(log10(cell2area*tmpfs),20)
set(gca,'XTick',[-1:3],'XTickLabel',10.^[-1:3])
title('fire size distribution','fontSize',20)
