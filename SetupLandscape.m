function matt=SetupLandscape(segfunc,szs,randkey,landprms,tradeofffunc,plotflag)
% by Yuval Zelnik (24/11/2023)
% Build landscape, including stand age&density, as well as soil-water-potential (SWP) and ,if relevant, altitude 
% INPUTS:
% segfunc is a function handle for a function that define stand shape and spatial organization
% szs is a vector of 3 length variables, the first 2 for the landscape, the third for average stand size (in units of length)
% randkey determines the randomization
% landprms is a vector of 9 parameters for defining the landscape (see below)
% tradeofffunc defines the relation between stand age and density that is assumed
% plotflag>0 will plot a visualization of the landscape
% OUTPUTS:
% matt is a 3-d matrix of the landscape, with the third dimension as: [age density SWP elevation]
if(nargin<6) plotflag=0; end;

standarrprms = landprms(1:3); % first 3 values are for proportion baseline value, and distribution width for area and proportion
moreprms = landprms(4:end); % the rest of the landscape parameters

if(length(tradeofffunc)<4) tradeofffunc(4)=0; end; % fraction of different species (e.g. deciduous)
if(length(moreprms)<5) moreprms(5)=0; end; % default of angle=0 for soil-water gradient
if(length(moreprms)<6) moreprms(6)=0; end; % default of no altitude difference

denstwdth=moreprms(1); % width of stand density randomzation
agemax   =moreprms(2); % max stand age
swpwdth  =moreprms(3); % width of soil water potential gradient
swpbase  =moreprms(4); % avg level of soil water potential
swpangl  =moreprms(5); % angle of gradient in degrees (-1 for random)
altdiff  =moreprms(6); % difference in altitude between high and low points (def=0)

% vectors of soil-water potential and altitude 
vecswp = [0.5:szs(2)/2 szs(2)/2-0.5:-1:0]*swpwdth*2/szs(2)+(swpbase-swpwdth/2);
vecalt = altdiff-altdiff*[0.5:szs(2)/2 szs(2)/2-0.5:-1:0]/(szs(2)*0.5);

if(tradeofffunc(1)==0) % power-law
    tof=@(x,prms) prms(2)*x.^prms(3);     
elseif(tradeofffunc(1)==1) % linear
    tof=@(x,prms) prms(2)+x*prms(3); 
else
    error('undefined function type');
end;

segmat=segfunc(szs,randkey,standarrprms);
unqvals=unique(segmat(:));

rng(randkey+1) % randomizing again, from somewhere else
% get stand ages from a uniform distribution
standages=randi(agemax,length(unqvals),1);
% get stand density using the tradeoff function multipled by a random factor
standdnss=tof(standages,tradeofffunc(1:3)).*(2.^(2*denstwdth*(rand(length(unqvals),1)-0.5)));

if(tradeofffunc(4)>0) % do we have a second tree species/type?
    typevec  = ceil(rand(length(unqvals),1)-tradeofffunc(4))*2-1; % fraction of landprms(4) is turned into -1
    standages= standages.*typevec;
end;

% stand age and density
matta=zeros(szs(1:2));
mattd=zeros(szs(1:2));

for ii=1:length(unqvals)   % for each segment (forest stand) 
    matta(segmat==unqvals(ii))=standages(ii);
    mattd(segmat==unqvals(ii))=standdnss(ii);   
end;

% replicate vectors of soil-water potential and altitude
mattswp=repmat(vecswp,szs(1),1);
mattalt=repmat(vecalt,szs(1),1);
if(swpangl) % if angle (on x-y plane) is not 0
    if(swpangl<0) % random angle
        rng(randkey+2) % randomizing again, from somewhere else
        swpangl=rand(1)*180;
    end;
    % rotate soil-water-potential
    tmpimg=ImgRotate([mattswp;mattswp],swpangl);
    mattswp=tmpimg(szs(1)/2+(1:szs(1)),:);
    mattswp(mattswp==0)=min(vecswp);
    
    % rotate altitude
    tmpimg=ImgRotate([mattalt;mattalt],swpangl);
    mattalt=tmpimg(szs(1)/2+(1:szs(1)),:);
    mattalt(mattalt==0)=max(vecalt);
    
end;
% add them all up together
matt=cat(3,matta,mattd,mattswp,mattalt);


if(plotflag==1) % plot things out?
  imagesc(cat(3,matt(:,:,1)/80,matt(:,:,2)/2500,(matt(:,:,3)+1.6)*0.9)-0.5);
elseif(plotflag==2)
  imagesc(matt(:,:,1));
end;

end