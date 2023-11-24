function segmat=SegmentAggregated(szs,randkey,moreprms)
% by Yuval Zelnik (24/11/2023)
% segment landscape using an aggregation algorithm
if(nargin<3) moreprms=[]; end; 

% Put default values into moreprms = [baseprop,varb-size,varb-prop,patchext]=[1,0,0,0]
if(length(moreprms)<1)   baseprop=1; else   baseprop=moreprms(1); end;
if(length(moreprms)<2)   varb(1)=0;  else   varb(1)=moreprms(2);  end;
if(length(moreprms)<3)   varb(2)=0;  else   varb(2)=moreprms(3);  end;
if(length(moreprms)<4)   patchext=0; else   patchext=moreprms(4); end;

% landscape size
basesz=szs(1:2);

%number of patches and their size
patchnum   = ceil(prod(basesz)/(szs(3)^2));
if(patchext)
    patchavgarea=patchext^2;
else
    patchavgarea=szs(3)^2;
end;

rng(randkey); % randomize
% define patches: corners, areas, proportions
pacorners = [randi(szs(1),patchnum,1) randi(szs(2),patchnum,1)];
patcharea = patchavgarea.*(1+2*varb(1)*(rand(patchnum,1)-0.5));
patchprop = baseprop.*(2.^(2*varb(2)*(rand(patchnum,1)-0.5)));

% use these to define the dimensions (width, length) of the patches
patchdims = ceil(repmat(sqrt(patcharea),1,2).*[patchprop 1./patchprop]);
szmax=max(patchdims(:)); % maximal dimensions of patches

segmat=zeros(szs(1:2)+szmax); % define the matrix, with padding
% run through patches/segments (skip over the first one, it is reserved for the baseline)
for ii=2:patchnum  
    % location and size of current patch/segment
    tmploc=pacorners(ii,:);   
    tmpsz=patchdims(ii,:);   
    % put value into this patch/segment
    segmat(tmploc(1)+(1:tmpsz(1)),tmploc(2)+(1:tmpsz(2)))=ii;
end;
segmat=segmat(2:basesz(1)+1,2:basesz(2)+1,:); % make sure we have the right size in the end

end