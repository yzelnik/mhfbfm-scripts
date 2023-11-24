function segmat=SegmentParceled(szs,randkey,moreprms)
% by Yuval Zelnik (24/11/2023)
% segment landscape into standard size parcels
if(nargin<3) moreprms=1; end; % proportions of x and y sizes

basesz=szs(1:2);
% define the parcel parameters: area, length-dimensions, number 
parcelarea = szs(3)^2;
parcelsz   = ceil(sqrt(parcelarea)*[moreprms(1) 1./moreprms(1)]);
parcelnum  = ceil(basesz./parcelsz);

rng(randkey); % randomize

[~,sortind]=sort(rand(prod(parcelnum),1)); % choose random order

% go through each parcel/segment, running on x and y axes
segmat=zeros(szs(1:2));
for ii=1:parcelnum(1) 
  for jj=1:parcelnum(2)
    % define location of parcel
    tmploc=[(ii-1)*parcelsz(1) (jj-1)*parcelsz(2)];
    % put in value for this parcel/segment
    segmat(tmploc(1)+(1:parcelsz(1)),tmploc(2)+(1:parcelsz(2)))=sortind(ii+(jj-1)*parcelnum(1));
  end;
end;
% make sure we have the correct size (cut out extra bits, if neccessary)
segmat=segmat(1:basesz(1),1:basesz(2),:);

end