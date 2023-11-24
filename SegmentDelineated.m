function segmat=SegmentDelineated(szs,randkey,moreprms)
% by Yuval Zelnik (24/11/2023)
% segment landscape using an deliniated algorithm
if(nargin<3) moreprms=[]; end;

maxiter=1e4;

% Put default values into moreprms = [baseprop,varb-size,varb-prop]=[1,0,0,0]
if(length(moreprms)<1)   baseprop=1; else   baseprop=moreprms(1); end;
if(length(moreprms)<2)   varb(1)=0;  else   varb(1)=moreprms(2);  end;
if(length(moreprms)<3)   varb(2)=0;  else   varb(2)=moreprms(3);  end;

% average area and number of segments (given the average area)
segavgarea = szs(3)^2;
segnum  = ceil(szs(1)*szs(2)/segavgarea);

rng(randkey); % randomize
% define centers and actual area for segments
cents = [randi(szs(1),segnum,1) randi(szs(2),segnum,1)];
segarea = segavgarea.*(1+2*varb(1)*(rand(segnum,1)-0.5));
% fix-up of segment area
segarea = round(segarea.*(prod(szs(1:2)))/sum(segarea));
toosmall=prod(szs(1:2))-sum(segarea);
segarea(1:min(toosmall,segnum))=segarea(1:min(toosmall,segnum))+1;

% define proportions of sements (width vs length)
segprop = baseprop.*(2.^(2*varb(2)*(rand(segnum,1)-0.5)));
% and now define (initial) dimensions of segments
segdims = ceil(repmat(sqrt(segarea),1,2).*[segprop 1./segprop]);

% calculate bounds for each segment (making sure we're not below 1, and not above landscape dimensions
segmins = max(cents+1-floor(segdims/2),1);
segmaxs = cents+ceil(segdims/2);
segmaxs(:,1)=min(segmaxs(:,1),szs(1)); 
segmaxs(:,2)=min(segmaxs(:,2),szs(2));

% go through segments
segmat=zeros(szs(1:2));
for ii=1:size(cents,1)
    ind=1;
    pixleft=segarea(ii);
    % repeatedly enlrage region of interest, until enough empty cells are found
    while(pixleft>0 && ind<maxiter)
        % how many empty cells do we have?
        tmppix = find(segmat(segmins(ii,1):segmaxs(ii,1),segmins(ii,2):segmaxs(ii,2))==0);
        tmpdim = segmaxs(ii,1)-segmins(ii,1)+1; % should be just segdims(1), but due to boundry issues
        
        newpix = segmins(ii,1)+mod(tmppix-1,tmpdim) + (floor((tmppix-1)/tmpdim)+segmins(ii,2)-1)*szs(1);
       
        if(length(newpix)>pixleft)
            % we got more cells than we need. just some, and finish up
            segmat(newpix(1:pixleft))=ii;
            pixleft=0;
        else
            % we either got exactly what we needed, or not enough yet. 
            segmat(newpix)=ii; % in any case, mark the ones we found
            pixleft=pixleft-length(newpix);     % update the number of cells we still need
            % the next 2 lines are relevant if we still didn't find enough cells for our segment
            % if we found exactly, they do not matter anymore as we're out of the loop just after this.
            segmins(ii,:)=max(segmins(ii,:)-1,[1 1]);
            segmaxs(ii,:)=min(segmaxs(ii,:)+1,szs(1:2));
        end;
        ind=ind+1;
    end;
end;

end