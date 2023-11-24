function [firemat,burntime]=RunFire(matt,startloc,st)
% by Yuval Zelnik (24/11/2023)
% run fire simulation, given the matrix of landscape features (matt), 
% starting location(s) of fire (startloc), and state-structure of other conditions (st)
% returns firemat (matrix of burnt cells) and burntime (iteration number overwhich fire has spread)

% size of landscape (matrix)
szs=[size(matt,1),size(matt,2)];
% define initial fire matrix
firemat=zeros(szs);
for ii=1:size(startloc,2)
   firemat(startloc(1,ii),startloc(2,ii))=1;
end;

% effective wind in 4 firections (1=east, 2=south, 3=west, 4=north)
for ii=1:4
  factorwindslope{ii}=st.fw(st.ws,(ii-1)*pi/2-st.wa,st.prmw);
end;


% Is there slope?
if(isfield(st,'cellsz') && size(matt,3)>3) 
  for ii=1:4
    if(ii==1) % calculate slope along x-axis
      slope=diff(matt(:,:,4)/st.cellsz,[],2);
      slope=[slope(:,1)  (slope(:,1:end-1)+slope(:,2:end))/2  slope(:,end)];
      tmpslope{ii}=slope;
    elseif(ii==2) % calculate slope along y-axis
      slope=diff(matt(:,:,4)/st.cellsz,[],1);
      slope=[slope(1,:) ; (slope(1:end-1,:)+slope(2:end,:))/2 ; slope(end,:)];
      tmpslope{ii}=slope;
    else % the other 2 directions just take as the negative of previous two
      slope=-tmpslope{ii-2}; 
    end;
   % size(slope)
    factorwindslope{ii}=factorwindslope{ii}.*st.fs(slope,st.prms);
  end;
end;
    
% Do we have a second cover type (e.g. non conifer)?
if(min(min(matt(:,:,1)))<0) % do we have a second type of trees (e.g. deciduous)?
    % NOTE: THIS IS A AD-HOC WORK AROUND, WITH NEGATIVE TREE AGE REPRESENTING SECOND TYPE OF TREES 
    tmpp1 = st.p1(1) + (st.p1(2)-st.p1(1))*(matt(:,:,1)<0); % p1 that works for two species
    % matrix of stand effect
    factorstandprop = st.p0+(tmpp1-st.p0).*st.fa(abs(matt(:,:,1)),st.prma).*st.fd(matt(:,:,2),st.prmd);
else % simpler case, only one type
  % matrix of stand effect
  factorstandprop = st.p0+(st.p1(1)-st.p0).*st.fa(matt(:,:,1),st.prma).*st.fd(matt(:,:,2),st.prmd);
end;

% matrix of fuel moisture effect
factormoisture = st.fm(st.fmd,st.fml,st.prmm);


% iterations/burntime
maxiter=1e6;
burntime=1;
% continue while the fire can spread further and we haven't gone over our max iterations
while(max(firemat(:))>0 && burntime<maxiter)
    firetmp=zeros(szs);
    
    % find the x,y locations of burnt spots
    [xx,yy]=ind2sub(szs,find(firemat>0));
    % find neighbors of burning cells (AUX FUNC)
    [tmplocs,tmpdirect]=neigbs([xx';yy'],szs,firemat);
    
    tmpfactorwindslope=zeros(size(tmpdirect));
    for ii=1:4 % run through 4 directions
        tmpinds=(tmpdirect==ii);
        tmpfactorwindslope(tmpinds)=factorwindslope{ii}(tmplocs(tmpinds)); 
    end;
    % calculate fire spread probability
    pp = (1-  (1-factorstandprop(tmplocs)).^tmpfactorwindslope) .*factormoisture(tmplocs);
    % random uniform distribution
    randvec=rand(size(pp));
   
    % enact burnt cells of this iteration
    firetmp(tmplocs(pp>randvec))=1;
    
    % update fire matrix
    firemat(firemat==1)=-1;
    firemat(firetmp==1)=1;
    
    burntime=burntime+1;
end;

if(burntime==maxiter) error('simulation did not finish its run properly!'); end;

end



%%% AUX FUNC %%% 
% find neighboring cells
function [inds,direct]=neigbs(foclocs,szs,mask)
  % find neighboring cells
  focalnum=size(foclocs,2);
  % use a stencil for the 4 neighboring cells around focal locations
  % stencil being: [0 1 0 -1; 1 0 -1 0]
  stencil = [zeros(1,focalnum) ones(1,focalnum) zeros(1,focalnum) -ones(1,focalnum) ; ones(1,focalnum) zeros(1,focalnum) -ones(1,focalnum) zeros(1,focalnum)];
  locs=repmat(foclocs,1,4)+stencil;
  direct=ceil((1:focalnum*4)/focalnum); % direction (values of 1 through 4)
  
  % ignore things out of bounds
  chosen=find((locs(1,:)>0 & locs(2,:)>0 & locs(1,:)<=szs(1) & locs(2,:)<=szs(2)));
  
  % change from (x,y) to single index
  inds=sub2ind(szs,locs(1,chosen),locs(2,chosen));
  % take out locs that were in the mask (i.e. already burnt)
  chosen(~~mask(inds))=[];
  
  % update inds (to avoid burnt cells as well)
  inds=sub2ind(szs,locs(1,chosen),locs(2,chosen));
 
  % now take only the chosen ones (inside bounds and not burnt)
  direct=direct(:,chosen);
end
