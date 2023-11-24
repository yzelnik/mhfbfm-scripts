function newimg=ImgRotate(orgimg,rotangle)
% use to rotate 2-d image

szs=size(orgimg);
loc=szs/2;

extrat=1.41;
tmpsz = ceil(extrat*szs/2)*2;

% add padding to avoid errors
if(loc(1)>size(orgimg,1)-tmpsz(1))
    orgimg(loc(1)+tmpsz(1),1,1)=0;
end;
if(loc(2)>size(orgimg,2)-tmpsz(2))
    orgimg(1,loc(2)+tmpsz(2),1)=0;
end;
if(loc(1)<tmpsz(1))
    orgimg=[zeros([tmpsz(1),size(orgimg,2),size(orgimg,3)]); orgimg];
    loc(1)=loc(1)+tmpsz(1);
end;
if(loc(2)<tmpsz(2))
    orgimg=[zeros([size(orgimg,1),tmpsz(2),size(orgimg,3)]) orgimg];
    loc(2)=loc(2)+tmpsz(2);
end;

imagepad = orgimg(loc(1)-tmpsz(1)/2+(1:tmpsz(1)),loc(2)-tmpsz(2)/2+(1:tmpsz(2)),:);

[nrows,ncols,~] = size(imagepad);
midx=ceil((ncols+1)/2);
midy=ceil((nrows+1)/2);

piangle = rotangle*pi/180;
Mr = [cos(piangle) sin(piangle); -sin(piangle) cos(piangle)]; 

% rotate about center
[X Y] = meshgrid(1:ncols,1:nrows);
XYt = [X(:)-midx Y(:)-midy]*Mr;
XYt = bsxfun(@plus,XYt,[midx midy]);

xout = round(XYt(:,1)); yout = round(XYt(:,2)); % nearest neighbor!
outbound = yout<1 | yout>nrows | xout<1 | xout>ncols;
if(size(orgimg,3)==3)
    zout=repmat(cat(3,1,2,3),nrows,ncols,1); 
    zout=zout(:);
else
    zout=repmat(1,nrows,ncols,1);
    zout=zout(:);
end;

xout(xout<1) = 1; xout(xout>ncols) = ncols;
yout(yout<1) = 1; yout(yout>nrows) = nrows;
xout = repmat(xout,[size(orgimg,3) 1]); yout = repmat(yout,[size(orgimg,3) 1]);
imagerot = imagepad(sub2ind(size(imagepad),yout,xout,zout(:))); % lookup
imagerot = reshape(imagerot,size(imagepad));
imagerot(repmat(outbound,[1 1 size(orgimg,3)])) = 0; % set background value to [0 0 0] (black)

newimg = imagerot((tmpsz(1)-szs(1))/2+(1:szs(1)),(tmpsz(2)-szs(2))/2+(1:szs(2)),:);

end
