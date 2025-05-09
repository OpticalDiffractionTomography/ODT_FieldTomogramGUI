function [goodp2,coefficients,compensationMap]=phaseCompensation(varargin)
p2=varargin{1};
n=varargin{2};

if length(varargin)==3
    mask=varargin{3};
else
    mask=ones(size(p2));
end

[imY, imX]=size(p2);
[XX, YY]=meshgrid(1:imX,1:imY);

p2mask=p2.*mask;
%                 p2mask=omit_outliers(p2mask);
p2mask=p2mask(:);

X=zeros(sum(p2mask~=0),n);Y=X;
for ii=1:n
    XXX=XX.^ii;XXX=XXX(:);XXX(p2mask==0)=[];X(:,ii)=XXX;
    YYY=YY.^ii;YYY=YYY(:);YYY(p2mask==0)=[];Y(:,ii)=YYY;
end
p2mask(p2mask==0)=[];
E=ones(sum(p2mask~=0),1);AA=[X,Y,E];
coefficients=(AA'*AA)\(AA'*p2mask);

compensationMap=coefficients(end).*ones(imY,imX);
%             goodp2=p2-coefficients(end).*ones(imY,imX);
for ii=1:n
    compensationMap=compensationMap+coefficients(ii).*XX.^ii;
    compensationMap=compensationMap+coefficients(n+ii).*YY.^ii;
end
%             subplot(2,3,1),imagesc(p2),axis image
%             subplot(2,3,2),imagesc(p2-compensationMap),axis image
%             subplot(2,3,3),imagesc(p2.*mask),axis image
goodp2=p2-compensationMap;