function [Reconimg, ORytov]= ODTReconstruction_Berlin(retAmplitude,retPhase,TP)
ORytov=gpuArray(single(zeros(TP.ZP2,TP.ZP2,TP.ZP3)));
Count=(single(zeros(TP.ZP2,TP.ZP2,TP.ZP3)));
for kk=TP.frameList
    %     Hfilter=ones(TP.ZP,TP.ZP);
    Hfilter=mk_ellipse(8,6,TP.ZP,TP.ZP);
    Hfilter=circshift(Hfilter,[44 -8]).*circshift(Hfilter,-[42 -10]);
    %     subplot(221),imagesc(retPhase(:,:,kk),[-3 3]),axis image,colormap('jet')
    FRytov=squeeze(log(retAmplitude(:,:,kk))+1i*retPhase(:,:,kk));
    FRytov=gpuArray(padarray(FRytov,[round((TP.ZP-TP.xx)/2) round((TP.ZP-TP.yy)/2)],'replicate'));
    UsRytov=fftshift(fft2(FRytov)).*(TP.res2)^2.*Hfilter;   
    UsRytov=circshift(UsRytov,[round(TP.f_dx2(kk)) round(TP.f_dy2(kk))]);
     %     subplot(222),imagesc(log10(abs(UsRytov))),axis image
    xr=(TP.ZP*TP.res2*TP.NA/TP.lambda);
    UsRytov=UsRytov.*~mk_ellipse(xr,xr,TP.ZP,TP.ZP);
    [ky kx]=meshgrid(TP.kres*(-floor(TP.ZP/2)+1:floor(TP.ZP/2)),TP.kres*(-floor(TP.ZP/2)+1:floor(TP.ZP/2)));
    kz=real(sqrt((TP.n_m*TP.k0)^2-kx.^2-ky.^2));
    Kx=kx-TP.k0_x(kk);Ky=ky-TP.k0_y(kk);Kz=kz-TP.k0_z(kk);
    Uprime=1i.*2*2*pi*kz.*UsRytov;
    xind=find((kz>0).*~mk_ellipse(xr,xr,TP.ZP,TP.ZP)...
        .*(Kx>(TP.kres*(-floor(TP.ZP2/2)+1)))...
        .*(Ky>(TP.kres*(-floor(TP.ZP2/2)+1)))...
        .*(Kz>(TP.kres*(-floor(TP.ZP3/2)+1)))...
        .*(Kx<(TP.kres*(floor(TP.ZP2/2))))...
        .*(Ky<(TP.kres*(floor(TP.ZP2/2))))...
        .*(Kz<(TP.kres*(floor(TP.ZP3/2)))));
    
    Uprime=Uprime(xind);
    Kx=Kx(xind);
    Ky=Ky(xind);
    Kz=Kz(xind);
    
    Kx=round(Kx/TP.kres+TP.ZP2/2);Ky=round(Ky/TP.kres+TP.ZP2/2);Kz=round(Kz/TP.kres+TP.ZP3/2);
    KTP.ZP=(Kz-1)*TP.ZP2^2+(Ky-1)*TP.ZP2+Kx;
    
    [trash, idx]=unique(KTP.ZP,'first');
    KTP.ZP=(KTP.ZP(sort(idx)));    Uprime=Uprime(sort(idx));
    temp=ORytov(KTP.ZP);
    ORytov(KTP.ZP)=temp+Uprime;
    %     subplot(223),imagesc(log10(abs(ORytov(:,:,end/2)))),axis image
    %     subplot(224),imagesc(squeeze(log10(abs(ORytov(:,end/2,:))))),axis image
    
    Count(KTP.ZP)=Count(KTP.ZP)+1;
    %     pause()
end
% clearvars -except ORytov Count TP.ZP2 TP.n_m TP.lambda res3
ORytov=gather(ORytov);
ORytov(Count>0)=ORytov(Count>0)./Count(Count>0);
ORytov=gpuArray(ORytov);
%%
Reconimg=ifftn(ifftshift(ORytov))./(TP.res3^2*TP.res4);
Reconimg=TP.n_m*sqrt(1-Reconimg.*(TP.lambda/(TP.n_m*2*pi))^2);
Reconimg=(gather(Reconimg));
return;