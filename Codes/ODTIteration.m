function Reconimg= ODTIteration(Reconimg,ORytov,TP,IterNum)
ORytov_index=(find(abs(ORytov)>0));
ORytov=ORytov(ORytov_index);
normFact=1./(TP.res3^2*TP.res4);
for mm = 1:IterNum
    if isempty(TP.signPhase)||(TP.signPhase==1)
        id = (real(Reconimg)<TP.n_m);
    else
        id = (real(Reconimg)>TP.n_m);
    end
    %     subplot(223),imagesc(real(squeeze(Reconimg(:,:,1))),[1.33 1.35]),axis image
    Reconimg(id)=TP.n_m-1i*imag(Reconimg(id));
    %     subplot(224),imagesc(real(squeeze(Reconimg(:,:,1))),[1.33 1.35]),axis image
    
    clear id;
    Reconimg=-(2*pi*TP.n_m/TP.lambda)^2.*(Reconimg.^2/TP.n_m^2-1);
    ORytov_new=fftshift(fftn(Reconimg))/normFact;
    %         ORytov_new=ORytov_new.*ORytov_index+ORytov;
    %     subplot(221),imagesc(log10(abs(squeeze(ORytov_new(:,:,end/2))))),axis image
    ORytov_new(ORytov_index)=ORytov;
    %     subplot(222),imagesc(log10(abs(squeeze(ORytov_new(:,:,end/2))))),axis image
    Reconimg=(ifftn(ifftshift(ORytov_new)))*normFact;
    Reconimg=TP.n_m*sqrt(1-Reconimg.*(TP.lambda/(TP.n_m*2*pi))^2);
end
Reconimg=(gather(fftshift(Reconimg,3)));
return