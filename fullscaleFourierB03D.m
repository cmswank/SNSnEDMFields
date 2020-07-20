function [bxx,Bx,By,Bz,bx,model]=fullscaleFourierB03D(ang,Bmsr,ferroNum)
%function that determines the fourier weight (proportinal to relaxation
%rate)


if nargin < 3 
    ferroNum=1;
end
if nargin <2
    Bmsr=1;
end

%Fourier coordinates, (discrete transform fourier coefficient positions)
[coords] = fullscaleCoordinatesB0SYM3D();

%Find magnetic field from comsol given angle of wires
%also coordinates points

%#2 has the design upset. 
[model,ang,mesherror] = fullscaleComsol3DSYMFerroCloak2FerroOpt(ang,Bmsr,ferroNum);
%[model,ang,mesherror] = fullscaleComsol3D_romerXS(ang,Bmsr,ferroNum);

if mesherror==true
    ang=CosCoilGeneratorAlongX(length(ang));
    [model] = fullscaleComsol3DSYMFerroCloak2FerroOpt(ang,Bmsr,ferroNum);
end
info=mphgetselection(model.selection('geom1_csel2_bnd'));
figure(82)
mphviewselection(model,'geom1',info.entities,'boundary','facecolorselected',[1 1 0],'facealpha',0.5)
drawnow

%Make a goofy movie (Larry says it needs goofy audio)
%load('F2File.mat');
%F2(f2n)=getframe(figure(82));
%f2n=f2n+1;
%save('F2File','F2','f2n');

%  myVideo=VideoWriter('heatflushmovie.avi');
%  myVideo.FrameRate=2.5;
%  open(myVideo);
%  writeVideo(myVideo,F2);
%  close(myVideo);
%T2 and frequency shifts
%[Bxp,By,Bz] = mphinterp(model,{'mf.Bx','mf.By','mf.Bz'},'coord',coords');

%just T2. 
[Bxp,Byp,Bzp] = mphinterp(model,{'mf.Bx' 'mf.By' 'mf.Bz'},'coord',coords');

%create a matrix for Bx,By

%Bx=zeros(17);
Bz=zeros(13,7,7);
By=zeros(13,7,7);
Bx=zeros(13,7,7);
i4=1;
for iii = 1:49:637
Bxpp=zeros(7,7);
Bypp=zeros(7,7);
Bzpp=zeros(7,7);
ii=1;
for i = 1:7:49 %1x14 and 2x7 half coefficients, equivelent to 14x14x14 due to symmetry 
    Bxpp(ii,:)=Bxp((iii-1)+i:(i+6)+(iii-1));
    Bypp(ii,:)=Byp((iii-1)+i:(i+6)+(iii-1));
    Bzpp(ii,:)=Bzp((iii-1)+i:(i+6)+(iii-1));
ii=ii+1;
end
Bx(i4,:,:)=Bxpp;
By(i4,:,:)=Bypp;
Bz(i4,:,:)=Bzpp;
i4=i4+1;
end

%%%no symmetry in x actually. 
%%%Symmetry in y and z. 

%create full cell (take away symmetry (still symmetric obviously))

Bxcell=zeros(13,13,13);
Bycell=zeros(13,13,13);
Bzcell=zeros(13,13,13);
for i = 1:size(Bx,1)
Bxcell(i,:,:)=[squeeze(flip(flip(Bx(i,:,:),2),3)),squeeze(flip(Bx(i,:,2:end),2));squeeze(flip(Bx(i,2:end,:),3)),squeeze(Bx(i,2:end,2:end))];
Bycell(i,:,:)=[squeeze(flip(flip(By(i,:,:),2),3)),squeeze(flip(By(i,:,2:end),2));squeeze(flip(By(i,2:end,:),3)),squeeze(By(i,2:end,2:end))];
Bzcell(i,:,:)=[squeeze(flip(flip(Bz(i,:,:),2),3)),squeeze(flip(Bz(i,:,2:end),2));squeeze(flip(Bz(i,2:end,:),3)),squeeze(Bz(i,2:end,2:end))];
end

%By=[flip(flip(By,1),2),flip(By(:,2:end),1);flip(By(2:end,:),2),By(2:end,2:end)];

%create mirror images (not needed, theory is the same for the real part. ie
%T2, what we are minimizing... also unsymmetric contributions are
%identically zero and only introduced by fabrication defects. 
%Bx=[Bx,flip(Bx(:,1:end-1),2);flip(Bx(1:end-1,:),1),flip(flip(Bx(1:end-1,1:end-1),1),2)];
%By=[By,flip(By(:,1:end-1),2);flip(By(1:end-1,:),1),flip(flip(By(1:end-1,1:end-1),1),2)];


%this is a B0 coil in an MSR so we CANNOT normalize. 
Bxnorm=0.03/mean(mean(mean(Bxcell)));
disp(Bxnorm);
Bx=Bxcell*Bxnorm;
By=Bycell*Bxnorm;
Bz=Bzcell*Bxnorm;
%By=By*Bxnorm;
figure(83)
%load('F1File.mat');
%surf(Bx(:,:,1))

%curraxis=gca;
%if f1n==1

%    zaxis=curraxis.ZLim;
%end

%curraxis.ZLim=zaxis;

drawnow

% F1(f1n)=getframe(figure(83));
% f1n=f1n+1;
% save('F1File','F1','f1n','zaxis');
%compute fourier spectrum at q=n*Pi/L.

%T2 per component
bx=8*abs(fftn(Bx-mean(mean(mean(Bx))))/size(Bx,1)/size(Bx,2)/size(Bx,3));

%phase shift (ish)
%by=4*fftn(By-mean(mean(By)))/size(By,1)/size(By,2);

%remember its exp(-i*q(x-x0)) so we multiply by its conjugate to get the
%unweighted field correlation coefficients. (this is correct for us because
%of symmetry, the true way is to do is real(dct(Bx).*fft(Bx)) but our coils
% must be symmetric in x and y. 

%bx=real(bx);
%by=real(by);

%add conditional density factor for better optimization. 
% for i = 1:size(bx,1)
%     for ii = 1 :size(bx,2)
%     
%         if i < 2 && ii <2
%             relaxNorm(i,ii)=1;
%         else
%             relaxNorm(i,ii)=1/sqrt( (i<=6)*(i-1)^2+(i>6)*(13-i)^2+(ii<=6)*(ii-1)^2+(ii>6)*(13-ii)^2);
%         end
%        
%        
%         
%     end
%     
% end



%optimize for neutrons;
%13x13x13  
% relaxNormpp=linspace(-6,6,13);
% relaxNormp=meshgrid(relaxNormpp,relaxNormpp);

%bxx is what we wish to minimize 

%bx=bx.*relaxNorm;


%T2, if you want conditional probability factor include .*relaxNorm  (or not in this case). 
bxx=sum(sum(sum((bx.^2)))); %% 1E-7 G/cm (40 cm) is bxx=2.195E-11 (this number is for 2D)! (GOAL) this is 4.7E3 seconds T2



%without conditional probabilty factor. 
%bxx=sum(sum(sum((bx.^2)))); 

 bxxp=bxx;
 angp = ang;
 
%save best results per solution
if exist(['AngleFinalB0FullscaleUpset1',num2str(length(ang)),'.mat'],'file')
load(['AngleFinalB0FullscaleUpset1',num2str(length(ang)),'.mat']);
if bxxp<bxx
ang=angp;
bxx=bxxp;
save(['AngleFinalB0FullscaleUpset1',num2str(length(ang))],'ang','bxx');
else 
ang=angp;
bxx=bxxp;
end


else
    save(['AngleFinalB0FullscaleUpset1',num2str(length(ang))],'ang','bxx');
end

%disply results per run
%angs='';
%for i = 1:length(ang)
%angs=[angs,num2str(ang(i)), ', '];  
%end
%disp(['Angles = ',angs]);
disp(['bxx =  ',num2str(bxx),', Goal < 2.73E-8']);  
%Goal is equivielent of T2 for a linear gradient at 1.5 ppmB0/cm over 40 cm.(3D) 

