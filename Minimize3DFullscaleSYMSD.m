%Simplex Brute force search minimization
% 
% for i = 15:18
% ang=CosCoilGenerator(i);
% if i==15
%     load('PulstarB015.mat');
%     ang=AngleFinal;
% end
% 
% options = optimset('MaxIter',1000,'MaxFunEvals',1000,'TolFun',1E-3,'TolX',1E-5);
% minimize=@(angle)PulstarComsolB0FC3D(angle);
% x0=ang;
% [AngleFinal,BxxFinal,Bxxexitflag,BxxGOF]=fminsearch(minimize,x0,options);
% end

%[bxx,Bx,bx,model]=fullscaleFourierB03D(ang,Bmsr,ferroNum)

for i=25

%if i~=15
 ang=CosCoilGeneratorAlongX(i);
%end
 %[~,~,~,modelcos]=fullscaleFourierB03D(ang);
%save(['CosineCoilB0FullscaleUpset',num2str(length(ang))],'modelcos');
 %end
%ang=ang(2:end);
%for i = 20:35

%[ang] = CosCoilGeneratorAlongX(n);

options = optimset('MaxIter',1000,'MaxFunEvals',1000,'TolFun',1E-12,'TolX',1E-4);
minimize=@(angle)fullscaleFourierSD3D(angle);
%bxx16(i)=
x0=ang;
[AngleFinalSD,BxxFinalSD,BxxexitflagSD,BxxGOFSD]=fminsearch(minimize,x0,options);


end


