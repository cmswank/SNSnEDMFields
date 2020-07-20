%i is the number of winds per 1/2 coil. 

[bxx,Bx,By,Bz,bx,model]=fullscaleFourierB03D(ang,Bmsr,ferroNum)
for i=7:13


%ang is the initial otpimization parameters. For this case it is wire
%positions. 
ang=AngleFinalSD;

%mimization options.
options = optimset('MaxIter',1000,'MaxFunEvals',1000,'TolFun',1E-10,'TolX',1E-4);

%minimize is the function who's inputs are solely the parameters to be
%optmizied. 

minimize=@(angle)fullscaleFourierB03D(angle,1,1); %,0.001);

%x0 is inital parameters (copy of ang used to iput to fminsearch).  
x0=ang;
[AngleFinalSD,BxxFinalSD,BxxexitflagSD,BxxGOFSD]=fminsearch(minimize,x0,options);

                                                 %use fminunc for migrad
                                                 %use fmincon for
                                                 %constrained optimization
end



