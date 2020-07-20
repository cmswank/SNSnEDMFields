function [coords] = fullscaleCoordinatesB0SYM3D()
%outputs the coordinates for the fourier points for a 33 term sum. 
%ignoring the zeroth order terms.
%32 is probably much larger than needed for this purpose. 

xp=linspace(.05,0.1260,13);
yp=linspace(0,0.2,7);%-0.2:0.4/32:0.2;
zp=linspace(0,0.051,7);

coords=zeros(length(xp)*length(yp)*length(zp),3);
iiii=1;
for i= 1:length(xp)
for ii = 1:length(yp)
    for iii = 1:length(zp)
    coords(iiii,:)=[xp(i),yp(ii),zp(iii)];
        iiii=iiii+1;
    end
end
end
%end

