function [out,ang,mesherror] = fullscaleComsol3DSYMFerroCloak2FerroOpt(ang,Bmsr,ferroNum,onum)
%ferroNum is now an array of the z position of the bottom of the metglas strips 
% FullScale3DSYMFerroCloak.m
%
% Model exported on Feb 7 2018, 16:50 by COMSOL 5.3.0.260.
if nargin < 4
 onum=false;
end

if onum
    zpos(1)=-2.54/100;
end

import com.comsol.model.*
import com.comsol.model.util.*
mesherror=false;
model = ModelUtil.create('Model');

ModelUtil.showProgress(true);

model.modelPath('/data1/cmswank/comsol/matlab/COMSOLScripts');

model.label('FullScale3DSYMFerroCloak.mph');

model.comments(['Untitled\n\n']);

model.component.create('comp1', false);

model.component('comp1').geom.create('geom1', 3);

model.result.table.create('evl3', 'Table');

model.component('comp1').mesh.create('mesh1');
blockheight=2.8;
blockwidth=2;
model.component('comp1').geom('geom1').selection.create('csel1', 'CumulativeSelection');
model.component('comp1').geom('geom1').selection('csel1').label('Outer Boundary');
model.component('comp1').geom('geom1').selection.create('csel2', 'CumulativeSelection');
model.component('comp1').geom('geom1').selection('csel2').label('Lead Shield');
model.component('comp1').geom('geom1').selection.create('csel3', 'CumulativeSelection');
model.component('comp1').geom('geom1').selection('csel3').label('Lead Endcap');
model.component('comp1').geom('geom1').selection.create('csel4', 'CumulativeSelection');
model.component('comp1').geom('geom1').selection('csel4').label('Wire');
model.component('comp1').geom('geom1').selection.create('csel40', 'CumulativeSelection');
model.component('comp1').geom('geom1').selection('csel40').label('DoubleCurrentWire');
model.component('comp1').geom('geom1').selection.create('csel5', 'CumulativeSelection');
model.component('comp1').geom('geom1').selection('csel5').label('Cell');
model.component('comp1').geom('geom1').create('blk1', 'Block');
model.component('comp1').geom('geom1').feature('blk1').set('contributeto', 'csel1');
model.component('comp1').geom('geom1').feature('blk1').set('pos', [0 0 0]);
model.component('comp1').geom('geom1').feature('blk1').set('size', [blockwidth blockwidth blockheight]);
model.component('comp1').geom('geom1').create('blk2', 'Block');
model.component('comp1').geom('geom1').feature('blk2').set('contributeto', 'csel1');
model.component('comp1').geom('geom1').feature('blk2').set('pos', [0 0 0]);
model.component('comp1').geom('geom1').feature('blk2').set('size', [blockwidth blockwidth blockheight]);
model.component('comp1').geom('geom1').create('cyl1', 'Cylinder');
model.component('comp1').geom('geom1').feature('cyl1').set('type', 'surface');
model.component('comp1').geom('geom1').feature('cyl1').set('r', '(40.889)*2.54/100');
model.component('comp1').geom('geom1').feature('cyl1').set('h', '94*2.54/100/2');
%model.component('comp1').geom('geom1').create('ic8', 'InterpolationCurve');
%model.component('comp1').geom('geom1').feature('ic8').set('contributeto', 'csel4');
%model.component('comp1').geom('geom1').feature('ic8').set('table', {'36.5*2.54/100*cos(pi/3)' '36.5*2.54/100*sin(pi/3)' '94*2.54/100/2+2*2.54/100'; '(36.5-0.9-0.6)*2.54/100*cos(pi/3)' '(36.5-0.9-0.6)*2.54/100*sin(pi/3)' '94*2.54/100/2+2*2.54/100'});
model.component('comp1').geom('geom1').create('int1', 'Intersection');
model.component('comp1').geom('geom1').feature('int1').set('contributeto', 'csel2');
model.component('comp1').geom('geom1').feature('int1').selection('input').set({'blk2' 'cyl1'});
model.component('comp1').geom('geom1').create('blk3', 'Block');
model.component('comp1').geom('geom1').feature('blk3').set('contributeto', 'csel1');
model.component('comp1').geom('geom1').feature('blk3').set('pos', [0 0 0]);
model.component('comp1').geom('geom1').feature('blk3').set('size', [blockwidth blockwidth blockheight]);
model.component('comp1').geom('geom1').create('wp1', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp1').set('quickz', '94*2.54/100/2');
model.component('comp1').geom('geom1').feature('wp1').set('unite', true);
model.component('comp1').geom('geom1').feature('wp1').geom.create('c1', 'Circle');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('c1').set('r', '(39.731)*2.54/100');
model.component('comp1').geom('geom1').feature('wp1').geom.create('c2', 'Circle');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('c2').set('r', '(31.259-1.4)*2.54/100');
model.component('comp1').geom('geom1').feature('wp1').geom.create('c3', 'Circle');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('c3').set('r', '31.259*2.54/100');
% model.component('comp1').geom('geom1').feature('wp1').geom.create('c4', 'Circle');
% model.component('comp1').geom('geom1').feature('wp1').geom.feature('c4').set('r', '(35.9-1.8)*2.54/100');
% model.component('comp1').geom('geom1').feature('wp1').geom.create('c5', 'Circle');
% model.component('comp1').geom('geom1').feature('wp1').geom.feature('c5').set('r', '(35.9-1.8-1.06)*2.54/100');
% model.component('comp1').geom('geom1').feature('wp1').geom.create('c6', 'Circle');
% model.component('comp1').geom('geom1').feature('wp1').geom.feature('c6').set('r', '(35.9-1.8-1.06-1.8)*2.54/100');
% model.component('comp1').geom('geom1').feature('wp1').geom.create('c7', 'Circle');
% model.component('comp1').geom('geom1').feature('wp1').geom.feature('c7').set('r', '(35.9-1.8-1.06-1.8-1.4)*2.54/100');
model.component('comp1').geom('geom1').feature('wp1').geom.create('co1', 'Compose');

%%%This can be used to make lead annuli in the lead 
%model.component('comp1').geom('geom1').feature('wp1').geom.feature('co1').active(false);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('co1').set('formula', 'c1-c3+c2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


model.component('comp1').geom('geom1').create('int2', 'Intersection');
model.component('comp1').geom('geom1').feature('int2').set('contributeto', 'csel3');
model.component('comp1').geom('geom1').feature('int2').selection('input').set({'blk3' 'wp1'});
model.component('comp1').geom('geom1').create('blk5', 'Block');
model.component('comp1').geom('geom1').feature('blk5').set('contributeto', 'csel1');
model.component('comp1').geom('geom1').feature('blk5').set('pos', [0 0 0]);
model.component('comp1').geom('geom1').feature('blk5').set('size', [blockwidth blockwidth blockheight]);
model.component('comp1').geom('geom1').create('blk4', 'Block');
model.component('comp1').geom('geom1').feature('blk4').set('contributeto', 'csel5');
model.component('comp1').geom('geom1').feature('blk4').set('pos', {'.05' '0' '0'});
model.component('comp1').geom('geom1').feature('blk4').set('size', [0.076 0.2 0.051]);


%CREATE PHYSICS MODEL HERE!!!! THIS will allow us to easily create saddle
%currents!
model.component('comp1').physics.create('mf', 'InductionCurrents', 'geom1');

ang=sort(ang,'descend'); %counting current is easier for me to understand this way. 

 


%%%fix beyond bounds stuff. 

if any(ang<=0)
    fxzero=sum(ang<=0);
    for i = length(ang)-(fxzero+1):length(ang)
        ang(i)=1E-3*(length(ang)-i+1);
    end
end

if any(ang>=pi/2)
    fxpio2=sum(ang>=pi/2);
    for i = 1:fxpio2
        ang(i)=pi/2-1E-3*(i-1);
    end
end

%disp(ang);


%too close proximity correction (two methods)

%Push Away Method. 
%this code pushes the close proximity wires away from each other. 
% if any(diff(ang)<1E-2)
%   dbind=find(diff(ang)>-1E-2)+1;
%   for i = dbind
%       if i>round(length(ang)/2)
%       ang(i-1)=ang(i-1)+.075;
%       else
%           ang(i)=ang(i)-0.075;
%       end
%   end
% end

%Double Current Method (maybe better?!?)
%this code averages the angle and to be doubled later. doubles the current. stores the new as
%angout. 

% if any(diff(ang)<1E-2)
%     dbind=find(diff(ang)>-1E-2)+1;
%     
%     for i = dbind
%     ang(i-1)=ang(i)/2+ang(i-1)/2;
%     ang(i)=ang(i-1);
%     end
%    
% end
%%%% after much effort it was determined both of these methods sucked.

disp([Bmsr,ferroNum]);
ang=sort(ang,'descend');




%
%%%
%%%%%      Wire Geometry starts here
doublebool=false;
firstsaddle=true;
pbholes=['int2'];
for i = 1:length(ang)
    name1=['wire',num2str(i),'_1'];
    name2=['wire',num2str(i),'_2'];
    name3=['wire',num2str(i),'_3'];
    namesadI=['cselsad',num2str(i)];
    namedc=['edcsad',num2str(i)];
    namepbhole=['pbhole',num2str(i)];
    pbholes=[pbholes,'-',namepbhole];
   name1=['wire',num2str(i),'_1'];
model.component('comp1').geom('geom1').selection.create(namesadI, 'CumulativeSelection');
model.component('comp1').geom('geom1').selection(namesadI).label(['Saddle Current ', num2str(i)]);

%create lead cut outs. 
model.component('comp1').geom('geom1').create(namepbhole, 'Block');
model.component('comp1').geom('geom1').feature(namepbhole).set('pos', {['37.773*2.54/100*cos(',num2str(ang(i)),'), 37.773*2.54/100*sin(',num2str(ang(i)),'), 94*2.54/100/2']});
model.component('comp1').geom('geom1').feature(namepbhole).set('rot', ['(180/pi)*',num2str(ang(i))]);
model.component('comp1').geom('geom1').feature(namepbhole).set('base', 'center');
model.component('comp1').geom('geom1').feature(namepbhole).set('size', {'0.75*2.54/100' '0.25*2.54/100' '0.05'});





if i<length(ang)
%%%saddle return (parametric curve section)
    if ang(i+1)==ang(i)
   
    else
        if firstsaddle==true
        model.component('comp1').geom('geom1').create(name1, 'ParametricCurve');
        model.component('comp1').geom('geom1').feature(name1).set('parmin', num2str(ang(i+1)));
        model.component('comp1').geom('geom1').feature(name1).set('parmax', num2str(ang(i)));
        model.component('comp1').geom('geom1').feature(name1).set('coord', {'(37.773+2)*2.54/100*cos(s)' '(37.773+2)*2.54/100*sin(s)' '94*2.54/100/2+2*2.54/100'});
        model.component('comp1').geom('geom1').feature(name1).set('contributeto', namesadI); 
        %%%
    %%% saddle currents physics
        model.component('comp1').physics('mf').create(namedc, 'EdgeCurrent', 1);
        model.component('comp1').physics('mf').feature(namedc).selection.named(['geom1_',namesadI,'_edg']);
        model.component('comp1').physics('mf').feature(namedc).set('Ie', num2str(i));
        firstsaddle=false;
        else
        model.component('comp1').geom('geom1').create(name1, 'ParametricCurve');
        model.component('comp1').geom('geom1').feature(name1).set('parmin', num2str(ang(i+1)));
        model.component('comp1').geom('geom1').feature(name1).set('parmax', num2str(ang(i)));
        model.component('comp1').geom('geom1').feature(name1).set('coord', {'(37.773+2)*2.54/100*cos(s)' '(37.773+2)*2.54/100*sin(s)' '94*2.54/100/2+2*2.54/100'});
        model.component('comp1').geom('geom1').feature(name1).set('contributeto', namesadI); 
        %%%
    %%% saddle currents physics
        model.component('comp1').physics('mf').create(namedc, 'EdgeCurrent', 1);
        model.component('comp1').physics('mf').feature(namedc).selection.named(['geom1_',namesadI,'_edg']);
        model.component('comp1').physics('mf').feature(namedc).set('Ie', num2str(-i));
            
        end
    end
else
model.component('comp1').geom('geom1').create(name1, 'ParametricCurve');
model.component('comp1').geom('geom1').feature(name1).set('parmin','0');
model.component('comp1').geom('geom1').feature(name1).set('parmax', num2str(ang(i)));
model.component('comp1').geom('geom1').feature(name1).set('coord', {'(37.773+2)*2.54/100*cos(s)' '(37.773+2)*2.54/100*sin(s)' '94*2.54/100/2+2*2.54/100'});
model.component('comp1').geom('geom1').feature(name1).set('contributeto', namesadI); 
%%%
%%% saddle currents physics
model.component('comp1').physics('mf').create(namedc, 'EdgeCurrent', 1);
model.component('comp1').physics('mf').feature(namedc).selection.named(['geom1_',namesadI,'_edg']);
model.component('comp1').physics('mf').feature(namedc).set('Ie', num2str(-i));


end

    
  if i<length(ang)
%Bring wires to saddle winds
    if ang(i+1)==ang(i)  
    doublebool=true;
    else
    
        if doublebool==true
        model.component('comp1').geom('geom1').create(name2, 'InterpolationCurve');
        %%%%%%%%%%coutribute to double current  
        model.component('comp1').geom('geom1').feature(name2).set('contributeto', 'csel40');
        model.component('comp1').geom('geom1').feature(name2).set('table', {['37.773*2.54/100*cos(',num2str(ang(i)),')'],['37.773*2.54/100*sin(',num2str(ang(i)),')'], '0';...
                                                                    ['37.773*2.54/100*cos(',num2str(ang(i)),')'],['37.773*2.54/100*sin(',num2str(ang(i)),')'], '94*2.54/100/2+2*2.54/100'});
        %contribute to double current!
        model.component('comp1').geom('geom1').create(name3, 'InterpolationCurve');
        model.component('comp1').geom('geom1').feature(name3).set('contributeto', 'csel40');
        model.component('comp1').geom('geom1').feature(name3).set('table', {['37.773*2.54/100*cos(',num2str(ang(i)),')'],['37.773*2.54/100*sin(',num2str(ang(i)),')'],'94*2.54/100/2+2*2.54/100'; ...
                                                                            ['(37.773+2)*2.54/100*cos(',num2str(ang(i)),')'],['(37.773+2)*2.54/100*sin(',num2str(ang(i)),')'], '94*2.54/100/2+2*2.54/100'});
        doublebool=false;
        
        else
    
        model.component('comp1').geom('geom1').create(name2, 'InterpolationCurve');
        model.component('comp1').geom('geom1').feature(name2).set('contributeto', 'csel4');
        model.component('comp1').geom('geom1').feature(name2).set('table', {['37.773*2.54/100*cos(',num2str(ang(i)),')'],['37.773*2.54/100*sin(',num2str(ang(i)),')'], '0';...
                                                                            ['37.773*2.54/100*cos(',num2str(ang(i)),')'],['37.773*2.54/100*sin(',num2str(ang(i)),')'], '94*2.54/100/2+2*2.54/100'});
        model.component('comp1').geom('geom1').create(name3, 'InterpolationCurve');
        model.component('comp1').geom('geom1').feature(name3).set('contributeto', 'csel4');
        model.component('comp1').geom('geom1').feature(name3).set('table', {['37.773*2.54/100*cos(',num2str(ang(i)),')'],['37.773*2.54/100*sin(',num2str(ang(i)),')'],'94*2.54/100/2+2*2.54/100'; ...
                                                                            ['(37.773+2)*2.54/100*cos(',num2str(ang(i)),')'],['(37.773+2)*2.54/100*sin(',num2str(ang(i)),')'], '94*2.54/100/2+2*2.54/100'});
    
        
         end
                                                                
    end
    
  else
      if doublebool==true
           model.component('comp1').geom('geom1').create(name2, 'InterpolationCurve');
        %%%%%%%%%%coutribute to double current  
        model.component('comp1').geom('geom1').feature(name2).set('contributeto', 'csel40');
        model.component('comp1').geom('geom1').feature(name2).set('table', {['37.773*2.54/100*cos(',num2str(ang(i)),')'],['37.773*2.54/100*sin(',num2str(ang(i)),')'], '0';...
                                                                    ['37.773*2.54/100*cos(',num2str(ang(i)),')'],['37.773*2.54/100*sin(',num2str(ang(i)),')'], '94*2.54/100/2+2*2.54/100'});
        %contribute to double current!
        model.component('comp1').geom('geom1').create(name3, 'InterpolationCurve');
        model.component('comp1').geom('geom1').feature(name3).set('contributeto', 'csel40');
        model.component('comp1').geom('geom1').feature(name3).set('table', {['37.773*2.54/100*cos(',num2str(ang(i)),')'],['37.773*2.54/100*sin(',num2str(ang(i)),')'],'94*2.54/100/2+2*2.54/100'; ...
                                                                            ['(37.773+2)*2.54/100*cos(',num2str(ang(i)),')'],['(37.773+2)*2.54/100*sin(',num2str(ang(i)),')'], '94*2.54/100/2+2*2.54/100'});
        doublebool=false;
      
      else
          model.component('comp1').geom('geom1').create(name2, 'InterpolationCurve');
        model.component('comp1').geom('geom1').feature(name2).set('contributeto', 'csel4');
        model.component('comp1').geom('geom1').feature(name2).set('table', {['37.773*2.54/100*cos(',num2str(ang(i)),')'],['37.773*2.54/100*sin(',num2str(ang(i)),')'], '0';...
                                                                            ['37.773*2.54/100*cos(',num2str(ang(i)),')'],['37.773*2.54/100*sin(',num2str(ang(i)),')'], '94*2.54/100/2+2*2.54/100'});
        model.component('comp1').geom('geom1').create(name3, 'InterpolationCurve');
        model.component('comp1').geom('geom1').feature(name3).set('contributeto', 'csel4');
        model.component('comp1').geom('geom1').feature(name3).set('table', {['37.773*2.54/100*cos(',num2str(ang(i)),')'],['37.773*2.54/100*sin(',num2str(ang(i)),')'],'94*2.54/100/2+2*2.54/100'; ...
                                                                            ['(37.773+2)*2.54/100*cos(',num2str(ang(i)),')'],['(37.773+2)*2.54/100*sin(',num2str(ang(i)),')'], '94*2.54/100/2+2*2.54/100'});
      
      end
  end
  
  
end
%apply current to the wires
model.component('comp1').physics('mf').create('edcDouble1', 'EdgeCurrent', 1);   
model.physics('mf').feature('edcDouble1').selection.named('geom1_csel40_edg');
model.physics('mf').feature('edcDouble1').set('Ie', '2');

model.component('comp1').physics('mf').create('edc1', 'EdgeCurrent', 1);   
model.physics('mf').feature('edc1').selection.named('geom1_csel4_edg');
model.physics('mf').feature('edc1').set('Ie', '1');
%%%%% end of wires.
%%%
%

%%%%holes in the lead simulation
%%%add heat flush tube 
model.component('comp1').geom('geom1').create('hfhole', 'Cylinder');
model.component('comp1').geom('geom1').feature('hfhole').set('r', '(4)*2.54/100');
model.component('comp1').geom('geom1').feature('hfhole').set('h', '.05');
model.component('comp1').geom('geom1').feature('hfhole').set('pos', {'8.4853*2.54/100' '8.4853*2.54/100' '94*2.54/100/2-0.025'});



%this is a formula: endcap-sum(holes)  
    pbholes=[pbholes,'-','hfhole'];

%%%%put holes in endcap (still contributes fullscaleFourierB03D(angle,1,1);to endcap). 
model.component('comp1').geom('geom1').create('co2', 'Compose');
model.component('comp1').geom('geom1').feature('co2').set('contributeto', 'csel3');
model.component('comp1').geom('geom1').feature('co2').set('formula', pbholes);
%model.component('comp1').geom('geom1').run;




%Top lead wire cover
model.component('comp1').geom('geom1').create('cyl2', 'Cylinder');
model.component('comp1').geom('geom1').feature('cyl2').set('pos', {'0' '0' '94*2.54/100/2'});
model.component('comp1').geom('geom1').feature('cyl2').set('type', 'surface');
model.component('comp1').geom('geom1').feature('cyl2').set('r', '(40.889)*2.54/100');
model.component('comp1').geom('geom1').feature('cyl2').set('h', '4*2.54/100');
model.component('comp1').geom('geom1').create('cyl3', 'Cylinder');
model.component('comp1').geom('geom1').feature('cyl3').set('pos', {'0' '0' '94*2.54/100/2'});
model.component('comp1').geom('geom1').feature('cyl3').set('type', 'surface');
model.component('comp1').geom('geom1').feature('cyl3').set('r', '(31.259)*2.54/100');
model.component('comp1').geom('geom1').feature('cyl3').set('h', '4*2.54/100');
model.component('comp1').geom('geom1').create('wp2', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp2').set('quickz', '94*2.54/100/2+4*2.54/100');
model.component('comp1').geom('geom1').feature('wp2').set('unite', true);
model.component('comp1').geom('geom1').feature('wp2').geom.create('c1', 'Circle');
model.component('comp1').geom('geom1').feature('wp2').geom.feature('c1').set('r', '(40.889)*2.54/100');
model.component('comp1').geom('geom1').feature('wp2').geom.create('c2', 'Circle');
model.component('comp1').geom('geom1').feature('wp2').geom.feature('c2').set('r', '(31.259)*2.54/100');
model.component('comp1').geom('geom1').feature('wp2').geom.create('co1', 'Compose');
model.component('comp1').geom('geom1').feature('wp2').geom.feature('co1').set('formula', 'c1-c2');
model.component('comp1').geom('geom1').create('co1', 'Compose');
model.component('comp1').geom('geom1').feature('co1').set('formula', 'cyl2+cyl3+wp2');
model.component('comp1').geom('geom1').create('int3', 'Intersection');
model.component('comp1').geom('geom1').feature('int3').selection('input').set({'blk5' 'co1'});
model.component('comp1').geom('geom1').feature('int3').set('contributeto', 'csel2');



%metglas shield gee i'm a tree
model.component('comp1').geom('geom1').selection.create('cselmet','CumulativeSelection');
model.component('comp1').geom('geom1').selection('cselmet').label('metglas flux return');
model.component('comp1').geom('geom1').create('cyl4', 'Cylinder');
model.component('comp1').geom('geom1').feature('cyl4').set('type', 'surface');
model.component('comp1').geom('geom1').feature('cyl4').set('r', '(39.861)*2.54/100');
model.component('comp1').geom('geom1').feature('cyl4').set('h', '96*2.54/100/2');
model.component('comp1').geom('geom1').create('blk6', 'Block');
model.component('comp1').geom('geom1').feature('blk6').set('contributeto', 'csel1');
model.component('comp1').geom('geom1').feature('blk6').set('pos', [0 0 0]);
model.component('comp1').geom('geom1').feature('blk6').set('size', [blockwidth blockwidth blockheight]);
model.component('comp1').geom('geom1').create('int4', 'Intersection');
model.component('comp1').geom('geom1').feature('int4').selection('input').set({'blk6' 'cyl4'});
model.component('comp1').geom('geom1').feature('int4').set('contributeto', 'cselmet');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Metglas shield fix 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Metglas MSR  Shielding fix radius. 
%'(31.259-1.4)*2.54/100'
% model.component('comp1').geom('geom1').create('cylmetfix', 'Cylinder');
% model.component('comp1').geom('geom1').feature('cylmetfix').set('type', 'surface');
% model.component('comp1').geom('geom1').feature('cylmetfix').set('r', '(31.259-1.4)*2.54/100');
% model.component('comp1').geom('geom1').feature('cylmetfix').set('h', '10*2.54/100');
% model.component('comp1').geom('geom1').feature('cylmetfix').set('pos', [0 0 94*2.54/100/2]);
% model.component('comp1').geom('geom1').create('blkmetfix', 'Block');
% model.component('comp1').geom('geom1').feature('blkmetfix').set('contributeto', 'csel1');
% model.component('comp1').geom('geom1').feature('blkmetfix').set('pos', [0 0 0]);
% model.component('comp1').geom('geom1').feature('blkmetfix').set('size', [blockwidth blockwidth blockheight]);
% model.component('comp1').geom('geom1').create('intmetfix', 'Intersection');
% model.component('comp1').geom('geom1').feature('intmetfix').selection('input').set({'blkmetfix' 'cylmetfix'});
% model.component('comp1').geom('geom1').feature('intmetfix').set('contributeto', 'cselmet');



%Neutron beam hole. 
% model.component('comp1').geom('geom1').create('blk72hole', 'Block');
% model.component('comp1').geom('geom1').feature('blk72hole').set('contributeto', 'csel1');
% model.component('comp1').geom('geom1').feature('blk72hole').set('pos', [0.045 (39.861)*2.54/100-3*2.54/100 0]);
% model.component('comp1').geom('geom1').feature('blk72hole').set('size', [8.6/100 3*5.08/100 5.6/100]);
% 
% model.component('comp1').geom('geom1').create('metshield4', 'Difference');
% model.component('comp1').geom('geom1').feature('metshield4').selection('input').set({'int4'});
% model.component('comp1').geom('geom1').feature('metshield4').selection('input2').set({'blk72hole'});
% model.component('comp1').geom('geom1').run('metshield4');
% model.component('comp1').geom('geom1').feature('metshield4').set('contributeto', 'cselmet');
% 
% %metglas shield gee i'm a tree
% model.component('comp1').geom('geom1').selection.create('cselmethole','CumulativeSelection');
% model.component('comp1').geom('geom1').selection('cselmethole').label('metglas flux return thin section');
% model.component('comp1').geom('geom1').create('cyl4hole', 'Cylinder');
% model.component('comp1').geom('geom1').feature('cyl4hole').set('type', 'surface');
% model.component('comp1').geom('geom1').feature('cyl4hole').set('r', '(39.861)*2.54/100');
% model.component('comp1').geom('geom1').feature('cyl4hole').set('h', '96*2.54/100/2');
% model.component('comp1').geom('geom1').create('blk6hole', 'Block');
% model.component('comp1').geom('geom1').feature('blk6hole').set('contributeto', 'csel1');
% model.component('comp1').geom('geom1').feature('blk6hole').set('pos', [0 0 0]);
% model.component('comp1').geom('geom1').feature('blk6hole').set('size', [blockwidth blockwidth blockheight]);
% model.component('comp1').geom('geom1').create('int4hole', 'Intersection');
% model.component('comp1').geom('geom1').feature('int4hole').selection('input').set({'blk6hole' 'cyl4hole'});
% 
% 
% 
% 
% model.component('comp1').geom('geom1').create('blk7hole', 'Block');
% model.component('comp1').geom('geom1').feature('blk7hole').set('contributeto', 'csel1');
% model.component('comp1').geom('geom1').feature('blk7hole').set('pos', [0 (39.861)*2.54/100-3*2.54/100 0]);
% model.component('comp1').geom('geom1').feature('blk7hole').set('size', [13.1/100 3*5.08/100 5.6/100]);
% 
% 
% model.component('comp1').geom('geom1').create('int5hole', 'Intersection');
% model.component('comp1').geom('geom1').feature('int5hole').selection('input').set({'blk7hole' 'int4hole'});
% %model.component('comp1').geom('geom1').feature('int5hole').set('contributeto', 'cselmethole');
% 
% model.component('comp1').geom('geom1').create('blk8hole', 'Block');
% model.component('comp1').geom('geom1').feature('blk8hole').set('contributeto', 'csel1');
% model.component('comp1').geom('geom1').feature('blk8hole').set('pos', [0 (39.861)*2.54/100-3*2.54/100 0]);
% model.component('comp1').geom('geom1').feature('blk8hole').set('size', [13.1/100 3*5.08/100 0.05/100]);
% 
% model.component('comp1').geom('geom1').create('metshield42', 'Difference');
% model.component('comp1').geom('geom1').feature('metshield42').selection('input').set({'int5hole'});
% model.component('comp1').geom('geom1').feature('metshield42').selection('input2').set({'blk8hole'});
% model.component('comp1').geom('geom1').run('metshield42');
% model.component('comp1').geom('geom1').feature('metshield42').set('contributeto', 'cselmethole');

%Radius of Metglas shield fix
%'(31.259-1.4)*2.54/100'


%metglas Shield Physics
model.component('comp1').physics('mf').create('ms1', 'MagneticShielding', 2);
model.component('comp1').physics('mf').feature('ms1').selection.named('geom1_cselmet_bnd');
model.component('comp1').physics('mf').feature('ms1').set('ds', '.116[mm]');
model.component('comp1').physics('mf').feature('ms1').set('murbnd', 60000);
model.component('comp1').physics('mf').feature('ms1').set('coordinateSystem', 'sys1');
model.component('comp1').physics('mf').feature('ms1').set('murbnd_mat', 'userdef');

%metglas Shield Physics
% model.component('comp1').physics('mf').create('mshole', 'MagneticShielding', 2);
% model.component('comp1').physics('mf').feature('mshole').selection.named('geom1_cselmethole_bnd');
% model.component('comp1').physics('mf').feature('mshole').set('ds', '.015[mm]');
% model.component('comp1').physics('mf').feature('mshole').set('murbnd', 1);
% model.component('comp1').physics('mf').feature('mshole').set('coordinateSystem', 'sys1');
% model.component('comp1').physics('mf').feature('mshole').set('murbnd_mat', 'userdef');





%Copper shield gee i'm a tree use when sd is active. %%%MUST CREATE PHYS
% model.component('comp1').geom('geom1').selection.create('cselcop','CumulativeSelection');
% model.component('comp1').geom('geom1').selection('cselcop').label('copper rf flux return');
% model.component('comp1').geom('geom1').create('cyl48', 'Cylinder');
% model.component('comp1').geom('geom1').feature('cyl48').set('type', 'surface');
% model.component('comp1').geom('geom1').feature('cyl48').set('r', '(34.126)*2.54/100');
% model.component('comp1').geom('geom1').feature('cyl48').set('h', '94*2.54/100/2');
% model.component('comp1').geom('geom1').create('blk68', 'Block');
% model.component('comp1').geom('geom1').feature('blk68').set('contributeto', 'csel1');
% model.component('comp1').geom('geom1').feature('blk68').set('pos', [0 0 0]);
% model.component('comp1').geom('geom1').feature('blk68').set('size', [blockwidth blockwidth blockheight]);
% model.component('comp1').geom('geom1').create('int48', 'Intersection');
% model.component('comp1').geom('geom1').feature('int48').selection('input').set({'blk68' 'cyl48'});
% model.component('comp1').geom('geom1').feature('int48').set('contributeto', 'cselcop');


%metglas Shield Physics (USE TO MAKE COPPER PHYSICS?) (NO PROBABLY USE MIBC
% model.component('comp1').physics('mf').create('ms1', 'MagneticShielding', 2);
% model.component('comp1').physics('mf').feature('ms1').selection.named('geom1_cselmet_bnd');
% model.component('comp1').physics('mf').feature('ms1').set('ds', '.012[mm]');
% model.component('comp1').physics('mf').feature('ms1').set('murbnd', 60000);
% model.component('comp1').physics('mf').feature('ms1').set('coordinateSystem', 'sys1');
% model.component('comp1').physics('mf').feature('ms1').set('murbnd_mat', 'userdef');






%%%ferrocloak geometry! (radius correct for upset design. )
L=98*2.54/100;
%spacez=L-2.54/100*2*ferroNum;
%zspace=spacez/(2*ferroNum-1);
model.component('comp1').geom('geom1').selection.create('cselferro','CumulativeSelection');
model.component('comp1').geom('geom1').selection('cselferro').label('ferro cloak');

for i = 1:length(ferroNum)
zpos=ferroNum(i);

nameblock=['block',num2str(i)];
cylname=['ferromet',num2str(i)];
intname=['ferroint',num2str(i)];

model.component('comp1').geom('geom1').create(cylname, 'Cylinder');
model.component('comp1').geom('geom1').feature(cylname).set('type', 'surface');
model.component('comp1').geom('geom1').feature(cylname).set('r', '(40.889)*2.54/100+0.005');
model.component('comp1').geom('geom1').feature(cylname).set('h', '2*2.54/100');
model.component('comp1').geom('geom1').feature(cylname).set('pos', {'0','0',num2str(zpos)});
model.component('comp1').geom('geom1').create(nameblock, 'Block');
model.component('comp1').geom('geom1').feature(nameblock).set('contributeto', 'csel1');
model.component('comp1').geom('geom1').feature(nameblock).set('pos', [0 0 0]);
model.component('comp1').geom('geom1').feature(nameblock).set('size', [blockwidth blockwidth blockheight]);
model.component('comp1').geom('geom1').create(intname, 'Intersection');
model.component('comp1').geom('geom1').feature(intname).selection('input').set({nameblock cylname});
model.component('comp1').geom('geom1').feature(intname).set('contributeto', 'cselferro');



end
%%%ferro cloak physics!
model.component('comp1').physics('mf').create('ms2', 'MagneticShielding', 2);
model.component('comp1').physics('mf').feature('ms2').selection.named('geom1_cselferro_bnd');
model.component('comp1').physics('mf').feature('ms2').set('ds', '.029[mm]');
model.component('comp1').physics('mf').feature('ms2').set('murbnd', 60000);
model.component('comp1').physics('mf').feature('ms2').set('coordinateSystem', 'sys1');
model.component('comp1').physics('mf').feature('ms2').set('murbnd_mat', 'userdef');

try

%model.component('comp1').geom('geom1').run;

model.component('comp1').material.create('mat1', 'Common');
model.component('comp1').material('mat1').propertyGroup('def').func.create('eta', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func.create('Cp', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func.create('rho', 'Analytic');
model.component('comp1').material('mat1').propertyGroup('def').func.create('k', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func.create('cs', 'Analytic');
model.component('comp1').material('mat1').propertyGroup.create('RefractiveIndex', 'Refractive index');

%Outer Geometry (so that our boundary conditions are defined!, this may
%have been the trouble with the previous version)

model.component('comp1').geom('geom1').selection.create('csel6', 'CumulativeSelection');
model.component('comp1').geom('geom1').selection('csel6').label('CoppperShield');
model.component('comp1').geom('geom1').selection.create('csel7', 'CumulativeSelection');
model.component('comp1').geom('geom1').selection('csel7').label('PerfectMagneticConductor');
model.component('comp1').geom('geom1').selection.create('csel8', 'CumulativeSelection');
model.component('comp1').geom('geom1').selection('csel8').label('InsulationBC');
model.component('comp1').geom('geom1').selection.create('csel9', 'CumulativeSelection');
model.component('comp1').geom('geom1').selection('csel9').label('BfieldMSR');

model.component('comp1').geom('geom1').create('wp3', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp3').set('contributeto', 'csel7');
model.component('comp1').geom('geom1').feature('wp3').set('quickplane', 'yz');
model.component('comp1').geom('geom1').feature('wp3').set('unite', true);
model.component('comp1').geom('geom1').feature('wp3').geom.create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp3').geom.feature('r1').set('size', [blockwidth blockheight]);
model.component('comp1').geom('geom1').create('wp4', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp4').set('contributeto', 'csel8');
model.component('comp1').geom('geom1').feature('wp4').set('quickplane', 'xz');
model.component('comp1').geom('geom1').feature('wp4').set('unite', true);
model.component('comp1').geom('geom1').feature('wp4').geom.create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp4').geom.feature('r1').set('size', [blockwidth blockheight]);
model.component('comp1').geom('geom1').create('wp5', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp5').set('contributeto', 'csel8');
model.component('comp1').geom('geom1').feature('wp5').set('unite', true);
model.component('comp1').geom('geom1').feature('wp5').geom.create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp5').geom.feature('r1').set('size', [blockwidth blockwidth]);
model.component('comp1').geom('geom1').create('wp6', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp6').set('contributeto', 'csel9');
model.component('comp1').geom('geom1').feature('wp6').set('quickz', blockheight);
model.component('comp1').geom('geom1').feature('wp6').set('unite', true);
model.component('comp1').geom('geom1').feature('wp6').geom.create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp6').geom.feature('r1').set('size', [blockwidth blockwidth]);
model.component('comp1').geom('geom1').create('wp7', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp7').set('contributeto', 'csel9');
model.component('comp1').geom('geom1').feature('wp7').set('quickplane', 'xz');
model.component('comp1').geom('geom1').feature('wp7').set('quicky', blockwidth);
model.component('comp1').geom('geom1').feature('wp7').set('unite', true);
model.component('comp1').geom('geom1').feature('wp7').geom.create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp7').geom.feature('r1').set('size', [blockwidth blockheight]);
model.component('comp1').geom('geom1').create('wp8', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp8').set('contributeto', 'csel9');
model.component('comp1').geom('geom1').feature('wp8').set('quickplane', 'yz');
model.component('comp1').geom('geom1').feature('wp8').set('quickx', blockwidth);
model.component('comp1').geom('geom1').feature('wp8').set('unite', true);
model.component('comp1').geom('geom1').feature('wp8').geom.create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp8').geom.feature('r1').set('size', [blockwidth blockheight]);
model.component('comp1').geom('geom1').run;






%%%%%%%%%%%%%%%%%%%%       PHYSICS               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%model.component('comp1').physics.create('mf', 'InductionCurrents', 'geom1');
%This was created eariler, above the wire creation: above is a copy. 

% This is done correctly now, no source of error. 

model.component('comp1').physics('mf').create('pmc1', 'PerfectMagneticConductor', 2);
model.component('comp1').physics('mf').feature('pmc1').selection.named('geom1_csel7_bnd');
model.component('comp1').physics('mf').create('mi2', 'MagneticInsulation', 2);
%model.component('comp1').physics('mf').feature('mi2').active(false);
model.component('comp1').physics('mf').feature('mi2').selection.named('geom1_csel3_bnd');
model.component('comp1').physics('mf').create('mi3', 'MagneticInsulation', 2);
model.component('comp1').physics('mf').feature('mi3').selection.named('geom1_csel2_bnd');
model.component('comp1').physics('mf').create('mi4', 'MagneticInsulation', 2);
model.component('comp1').physics('mf').feature('mi4').selection.named('geom1_csel8_bnd');
model.component('comp1').physics('mf').create('mfb1', 'MagneticFieldBoundary', 2);
model.component('comp1').physics('mf').feature('mfb1').selection.named('geom1_csel9_bnd');


% %%% I don't even know what this is!
% model.component('comp1').physics('mf').feature('al1').set('MsJA', {'1e6'; '0'; '0'; '0'; '1e6'; '0'; '0'; '0'; '1e6'});
% model.component('comp1').physics('mf').feature('al1').set('aJA', [200; 0; 0; 0; 200; 0; 0; 0; 200]);
% model.component('comp1').physics('mf').feature('al1').set('kJA', [300; 0; 0; 0; 300; 0; 0; 0; 300]);


model.component('comp1').physics('mf').feature('mfb1').set('H0', [Bmsr; 0; 0]);








% %%%%%%%%%%%%%%%%% BS
% model.result.table('evl3').label('Evaluation 3D');
% model.result.table('evl3').comments('Interactive 3D values');
% 
% model.capeopen.label('Thermodynamics Package');
% 
% model.component('comp1').view('view1').set('transparency', true);
% model.component('comp1').view('view2').axis.set('xmin', -1.5437811613082886);
% model.component('comp1').viewfullscaleFourierB03D(angle,3E-6,22);('view2').axis.set('xmax', 3.3022208213806152);
% model.component('comp1').view('view2').axis.set('ymin', -1.1096162796020508);
% model.component('comp1').view('view2').axis.set('ymax', 1.2305188179016113);
% model.component('comp1').view('view2').axis.set('abstractviewlratio', -0.22216369211673737);
% model.component('comp1').view('view2').axis.set('abstractviewrratio', 0.9067254662513733);
% model.component('comp1').view('view2').axis.set('abstractviewbratio', -0.18845641613006592);
% model.component('comp1').view('view2').axis.set('abstractviewtratio', 0.0710793137550354);
% model.component('comp1').view('view2').axis.set('abstractviewxscale', 0.005085375625640154);
% model.component('comp1').view('view2').axis.set('abstractviewyscale', 0.005085376091301441);
% model.component('comp1').view('view3').axis.set('xmin', -1.9861934185028076);
% model.component('comp1').view('view3').axis.set('xmax', 2.1574933528900146);
% model.component('comp1').view('view3').axis.set('ymin', -1.1401351690292358);
% model.component('comp1').view('view3').axis.set('ymax', 1.3114351034164429);
% model.component('comp1').view('view3').axis.set('abstractviewlratio', -0.42961961030960083);
% model.component('comp1').view('view3').axis.set('abstractviewrratio', 0.42961961030960083);
% model.component('comp1').view('view3').axis.set('abstractviewbratio', -0.05000002309679985);
% model.component('comp1').view('view3').axis.set('abstractviewtratio', 0.05000002309679985);
% model.component('comp1').view('view3').axis.set('abstractviewxscale', 0.004441250581294298);
% model.component('comp1').view('view3').axis.set('abstractviewyscale', 0.004441250581294298);
% %%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%     Materials!
model.component('comp1').material('mat1').label('Air');
model.component('comp1').material('mat1').set('family', 'air');
model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('pieces', {'200.0' '1600.0' '-8.38278E-7+8.35717342E-8*T^1-7.69429583E-11*T^2+4.6437266E-14*T^3-1.06585607E-17*T^4'});
model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('pieces', {'200.0' '1600.0' '1047.63657-0.372589265*T^1+9.45304214E-4*T^2-6.02409443E-7*T^3+1.2858961E-10*T^4'});
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('expr', 'pA*0.02897/8.314/T');
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('args', {'pA' 'T'});
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('dermethod', 'manual');
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('argders', {'pA' 'd(pA*0.02897/8.314/T,pA)'; 'T' 'd(pA*0.02897/8.314/T,T)'});
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('plotargs', {'pA' '0' '1'; 'T' '0' '1'});
model.component('comp1').material('mat1').propertyGroup('def').func('k').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('k').set('pieces', {'200.0' '1600.0' '-0.00227583562+1.15480022E-4*T^1-7.90252856E-8*T^2+4.11702505E-11*T^3-7.43864331E-15*T^4'});
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('expr', 'sqrt(1.4*287*T)');
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('args', {'T'});
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('dermethod', 'manual');
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('argders', {'T' 'd(sqrt(1.4*287*T),T)'});
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('plotargs', {'T' '0' '1'});
model.component('comp1').material('mat1').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.component('comp1').material('mat1').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.component('comp1').material('mat1').propertyGroup('def').set('dynamicviscosity', 'eta(T[1/K])[Pa*s]');
model.component('comp1').material('mat1').propertyGroup('def').set('ratioofspecificheat', '1.4');
model.component('comp1').material('mat1').propertyGroup('def').set('electricconductivity', {'0[S/m]' '0' '0' '0' '0[S/m]' '0' '0' '0' '0[S/m]'});
model.component('comp1').material('mat1').propertyGroup('def').set('heatcapacity', 'Cp(T[1/K])[J/(kg*K)]');
model.component('comp1').material('mat1').propertyGroup('def').set('density', 'rho(pA[1/Pa],T[1/K])[kg/m^3]');
model.component('comp1').material('mat1').propertyGroup('def').set('thermalconductivity', {'k(T[1/K])[W/(m*K)]' '0' '0' '0' 'k(T[1/K])[W/(m*K)]' '0' '0' '0' 'k(T[1/K])[W/(m*K)]'});
model.component('comp1').material('mat1').propertyGroup('def').set('soundspeed', 'cs(T[1/K])[m/s]');
model.component('comp1').material('mat1').propertyGroup('def').addInput('temperature');
model.component('comp1').material('mat1').propertyGroup('def').addInput('pressure');
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', '');
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('ki', '');
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('ki', {'0' '0' '0' '0' '0' '0' '0' '0' '0'});




%%%%%%%%%%%%%%%%  SOLVER
model.study.create('std1');
model.study('std1').create('stat', 'Stationary');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').create('i1', 'Iterative');
model.sol('sol1').feature('s1').feature('i1').create('mg1', 'Multigrid');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').create('so1', 'SOR');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').create('so1', 'SOR');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').create('ams1', 'AMS');
model.sol('sol1').feature('s1').feature.remove('fcDef');

%%% bs continued
% model.result.dataset.create('surf1', 'Surface');
% model.result.dataset.create('surf2', 'Surface');
% model.result.dataset('surf1').selection.set([1 2 4 5 9 17 27 28 29 30 31 32 33 35]);
% model.result.dataset('surf2').selection.set([6 8 10 11 12 14 16 20 23]);
% model.result.create('pg1', 'PlotGroup3D');
% model.result('pg1').create('surf1', 'Surface');
% model.result('pg1').create('surf2', 'Surface');


%%%%solver  %trying to stop the damn mesh errors. 
model.sol('sol1').attach('std1');
model.sol('sol1').feature('s1').feature('i1').set('linsolver', 'fgmres');
try
    
    model.component('comp1').mesh('mesh1').autoMeshSize(4);
    model.sol('sol1').runAll;
catch
    
    try
    model.component('comp1').mesh('mesh1').autoMeshSize(3);
    model.sol('sol1').runAll;
    catch
            try
            model.component('comp1').mesh('mesh1').autoMeshSize(2);
            model.sol('sol1').runAll;
            catch
            model.component('comp1').mesh('mesh1').autoMeshSize(5);
            model.sol('sol1').runAll;
                try
                model.component('comp1').mesh('mesh1').autoMeshSize(1);
                model.sol('sol1').runAll;
                catch
                mesherror=true;    
                end
                    
            end
    end 
end

catch
    mesherror=true;
end
%%% more bs
% model.result('pg1').label('Magnetic Flux Density Norm (mf)');
% model.result('pg1').feature('surf1').set('data', 'surf1');
% model.result('pg1').feature('surf1').set('expr', 'mf.Bx');
% model.result('pg1').feature('surf1').set('descr', 'Magnetic flux density, x component');
% model.result('pg1').feature('surf1').set('resolution', 'normal');
% model.result('pg1').feature('surf2').set('data', 'surf2');
% model.result('pg1').feature('surf2').set('expr', 'mf.normJs');
% model.result('pg1').feature('surf2').set('unit', 'A/m');
% model.result('pg1').feature('surf2').set('descr', 'Surface current density norm');
% model.result('pg1').feature('surf2').set('resolution', 'normal');

out = model;
%[Bx,By,Bz] = mphinterp(model,{'mf.Bx','mf.By','mf.Bz'},'coord',coords');