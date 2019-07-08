
% Script to look at QZS hith cuboids and triangular prisms
%
% James O'Connell 5th June 2019

clc;
% close all;
clear;

figure;

Height = 0.04;

cuboid = @(l,b,h) [l/2,b/2,h/2;l/2,-b/2,h/2;-l/2,b/2,h/2;-l/2,-b/2,h/2;l/2,b/2,-h/2;l/2,-b/2,-h/2;-l/2,b/2,-h/2;-l/2,-b/2,-h/2];
triprism = @(l,b,h) [0,b/2,h/2;0,-b/2,h/2;l/2,b/2,-h/2;l/2,-b/2,-h/2;-l/2,b/2,-h/2;-l/2,-b/2,-h/2];
pyramid = @(l,b,h) [0,0,h/2;l/2,b/2,-h/2;l/2,-b/2,-h/2;-l/2,b/2,-h/2;-l/2,-b/2,-h/2];    
frustum = @(bnominal,h,theta) [(bnominal-h/tand(theta))/2,(bnominal-h/tand(theta))/2,h/2;...
    -(bnominal-h/tand(theta))/2,(bnominal-h/tand(theta))/2,h/2;...
    (bnominal-h/tand(theta))/2,-(bnominal-h/tand(theta))/2,h/2;...
    -(bnominal-h/tand(theta))/2,-(bnominal-h/tand(theta))/2,h/2;...
    (bnominal+h/tand(theta))/2,(bnominal+h/tand(theta))/2,-h/2;...
    -(bnominal+h/tand(theta))/2,(bnominal+h/tand(theta))/2,-h/2;...
    (bnominal+h/tand(theta))/2,-(bnominal+h/tand(theta))/2,-h/2;...
    -(bnominal+h/tand(theta))/2,-(bnominal+h/tand(theta))/2,-h/2];

bc = 0.02;
hc = 0.02;
magnetc{1} = cuboid(bc,bc,hc) + repmat([0,0,-hc/2],8,1);
magnetc{2} = cuboid(bc,bc,hc) + repmat([0,0,Height+hc/2],8,1);

b = 0.008;
h = 0.005;
theta = 130;
% magnetp{1} = frustum(b,h,theta) + repmat([0,0,-h/2],8,1);
% magnetp{2} = -frustum(b,h,theta) + repmat([0,0,0.04+h/2],8,1);

h = hc
b = sqrt(3*bc^2*hc/h)
hcc = 0.00;
V = pyramid(b,b,h);
magnetp{1} = [-V] + repmat([0,0,-((6*hcc^2+4*h*hcc+h^2)/(12*hcc+4*h))-hc/4],5,1);
magnetp{2} = [V] + repmat([0,0,Height],5,1) + repmat([0,0,hc/4+((6*hcc^2+4*h*hcc+h^2)/(12*hcc+4*h))],5,1);
% magnetp{2} = pyramid(b,b,h) + repmat([0,0,Height+hc/2+h/4],5,1);

bf = bc;
hfc = hc;
magnetfc = cuboid(bf,bf,hfc) + repmat([0,0,hfc/2],8,1);
dc = linspace(0,Height-hfc,201);
dc = dc(2:end-1);
dc = [zeros(length(dc),2),dc'];

bf = bc;
hfp = hc*bc^2/bf^2;
magnetfp = cuboid(bf,bf,hfp) + repmat([0,0,hfp/2],8,1);
dp = linspace(0,Height-hfp,201);
dp = dp(2:end-1);
dp = [zeros(length(dp),2),dp'];

% magnetp{1}(1,3) = -0.04;

sc1 = alphaShape(magnetc{1},inf);
sp1 = alphaShape(magnetp{1},inf);
sc2 = alphaShape(magnetc{2},inf);
sp2 = alphaShape(magnetp{2},inf);
sfc = alphaShape(magnetfc,inf);
sfp = alphaShape(magnetfp,inf);

volume(sc1)
volume(sp1)

magup = [0,0,1.3];
magdown = [0,0,-1.3];

Fctemp = polyhedronForce(magnetc,magnetfc,{magdown,magup},magup,16,mean(magnetfc),dc);
Fc = Fctemp(:,3);
Fptemp = polyhedronForce(magnetp,magnetfp,{magdown,magup},magup,16,mean(magnetfp),dp);
Fp = Fptemp(:,3);


% for i = 1:length(dc)
%     magnetFc = magnetfc + repmat([0,0,dc(i)],8,1);
%     
%     Fctemp = polyhedronForce(magnetc,magnetFc,{magdown,magup},magup,16,mean(magnetfc));
%     Fc(i) = Fctemp(3);
%     
% end
% 
% for i = 1:length(dp)
%     
%     magnetFp = magnetfp + repmat([0,0,dp(i)],8,1);
%     
%     Fptemp = polyhedronForce(magnetp,magnetFp,{magdown,magup},magup,16,mean(magnetfp));
%     Fp(i) = Fptemp(3);
% end

dFE = 0.002:0.002:0.018;
FcFE = [154.79,120.5,100.98,90.658,87.394,90.668,100.98,120.51,154.76];
FpFE = [130.1,113.29,101.49,94.561,92.28,94.582,101.49,113.3,130.1];
dFE = dFE(1:length(FpFE));

% Fc
% Fp

dd = dc(2,3)-dc(1,3);
Kc = -[(Fc(2)-Fc(1))/dd;(Fc(3:end)-Fc(1:end-2))/(2*dd);(Fc(end)-Fc(end-1))/dd];
Kp = -[(Fp(2)-Fp(1))/dd;(Fp(3:end)-Fp(1:end-2))/(2*dd);(Fp(end)-Fp(end-1))/dd];

dc = dc(:,3);
dp = dp(:,3);

subplot(2,3,[1,2]);
plot(dc,Fc,dp,Fp,dFE,FcFE,'k.',dFE,FpFE,'k.');
grid on;
legend('Cuboid','Pyramid','FEA results');
title('Force-displacement curve');
xlabel('Vertical displacement (m)');
ylabel('Force (N)');
axis tight;
V = axis;
% axis([0,V(1)+V(2),0,1.1*V(4)]);

sfc.Points = sfc.Points + repmat([0,0,mean(dc)],length(sfc.Points),1);
subplot(2,3,3);
plot(sc1);
hold on;
plot(sc2);
plot(sfc);
title('Cuboidal geometry');

subplot(2,3,[4,5]);
plot(dc,Kc,dp,Kp);
hold on
plot((dFE(2:end)+dFE(1:end-1))/2,-diff(FcFE)/0.002,'k.');
plot((dFE(2:end)+dFE(1:end-1))/2,-diff(FpFE)/0.002,'k.');
grid on;
legend('Cuboid','Pyramid','FEA results');
title('Stiffness-displacement curve');
xlabel('Vertical displacement (m)');
ylabel('Stiffness (N/m)');
axis tight;
V = axis;
% axis([0,V(1)+V(2),V(3)/1.1,V(4)*1.1]);

sfp.Points = sfp.Points + repmat([0,0,mean(dc)],length(sfp.Points),1);
subplot(2,3,6);
plot(sp1);
hold on;
plot(sp2);
plot(sfp);
title('Pyramid geometry');
