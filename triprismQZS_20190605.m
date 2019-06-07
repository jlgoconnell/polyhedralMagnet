
% Script to look at QZS hith cuboids and triangular prisms
%
% James O'Connell 5th June 2019

clc;
% close all;
clear;

figure;

Height = 0.2;

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

hc = 0.04;
magnetc{1} = cuboid(hc,hc,hc) + repmat([0,0,-hc/2],8,1);
magnetc{2} = cuboid(hc,hc,hc) + repmat([0,0,Height+hc/2],8,1);

b = 0.005;
h = 0.01;
theta = 130;
% magnetp{1} = frustum(b,h,theta) + repmat([0,0,-h/2],8,1);
% magnetp{2} = -frustum(b,h,theta) + repmat([0,0,0.04+h/2],8,1);

h = hc;%2/3*hc;
b = sqrt(3*hc^3/h);
magnetp{1} = -pyramid(b,b,h) + repmat([0,0,-h/4-hc/2],5,1);
magnetp{2} = pyramid(b,b,h) + repmat([0,0,Height+hc/2+h/4],5,1);

hf = 0.4*hc;
magnetf = cuboid(hc,hc,hf) + repmat([0,0,hf/2],8,1);

sc1 = alphaShape(magnetc{1},inf);
sp1 = alphaShape(magnetp{1},inf);
sc2 = alphaShape(magnetc{2},inf);
sp2 = alphaShape(magnetp{2},inf);
sf = alphaShape(magnetf,inf);

volume(sc1)
volume(sp1)

mag = [0,0,1.3];

d = linspace(0,Height-hf,202);
d = d(2:end-1);

for i = 1:length(d)
    magnetF = magnetf + repmat([0,0,d(i)],8,1);
    
    FFc = polyhedronForce(magnetc{1},magnetF,-mag,mag,12,mean(magnetf));
    FFc = FFc + polyhedronForce(magnetc{2},magnetF,mag,mag,16,mean(magnetf));
    Fc(i) = FFc(3);
    
    FFp = polyhedronForce(magnetp{1},magnetF,-mag,mag,12,mean(magnetf));
    FFp = FFp + polyhedronForce(magnetp{2},magnetF,mag,mag,16,mean(magnetf));
    Fp(i) = FFp(3);
    
    
    
    
    
end

Fc
Fp

dd = d(2)-d(1);
Kc = [(Fc(2)-Fc(1))/dd,(Fc(3:end)-Fc(1:end-2))/(2*dd),(Fc(end)-Fc(end-1))/dd];
Kp = [(Fp(2)-Fp(1))/dd,(Fp(3:end)-Fp(1:end-2))/(2*dd),(Fp(end)-Fp(end-1))/dd];

subplot(2,3,[1,2]);
plot(d,Fc,d,Fp);
grid on;
legend('Cuboid','Pyramid');
title('Force-displacement curve');
xlabel('Vertical displacement (m)');
ylabel('Force (N)');
axis tight;

sf.Points = sf.Points + repmat([0,0,mean(d)],length(sf.Points),1);
subplot(2,3,3);
plot(sc1);
hold on;
plot(sc2);
plot(sf);
title('Cuboidal geometry');

subplot(2,3,[4,5]);
plot(d,Kc,d,Kp);
grid on;
legend('Cuboid','Pyramid');
title('Stiffness-displacement curve');
xlabel('Vertical displacement (m)');
ylabel('Stiffness (N/m)');
axis tight;

subplot(2,3,6);
plot(sp1);
hold on;
plot(sp2);
plot(sf);
title('Pyramid geometry');
