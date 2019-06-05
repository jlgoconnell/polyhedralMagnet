
% Script to look at QZS hith cuboids and triangular prisms
%
% James O'Connell 5th June 2019

clc;
close all;

cuboid = @(l,b,h) [l/2,b/2,h/2;l/2,-b/2,h/2;-l/2,b/2,h/2;-l/2,-b/2,h/2;l/2,b/2,-h/2;l/2,-b/2,-h/2;-l/2,b/2,-h/2;-l/2,-b/2,-h/2];
triprism = @(l,b,h) [0,b/2,h/2;0,-b/2,h/2;l/2,b/2,-h/2;l/2,-b/2,-h/2;-l/2,b/2,-h/2;-l/2,-b/2,-h/2];
pyramid = @(l,b,h) [0,0,h/2;l/2,b/2,-h/2;l/2,-b/2,-h/2;-l/2,b/2,-h/2;-l/2,-b/2,-h/2];

magnetc{1} = cuboid(0.01,0.01,0.01) + repmat([0,0,-0.005],8,1);
magnetc{2} = cuboid(0.01,0.01,0.01) + repmat([0,0,0.045],8,1);

b = 0.02;
h = 0.03;
magnetp{1} = pyramid(b,b,h) + repmat([0,0,-h/2],5,1);
magnetp{2} = -pyramid(b,b,h) + repmat([0,0,0.04+h/2],5,1);

magnetf = cuboid(0.01,0.01,0.01) + repmat([0,0,0.005],8,1);

mag = [0,0,1];

d = linspace(0,0.03,11);
d = d(2:end-1);

for i = 1:length(d)
    magnetF = magnetf + repmat([0,0,d(i)],8,1);
    
    FFc = polyhedronForce(magnetc{1},magnetF,mag,mag,12,mean(magnetf));
    FFc = FFc + polyhedronForce(magnetc{2},magnetF,mag,mag,12,mean(magnetf));
    Fc(i) = FFc(3);
    
    FFp = polyhedronForce(magnetp{1},magnetF,mag,mag,12,mean(magnetf));
    FFp = FFp + polyhedronForce(magnetp{2},magnetF,mag,mag,12,mean(magnetf));
    Fp(i) = FFp(3);
    
    
    
    
    
end

Fc
Fp

dd = d(2)-d(1);
Kc = [(Fc(2)-Fc(1))/dd,(Fc(3:end)-Fc(1:end-2))/(2*dd),(Fc(end)-Fc(end-1))/dd];
Kp = [(Fp(2)-Fp(1))/dd,(Fp(3:end)-Fp(1:end-2))/(2*dd),(Fp(end)-Fp(end-1))/dd];

subplot(1,2,1);
plot(d,Fc,d,Fp);
grid on;
legend('Cuboid','Prism');

subplot(1,2,2);
plot(d,Kc,d,Kp);
grid on;
legend('Cuboid','Prism');