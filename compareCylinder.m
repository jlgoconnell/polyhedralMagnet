
% Script to compare the force between approximated cylinders using my code,
% to exact cylinders using Will's code and FEA
%
% James O'Connell 9th May 2019

clear;
close all;
clc;

magA = [0,0,-1];
magB = [0,0,1];
n = 33;
theta = linspace(0,2*pi,n);
theta = theta + diff(theta(1:2));
theta = theta(1:end-1)';
h = 0.02;
r = 0.02;
verticesA = [r*cos(theta),r*sin(theta),zeros(size(theta));r*cos(theta),r*sin(theta),-h*ones(size(theta))];

maga = struct('type','cylinder','magn',1,'dim',[r,h],'magdir',magA,'isring',0,'dir',[0,0,1]);
magb = struct('type','cylinder','magn',1,'dim',[r,h],'magdir',magB,'isring',0,'dir',[0,0,1]);

d = 0.001:0.001:0.02;

for i = 1:length(d)
    verticesB = [verticesA(:,1:2),-verticesA(:,3)+d(i)*ones(length(verticesA),1)];
    
    [Ftemp,~,tapprox(i)] = polyhedronForce(verticesA,verticesB,magA,magB,8,mean(verticesB));
    Fapprox(i,:) = Ftemp;
    
    
end

tapprox'

exactd = 0.001:0.0001:0.02;
Fexact = magnetforces(maga,magb,(h+exactd')*[0,0,1])';
Fexactd = magnetforces(maga,magb,(h+d')*[0,0,1])';

% FEA results:
% These are using two cylinders in the appropriate configuration with a
% percent error of 0.1%.
Ffea = [-0.01183,0.0098378,257.12;... % 98
    -0.038596,0.0064214,223.94;... % 100
    -0.019739,-0.0033181,198.31;... % 102
    -0.022984,0.0035233,177.19;... % 90
    0.041865,-0.0029169,159.43;... % 105
    0.0019984,-0.014451,144.12;... % 105
    0.006158,0.0059783,130.78;... % 100
    0.016622,-0.004528,119.11;... % 109
    -0.008102,0.017648,108.67;... % 94
    0.028015,-0.0073587,99.446;... % 123
    -0.0066658,-0.032518,91.089;... % 95
    -0.0095614,0.016148,83.711;... % 117
    0.0066655,0.000057487,76.929;... % 116
    -0.0011718,-0.0061586,70.855;... % 95
    -0.0016103,0.0094966,65.283;... % 89
    0.0062394,-0.0076713,60.244;... % 97
    -0.0039662,-0.016916,55.712;... % 105
    -0.0029971,-0.01511,51.517;... % 112
    0.0010841,0.018672,47.677;... % 91
    -0.0014862,-0.010364,44.164]; % 91

tfea = [98;100;102;90;105;105;100;109;94;123;95;117;116;95;89;97;105;112;91;91];

d = d*1000;
exactd = exactd*1000;

close figure 1;

figure;
plot(d,Fapprox(:,3),'ko');
hold on;
plot(d,Ffea(:,3),'k.');
plot(exactd,Fexact(:,3),'k-');
grid on;
legend('This work','FEA','Analytic solution');
xlabel('Separation distance, d (mm)');
ylabel('z-Force (N)');
title('z-Force between two cylindrical permanent magnets');

errorapprox = abs(Fapprox(:,3)-Fexactd(:,3));
errorfea = abs(Ffea(:,3)-Fexactd(:,3));
mean(errorapprox)
mean(errorfea)

figure
plot(d,errorapprox)
hold on;
plot(d,errorfea)
grid on;

errorapproxpc = errorapprox./Fexactd(:,3);
errorfeapc = errorfea./Fexactd(:,3);
figure;
plot(d,errorapproxpc);
hold on;
plot(d,errorfeapc);
grid on;





