
% Script to compare the results obtained from my code for Akoun and
% Yonnet's geometry
%
% James O'Connell 13th May 2019

clear;
close all;
clc;

cuboid = @(l,b,w) [l/2,b/2,w/2;l/2,-b/2,w/2;-l/2,b/2,w/2;-l/2,-b/2,w/2;l/2,b/2,-w/2;l/2,-b/2,-w/2;-l/2,b/2,-w/2;-l/2,-b/2,-w/2];

verticesA = cuboid(0.02,0.012,0.006);
verticesBinit = cuboid(0.012,0.02,0.006) + repmat([-0.004,-0.004,0.008],8,1);
magA = [0,0,0.38];
magB = magA;
torquepoint = mean(verticesBinit);
tolerance = 1e-3;
timeout = 1;

maga = magnetdefine('type','cuboid','dim',[0.02 0.012 0.006],'magn',0.38,'magdir','z');
magb = magnetdefine('type','cuboid','dim',[0.012 0.02 0.006],'magn',0.38,'magdir','z');

d = 0:0.002:0.03;

for i = 1:length(d)
    verticesB = verticesBinit + d(i)*[1,0,0];
    
    [Ftemp,~,tapprox(i)] = polyhedronForce(verticesA,verticesB,magA,magB,12,mean(verticesB));
    Fapprox(i,:) = Ftemp;
    
    
end

tapprox'

exactd = 0:0.0001:0.03;
Fexact = magnetforces(maga,magb,repmat([-0.004; -0.004; 0.008],1,length(exactd))'+exactd'*[1,0,0])';
Fexactd = magnetforces(maga,magb,repmat([-0.004; -0.004; 0.008],1,length(d))'+d'*[1,0,0])';

d = d*1000;
exactd = exactd*1000;

close figure 1;

figure;
p1 = plot(d,Fapprox,'ko');
hold on;
% plot(d,Ffea(:,3),'k.');
p2 = plot(exactd,Fexact,'k-');
grid on;
legend([p1(1),p2(1)],{'This work','Analytic solution'});
xlabel('Separation distance, d (mm)');
ylabel('z-Force (N)');
title('z-Force between two axially magnetised cylindrical permanent magnets');

errorapprox = abs(Fapprox-Fexactd);
% errorfea = abs(Ffea(:,3)-Fexactd(:,3));
mean(errorapprox)
% mean(errorfea)

% figure
% plot(d,errorapprox)
% hold on;
% % plot(d,errorfea)
% grid on;
% 
% errorapproxpc = errorapprox./Fexactd(:,3);
% % errorfeapc = errorfea./Fexactd(:,3);
% figure;
% plot(d,errorapproxpc);
% hold on;
% % plot(d,errorfeapc);
% grid on;