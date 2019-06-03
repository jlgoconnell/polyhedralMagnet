
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

d = 0:0.001:0.03;

for i = 1:length(d)
    verticesB = verticesBinit + d(i)*[1,0,0];
    
    [Ftemp,Ttemp,tapprox(i)] = polyhedronForce(verticesA,verticesB,magA,magB,16,mean(verticesB));
    Fapprox(i,:) = Ftemp;
    Tapprox(i,:) = Ttemp;
    
    
end

tapprox'

exactd = 0:0.0001:0.03;
Fexact = magnetforces(maga,magb,repmat([-0.004; -0.004; 0.008],1,length(exactd))'+exactd'*[1,0,0])';
Fexactd = magnetforces(maga,magb,repmat([-0.004; -0.004; 0.008],1,length(d))'+d'*[1,0,0])';
Texact = magnetforces(maga,magb,repmat([-0.004; -0.004; 0.008],1,length(exactd))'+exactd'*[1,0,0],'torque')';
Texactd = magnetforces(maga,magb,repmat([-0.004; -0.004; 0.008],1,length(d))'+d'*[1,0,0],'torque')';

d = d*1000;
exactd = exactd*1000;

close figure 1;

% Import FEA results
FEA = importdata('cuboidSimulationResults.mat');
% I ran simulations which calculated the torque about the wrong point, so
% I'm fixing it here:
FEA.Tfea(:,2) = FEA.Tfea(:,2)+FEA.Ffea(:,3).*d'/1000;
FEA.Tfea(:,3) = FEA.Tfea(:,3)-FEA.Ffea(:,2).*d'/1000;

i = 1:3;

figure;
p1 = plot(d,Fapprox(:,i),'ko');
hold on;
p2 = plot(d,FEA.Ffea(:,i),'k.');
p3 = plot(exactd,Fexact(:,i),'k-');
grid on;
legend([p1(1),p2(1),p3(1)],{'This work','FEA','Analytic solution'});
xlabel('Separation distance, d (mm)');
ylabel('Force (N)');
title('Force between two vertically magnetised cuboidal permanent magnets');

figure;
p1 = plot(d,1000*Tapprox(:,i),'ko');
hold on;
p2 = plot(d,1000*FEA.Tfea(:,i),'k.');
p3 = plot(exactd,1000*Texact(:,i),'k-');
grid on;
legend([p1(1),p2(1),p3(1)],{'This work','FEA','Analytic solution'});
xlabel('Separation distance, d (mm)');
ylabel('Torque (mNm)');
title('Torque between two vertically magnetised cuboidal permanent magnets');

errorapprox = abs(Fapprox-Fexactd);
errorfea = abs(FEA.Ffea-Fexactd);
mean(errorapprox)
mean(errorfea)

figure
semilogy(d,errorapprox,'ko-')
hold on;
semilogy(d,errorfea,'k.-')
grid on;

% errorapproxpc = errorapprox./Fexactd;
% errorfeapc = errorfea./Fexactd;
% figure;
% plot(d,errorapproxpc*100,'ko-');
% hold on;
% plot(d,errorfeapc*100,'ro-');
% grid on;
% ylim([-10,10]);
% legend('Approxx','Approxy','Approxz','FEAx','FEAy','FEAz');

figure;
thing = ['x','y','z'];
for i = 1:3
    subplot(2,3,i);
    p1 = plot(d,Fapprox(:,i),'ko');
    hold on;
    p2 = plot(d,FEA.Ffea(:,i),'k.');
    p3 = plot(exactd,Fexact(:,i),'k-');
    grid on;
%     legend([p1(1),p2(1),p3(1)],{'This work','FEA','Analytic solution'});
    xlabel('Separation distance, d (mm)');
    ylabel([thing(i),'-Force (N)']);
    
    subplot(2,3,i+3);
    p1 = plot(d,1000*Tapprox(:,i),'ko');
    hold on;
    p2 = plot(d,1000*FEA.Tfea(:,i),'k.');
    p3 = plot(exactd,1000*Texact(:,i),'k-');
    grid on;
%     legend([p1(1),p2(1),p3(1)],{'This work','FEA','Analytic solution'});
    xlabel('Separation distance, d (mm)');
    ylabel([thing(i),'-Torque (mNm)']);
end




