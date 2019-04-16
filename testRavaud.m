
% Script to test my results against those of Ravaud et al 2010

close all;
clc;
clear;

circleApproxN = 10;
theta = linspace(0,2*pi,circleApproxN+1)';
theta = theta(1:end-1);


%%

r1 = 0.0875;
r2 = 0.0875;
h1 = 0.025;
h2 = 0.025;
k1 = 8000;
k2 = 8000;
magA = [0,0,k1*4*pi*10^-7];
magB = [0,0,-k2*4*pi*10^-7];
verticesA = [r1*cos(theta),r1*sin(theta),zeros(size(theta));r1*cos(theta),r1*sin(theta),-h1*ones(size(theta))];

n = 100;
d = linspace(0,1,n+1);
d = d(2:end);

for i = 1:length(d)
    i
    verticesB = [r2*cos(theta),r2*sin(theta),zeros(size(theta));r2*cos(theta),r2*sin(theta),h2*ones(size(theta))]+repmat([0,0,d(i)],length(theta)*2,1);
    [F(i,:),~,t(i,:)] = polyhedronForceGauss(verticesA,verticesB,magA,magB,10,mean(verticesB));
end



%%

r1 = 0.0875;
r2 = 0.085;
h1 = 0.025;
h2 = 0.02;
magA = [0,0,1];
magB = [0,0,-1];
verticesA = [r1*cos(theta),r1*sin(theta),zeros(size(theta));r1*cos(theta),r1*sin(theta),-h1*ones(size(theta))];

n = 100;
d = linspace(0,0.5,n+1);
d = d(2:end);
for i = 1:length(d)
    i
    verticesB = [r2*cos(theta),r2*sin(theta),zeros(size(theta));r2*cos(theta),r2*sin(theta),h2*ones(size(theta))]+repmat([0,0,d(i)],length(theta)*2,1);
    [F(i,:),~,t(i,:)] = polyhedronForceGauss(verticesA,verticesB,magA,magB,10,mean(verticesB));
end