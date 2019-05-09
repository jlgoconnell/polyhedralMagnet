
% Script to compare the force between approximated cylinders using my code,
% to exact cylinders using Will's code
%
% James O'Connell 9th May 2019

clear;
close all;
clc;

magA = [0,0,-1];
magB = [0,0,1];
n = 17;
theta = linspace(0,2*pi,n);
theta = theta(1:end-1)';
h = 0.02;
r = 0.02;
verticesA = [r*cos(theta),r*sin(theta),zeros(size(theta));r*cos(theta),r*sin(theta),-h*ones(size(theta))];

maga = struct('type','cylinder','magn',1,'dim',[r,h],'magdir',magA,'isring',0,'dir',[0,0,1]);
magb = struct('type','cylinder','magn',1,'dim',[r,h],'magdir',magB,'isring',0,'dir',[0,0,1]);

d = 0.001:0.001:0.01;

for i = 1:length(d)
    verticesB = [verticesA(:,1:2),-verticesA(:,3)+d(i)*ones(length(verticesA),1)];
    
    [Ftemp,~,tapprox(i)] = polyhedronForce(verticesA,verticesB,magA,magB,6,mean(verticesB));
    Fapprox(i,:) = Ftemp;
    
    Ftemp = magnetforces(maga,magb,[0,0,h+d(i)]);
    Fexact(i,:) = Ftemp;
    
    
end