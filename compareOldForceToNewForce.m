
% Script to compare old and new methods for calculating forces
% 
% James O'Connell 1st April 2019

clear;
close all;
clc;

phi = (1+sqrt(5))/2;
A = 0.01*phi;
verticesA = A*[1,1,1;1,1,-1;1,-1,1;1,-1,-1;-1,1,1;-1,1,-1;-1,-1,1;-1,-1,-1;...
    0,phi,1/phi;0,phi,-1/phi;0,-phi,1/phi;0,-phi,-1/phi;...
    1/phi,0,phi;1/phi,0,-phi;-1/phi,0,phi;-1/phi,0,-phi;...
    phi,1/phi,0;phi,-1/phi,0;-phi,1/phi,0;-phi,-1/phi,0];
thetax = 31.715*pi/180;
R_x = [1,0,0;0,cos(thetax),-sin(thetax);0,sin(thetax),cos(thetax)];
verticesA = verticesA*R_x;
Sa = alphaShape(verticesA,Inf);
magA = [0,0,1];
verticesB = verticesA + repmat([0.002,0.003,0],length(verticesA),1);
verticesB = verticesB + repmat([0,0,0.001-min(verticesB(:,3))+max(verticesA(:,3))],length(verticesB),1);
Sb = alphaShape(verticesB,Inf);
magB = [0,0,-1];
torquept = mean(verticesB);

tic;
Ffft = polyhedronForceFFT(verticesA,verticesB,magA,magB,6)
timeout = toc;

[Fold,~,~,~] = polyhedronForce(verticesA,verticesB,magA,magB,torquept,1e-5,timeout)