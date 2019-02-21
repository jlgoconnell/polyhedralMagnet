
% Script to create two randomly generated tetrahedra with random
% magnetisations and calculate the force between them.
%
% James O'Connell 21st Feb 2019

close all;
clear all;
clc;

% Define parameters
verticesA = 0.1*rand(4,3);
verticesB = 0.1*rand(4,3);
magA = rand(1,3);
magA = magA/norm(magA);
magB = rand(1,3);
magB = magB/norm(magB);

% Make sure they don't intersect by having one above the other
zdiff = max(verticesA(:,3))-min(verticesB(:,3));
verticesB = verticesB + repmat([0,0,1.01*zdiff],length(verticesB),1);

% Actually, let's do Akoun and Yonnet's geometry
cuboid = @(l,b,w) [l/2,b/2,w/2;l/2,-b/2,w/2;-l/2,b/2,w/2;-l/2,-b/2,w/2;l/2,b/2,-w/2;l/2,-b/2,-w/2;-l/2,b/2,-w/2;-l/2,-b/2,-w/2];
verticesA = cuboid(0.02,0.012,0.006);
verticesB = cuboid(0.012,0.02,0.006) + repmat([-0.004,-0.004,0.008],8,1);
magA = [0,0,0.38];
magB = magA;
torquepoint = mean(verticesB);

phi = (1+sqrt(5))/2;
A = 0.01*phi;
verticesA = A*[1,1,1;1,1,-1;1,-1,1;1,-1,-1;-1,1,1;-1,1,-1;-1,-1,1;-1,-1,-1;...
    0,phi,1/phi;0,phi,-1/phi;0,-phi,1/phi;0,-phi,-1/phi;...
    1/phi,0,phi;1/phi,0,-phi;-1/phi,0,phi;-1/phi,0,-phi;...
    phi,1/phi,0;phi,-1/phi,0;-phi,1/phi,0;-phi,-1/phi,0];
thetax = 31.715*pi/180;
R_x = [1,0,0;0,cos(thetax),-sin(thetax);0,sin(thetax),cos(thetax)];
verticesA = verticesA*R_x;
verticesB = verticesA + repmat([0,0,0.05],length(verticesA),1);
torquepoint = mean(verticesB);
magA = [0,0,1];
magB = [1,0,0];

Sa = alphaShape(verticesA,Inf);
Sb = alphaShape(verticesB,Inf);
figure;
plot(Sa);
hold on;
plot(Sb);

[F,T] = polyhedronForce(verticesA,verticesB,magA,magB,mean(verticesB),1e-50,300)
[Fold,Told] = polyhedronForceOld(verticesA,verticesB,magA,magB,torquepoint)
% 
% figure(2);
% plot(1:size(F,1)-1,F(2:end,1),'r-');
% grid on;
% figure(3);
% plot(1:size(F,1)-1,F(2:end,2),'r-');
% grid on;
% figure(4);
% plot(1:size(F,1)-1,F(2:end,3),'r-');
% grid on;

% magnet_fixed.dim = [0.02 0.012 0.006];
% magnet_float.dim = [0.012 0.02 0.006];
% magnet_fixed.magn = 0.38;
% magnet_float.magn = 0.38;
% magnet_fixed.magdir = [0 0 1]; % z
% magnet_float.magdir = [0 0 1]; % z
% offset = [-0.004;-0.004;0.008];
% displ = 0;
% displ_range = offset+[1; 0; 0]*displ;
% 
% f1_xyz = magnetforces(magnet_fixed,magnet_float,displ_range);
% t1_xyz = magnetforces(magnet_fixed,magnet_float,displ_range,'torque');
