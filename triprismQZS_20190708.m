
% Script to mess with magnets using multi-magnet systems and displacement
% stuff implemented in my code

% close all;
clear;
clc;

cuboid = @(l,b,h) [l/2,b/2,h/2;l/2,-b/2,h/2;-l/2,b/2,h/2;-l/2,-b/2,h/2;l/2,b/2,-h/2;l/2,-b/2,-h/2;-l/2,b/2,-h/2;-l/2,-b/2,-h/2];
pyramid = @(l,b,h) [0,0,h/2;l/2,b/2,-h/2;l/2,-b/2,-h/2;-l/2,b/2,-h/2;-l/2,-b/2,-h/2];

magnetfloat = cuboid(0.02,0.02,0.02);
magnitudefloat = [0,0,1];

n = 201;
d = [zeros(n+2,2),linspace(-0.01,0.01,n+2)'];
d = d(2:end-1,:);

b = sqrt(3)*0.02;
b = 0.02;
V = pyramid(b,b,0.02);
% V = cuboid(0.02,0.02,0.02);
magnet1 = -V + repmat([0,0,-0.02],length(V),1) + repmat([0,0,-0.01],length(V),1);
magnitude1 = [0,0,-1];
magnet2 = V + repmat([0,0,0.02],length(V),1) + repmat([0,0,0.01],length(V),1);
magnitude2 = [0,0,1];
magnet3 = magnetfloat + repmat([0.03,0,0],length(magnetfloat),1);
magnitude3 = [0,0,1];
magnet4 = [-magnet3(:,1),magnet3(:,2:3)];
magnitude4 = magnitude3;


magnetfixed = {magnet1,magnet2,magnet3,magnet4};
magnitudefixed = {magnitude1,magnitude2,magnitude3,magnitude4};

% magnetfixed = {magnet1,magnet2};
% magnitudefixed = {magnitude1,magnitude2};
% 
% magnetfixed = {magnet2};
% magnitudefixed = {magnitude2};

[F,T,t] = polyhedronForce(magnetfixed,magnetfloat,magnitudefixed,magnitudefloat,8,mean(magnetfloat),d);

mf = alphaShape(magnetfloat,inf);
figure;
plot(mf);
hold on;
for i = 1:length(magnetfixed)
    m = alphaShape(magnetfixed{i},inf);
    plot(m);
end

figure;
plot(d(:,3),F(:,3));