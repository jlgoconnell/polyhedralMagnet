
% Script to try to get a linear magnet force
%
% James O'Connell, 3rd June 2019

% This uses a trapezial magnet at the bottom and a cuboid at the top

% -------
% |     |
% |  ^  |  B
% |  |  |
% -------
%
% |\   
% | \  
% |  \
% |   \    A
% | -> \
% -------
%    h

clear;
% close all;
clc;

% Set up magnets

cuboid = @(l,b,w) [l/2,b/2,w/2;l/2,-b/2,w/2;-l/2,b/2,w/2;-l/2,-b/2,w/2;l/2,b/2,-w/2;l/2,-b/2,-w/2;-l/2,b/2,-w/2;-l/2,-b/2,-w/2];

theta = 45;
h = 0.015;
verticesA = cuboid(h,0.01,0.01) + repmat([0,0,0.005],8,1);
verticesA(1:2,:) = [];
verticesA(1:2,3) = h*tand(theta);
verticesA = verticesA + repmat([0,0,-max(verticesA(:,3))],length(verticesA),1);
sa = alphaShape(verticesA,inf);
myverticesB = cuboid(0.01,0.01,0.01) + repmat([0,0,0.005],8,1);
sb = alphaShape(myverticesB,inf);
magA = [cosd(theta),0,-sind(theta)];
magB = [0,0,-1];

d = linspace(0,0.05,51);
d = d(2:end);

F = zeros(length(d),1);
for i = 1:length(d)
    verticesB = myverticesB + repmat([0,0,d(i)],length(myverticesB),1);
    
    FF = polyhedronForce(verticesA,verticesB,magA,magB,16,mean(verticesB));
    F(i,:) = FF(3);
    
end

% figure;
% plot(sa);
% hold on;
% plot(sb);

% Estimate stiffness
dd = d(2)-d(1);
K = zeros(size(F));
K(2:end-1) = (F(3:end)-F(1:end-2))/(2*dd);
K(1) = (F(2)-F(1))/dd;
K(end) = (F(end)-F(end-1))/dd;

figure;
subplot(2,1,1);
plot(d,F);
title(['h = ',num2str(h),' \theta = ',num2str(theta)]);
grid on;
subplot(2,1,2);
plot(d,K);
grid on;
