
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


verticesC = cuboid(0.02,0.02,0.01) + repmat([0,0,0.005],8,1);
verticesC = verticesC + repmat([0,0,-max(verticesC(:,3))],length(verticesC),1);
verticesD = cuboid(0.02,0.02,0.01) + repmat([0.033,0,0.03-0.005],length(verticesC),1);
verticesE = [-verticesD(:,1),verticesD(:,2:3)];

magnet{1} = verticesC;
magnet{2} = verticesD;
magnet{3} = verticesE;
dd = 0.06;
magnet{4} = magnet{1} + repmat([0,0,dd],length(magnet{1}),1);
magnet{5} = verticesE*[0,-1,0;1,0,0;0,0,1];
magnet{6} = verticesD*[0,-1,0;1,0,0;0,0,1];

myverticesB = cuboid(0.01,0.01,0.01) + repmat([0,0,0.005],8,1);
magB = [0,0,1];

magnitude{1} = [0,0,1];
magnitude{2} = [0,0,1];
magnitude{3} = [0,0,1];
magnitude{4} = [0,0,1];
magnitude{5} = [0,0,1];
magnitude{6} = [0,0,1];

d = linspace(0,dd-0.02,51);
d = d(2:end-1);

F = zeros(length(d),1);
for i = 1:length(d)
    verticesB = myverticesB + repmat([0,0,d(i)],length(myverticesB),1);
    
    FF = zeros(1,3);
    for j = 1:length(magnet)
        FF = FF + polyhedronForce(magnet{j},verticesB,magnitude{j},magB,8,mean(verticesB));
    end
    F(i,:) = FF(3);
    
end


% Estimate stiffness
dd = d(2)-d(1);
K = zeros(size(F));
K(2:end-1) = (F(3:end)-F(1:end-2))/(2*dd);
K(1) = (F(2)-F(1))/dd;
K(end) = (F(end)-F(end-1))/dd;

figure;
subplot(2,2,1);
plot(d,F);
ylabel('Force');
% title(['h = ',num2str(h),' \theta = ',num2str(theta)]);
grid on;
subplot(2,2,3);
plot(d,K);
ylabel('Stiffness');
grid on;

subplot(2,2,[2,4]);
s = alphaShape(verticesB,inf);
plot(s);
hold on;
for i = 1:length(magnet)
    s = alphaShape(magnet{i},inf);
    plot(s);
end