
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
frustum = @(bnominal,h,theta) [(bnominal-h/tand(theta))/2,(bnominal-h/tand(theta))/2,h/2;...
    -(bnominal-h/tand(theta))/2,(bnominal-h/tand(theta))/2,h/2;...
    (bnominal-h/tand(theta))/2,-(bnominal-h/tand(theta))/2,h/2;...
    -(bnominal-h/tand(theta))/2,-(bnominal-h/tand(theta))/2,h/2;...
    (bnominal+h/tand(theta))/2,(bnominal+h/tand(theta))/2,-h/2;...
    -(bnominal+h/tand(theta))/2,(bnominal+h/tand(theta))/2,-h/2;...
    (bnominal+h/tand(theta))/2,-(bnominal+h/tand(theta))/2,-h/2;...
    -(bnominal+h/tand(theta))/2,-(bnominal+h/tand(theta))/2,-h/2];
pyramid = @(l,b,h) [0,0,h/2;l/2,b/2,-h/2;l/2,-b/2,-h/2;-l/2,b/2,-h/2;-l/2,-b/2,-h/2];    

% verticesC = frustum(0.0194,0.01,130) + repmat([0,0,0.005],8,1); % Frustum
verticesC = cuboid(0.02,0.02,0.01) + repmat([0,0,0.005],8,1); % Cuboid
verticesC = verticesC + repmat([0,0,-max(verticesC(:,3))],length(verticesC),1);
h = 0.016;
verticesD = cuboid(0.02,0.02,h) + repmat([0.039,0,0.025],length(verticesC),1); % Cuboid
% verticesD = frustum(0.02,0.01,70)*[0,0,1;0,1,0;-1,0,0]+repmat([0.033,0,0.025],8,1); % Frustum
verticesE = [-verticesD(:,1),verticesD(:,2:3)];

h = 0.008;
b = 0.035;
verticesC = -pyramid(b,b,h) + repmat([0,0,-h/4-0.005],5,1);

magnet{1} = verticesC;
magnet{2} = verticesD;
magnet{3} = verticesE;
dd = 0.06;
magnet{4} = -magnet{1} + repmat([0,0,dd-0.01],length(magnet{1}),1);
magnet{5} = verticesE*[0,-1,0;1,0,0;0,0,1];
magnet{6} = verticesD*[0,-1,0;1,0,0;0,0,1];

myverticesB = cuboid(0.01,0.01,0.01) + repmat([0,0,0.005],8,1);
magB = [0,0,-1];

magnitude{1} = [0,0,1];
magnitude{2} = [0,0,1];
magnitude{3} = [0,0,1];
magnitude{4} = [0,0,1];
magnitude{5} = [0,0,1];
magnitude{6} = [0,0,1];

d = linspace(0.00,0.04,52);
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
dK = [(K(2)-K(1))/dd;(K(3:end)-K(1:end-2))/(2*dd);(K(end)-K(end-1))/dd];

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
% ylim([600,800]);
subplot(2,2,2);
s = alphaShape(verticesB,inf);
plot(s);
hold on;
for i = 1:length(magnet)
    s = alphaShape(magnet{i},inf);
    plot(s);
end
subplot(2,2,4);
plot(d,dK);
ylabel('d/dx (stiffness)');
grid on;

sum(abs(dK)<10000)/length(d)*100