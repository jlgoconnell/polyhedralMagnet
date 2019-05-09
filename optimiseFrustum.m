
% Script to assess forces between two frustum magnets with varying theta
% and h
%
% James O'Connell 6th May 2019

clear;
close all;
clc;

% Set up constants
V = 1e-4;       % Volume of each frustum
M = [0,0,1.3];  % Magnetisations in Teslas
d = 0.02;       % Distance between magnets

distances = 0.02;%0.001:0.001:0.06;
optH = zeros(length(distances),1);
optTheta = optH;
for k = 1:length(distances)
d = distances(k);

% Define varying parameters
theta = linspace(pi/4,pi/2,91);
theta = theta(2:end-1);
h = 0.02:0.001:0.033;

[Theta,H] = meshgrid(theta,h);
F = zeros(size(H));

% Set up for loops
for i = 1:length(h)
    for j = 1:length(theta)
        
        b = (sqrt(V/H(i,j)*tan(Theta(i,j))^2-H(i,j)^2/3)+H(i,j))/tan(Theta(i,j));
        bu = b-2*H(i,j)/tan(Theta(i,j));
        
        if b > 0 && b == real(b) && bu > 0 && bu == real(bu) % Make sure the frustum can physically exist
            verticesA = [-b/2,-b/2,0;-b/2,b/2,0;b/2,-b/2,0;b/2,b/2,0; ...
                -bu/2,-bu/2,-H(i,j);-bu/2,bu/2,-H(i,j);bu/2,-bu/2,-H(i,j);bu/2,bu/2,-H(i,j)];
            verticesB = -verticesA + repmat([0,0,d],size(verticesA,1),1);
            
            [myF,~,~] = polyhedronForce(verticesB,verticesA,M,M,12,mean(verticesB));
            F(i,j) = myF(3);
            
%             figure;
%             Sa = alphaShape(verticesA,inf);
%             Sb = alphaShape(verticesB,inf);
%             plot(Sa);
%             hold on;
%             plot(Sb);
        else
            F(i,j) = NaN;
        end
        
        
        
    end
end

figure;
surf(Theta*180/pi,H,F); xlabel('\theta'); ylabel('h'); zlabel('Force (N)');

[Max,I] = max(F(:));
fprintf('At a separation distance of %f cm, the maximum force is %f N, occuring when theta = %f degrees and h = %f cm\n',d*100,Max,Theta(I)*180/pi,H(I)*100);
optH(k) = H(I);
optTheta(k) = Theta(I)*180/pi;
end

figure;
yyaxis left
plot(distances*1000,optH*1000);
grid on
ylabel('Optimal height h (mm)');
yyaxis right
plot(distances*1000,optTheta);
title('Optimal height and angle to maximise force between two frusta');
ylabel('Optimal angle \theta (degrees)');
xlabel('Separation distance (mm)');