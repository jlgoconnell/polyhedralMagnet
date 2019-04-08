
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
verticesB = verticesA + repmat([0.00,0.00,0],length(verticesA),1);
verticesB = verticesB + repmat([0,0,0.001-min(verticesB(:,3))+max(verticesA(:,3))],length(verticesB),1);
verticesB = verticesB*[cosd(36),-sind(36),0;sind(36),cosd(36),0;0,0,1];
Sb = alphaShape(verticesB,Inf);
magB = [0,0,-1];
torquept = mean(verticesB);

tic;
for i = 1:7
Ffft(i,1:3) = polyhedronForceFFT(verticesA,verticesB,magA,magB,i)';
% t = 23.66
tfft(i,1) = toc;
end

tic;
[Fold,~,~,told] = polyhedronForce(verticesA,verticesB,magA,magB,torquept,1e-5,tfft(end))
toc;
Fold = Fold(2:end,:);

figure;
semilogx(tfft(3:end),Ffft(3:end,3),'r.-',told,Fold(:,3),'b.-');
grid on;
legend('FFT','Dumb');
title('z-Force between two magnets');
xlabel('Time taken (seconds)');
ylabel('Force (N)');

errorfft = abs(Ffft-Ffft(end,:));
errorold = abs(Fold-Fold(end,:));
figure;
semilogx(tfft(3:end),errorfft(3:end,3),'r.-',told,errorold(:,3),'b.-');
grid on;
legend('FFT','Dumb');
title('Error in force calculation');
xlabel('Time taken (seconds)');
ylabel('Error (N)');

Freal = [28.5137;42.7820;311.2703];

Flinear = Fold;
tlinear = told;

% mean(abs(Freal-Ffft))
% mean(abs(Freal-Fold(end,:)'))
% 
% for i = 3:8
%     tic;
%     Ffft(:,i) = polyhedronForceFFT(verticesA,verticesB,magA,magB,i);
%     t(i) = toc;
% end