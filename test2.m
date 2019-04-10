


close all
clear
clc

tic

% Akoun and Yonnet:
cuboid = @(l,b,w) [l/2,b/2,w/2;l/2,-b/2,w/2;-l/2,b/2,w/2;-l/2,-b/2,w/2;l/2,b/2,-w/2;l/2,-b/2,-w/2;-l/2,b/2,-w/2;-l/2,-b/2,-w/2];
verticesA = cuboid(0.02,0.012,0.006);
verticesB = cuboid(0.012,0.02,0.006) + repmat([-0.004,-0.004,0.008],8,1);
magA = [0,0,0.38];
magB = magA;
torquepoint = mean(verticesB);
magnet_fixed.dim = [0.02 0.012 0.006];
magnet_float.dim = [0.012 0.02 0.006];
magnet_fixed.magn = 0.38;
magnet_float.magn = 0.38;
magnet_fixed.magdir = [0 0 1]; % z
magnet_float.magdir = [0 0 1]; % z
offset = [-0.004;-0.004;0.008];
displ = 0;
displ_range = offset+[1; 0; 0]*displ;
f1_xyz = magnetforces(magnet_fixed,magnet_float,displ_range);
t1_xyz = magnetforces(magnet_fixed,magnet_float,displ_range,'torque');

% Dodecahedra:
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
verticesB = verticesB + repmat([0,0,0.01-min(verticesB(:,3))+max(verticesA(:,3))],length(verticesB),1);
verticesB = verticesB*[cosd(36),-sind(36),0;sind(36),cosd(36),0;0,0,1];
Sb = alphaShape(verticesB,Inf);
magB = [0,0,-1];


[Fold,~,~,told] = polyhedronForce(verticesA,verticesB,magA,magB,mean(verticesB),1e-50,60);
Fold(1,:) = [];

i = 1;
ttotal = 0;
while ttotal(end) < told(end)
    [Fgauss(i,:),tgauss(i)] = polyhedronForceGauss(verticesA,verticesB,magA,magB,i);
	ttotal = cumsum(tgauss);
    i = i + 1;
end
tgauss = ttotal;

figure;
plot(told,Fold(:,3),'r.-',tgauss,Fgauss(:,3),'b.-');
legend('Old','New');
grid on;

errorold = abs(Fold-Fold(end,:));
errorgauss = abs(Fgauss-Fgauss(end,:));
figure;
subplot(1,3,1);
loglog(told,errorold(:,1),'r.-',tgauss,errorgauss(:,1),'b.-');
grid on;
legend('Old','New');
subplot(1,3,2);
loglog(told,errorold(:,2),'r.-',tgauss,errorgauss(:,2),'b.-');
grid on;
legend('Old','New');
subplot(1,3,3);
loglog(told,errorold(:,3),'r.-',tgauss,errorgauss(:,3),'b.-');
grid on;
legend('Old','New');

% figure;
% semilogy(1:length(errorold),errorold,'r.-',1:length(errorgauss),errorgauss,'b.-');


