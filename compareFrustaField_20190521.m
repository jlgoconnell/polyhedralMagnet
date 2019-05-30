
% Calculating the field of several frusta and cuboids
% 
% 21st May 2019

clc;
% close all;
clear;

cuboid = @(l,b,w) [l/2,b/2,w/2;l/2,-b/2,w/2;-l/2,b/2,w/2;-l/2,-b/2,w/2;l/2,b/2,-w/2;l/2,-b/2,-w/2;-l/2,b/2,-w/2;-l/2,-b/2,-w/2];

magA = [0,0,1];
magB = [0,0,1];
magC = [0,0,-1];

h = 0.02;
hcube = 0.02;
V = 12e-6;
myytheta = 90;
mytheta = deg2rad(myytheta);

zmin = 0.0;
zmax = 0.06;

zz = linspace(0,zmax,51);
zz = zz(2:end)';
% zz = linspace(0.015,0.025,41)';
dz = zz(2)-zz(1);
pts = [zeros(length(zz),2),zz];

Fz = zeros(length(zz),length(mytheta));

for i = 1:length(mytheta)
    theta = mytheta(i);
    
    h = 0.01;
%     if i == 1
%         V = V*1.35;
%     else
%         V = V/1.35;
%     end

    b = sqrt(V/h-(h/tan(theta))^2/3)+h/tan(theta);
    bu = b-2*h/tan(theta);
    hcube = V/b^2;
    verticesA = [-b/2,-b/2,0;-b/2,b/2,0;b/2,-b/2,0;b/2,b/2,0; ...
        -bu/2,-bu/2,-h;-bu/2,bu/2,-h;bu/2,-bu/2,-h;bu/2,bu/2,-h];
%     verticesC = [verticesA(:,1:2),-verticesA(:,3)]+repmat([0,0,zmax+hcube],length(verticesA),1);
    verticesC = [verticesA(:,1:2),-verticesA(:,3)]+repmat([0,0,-3*h],length(verticesA),1);
    scaleA = 0.65;
    scaleC = 2;
    verticesA = cuboid(scaleA*b,scaleA*b,h) + repmat([0,0,-0.5*h],length(verticesA),1);
    verticesC = cuboid(scaleC*b,scaleC*b,h) + repmat([0,0,-0.5*h],length(verticesA),1);
    
    for j = 1:length(zz)
        
        verticesB = cuboid(b,b,h) + repmat([0,0,hcube/2+zz(j)],length(verticesA),1);
        
        F1 = 2*polyhedronForce(verticesA,verticesB,magA,magB,16,mean(verticesB));
        F2 = polyhedronForce(verticesC,verticesB,magC,magB,16,mean(verticesB));
        Fz(j,i) = F1(3)+F2(3);
    end

    Sa = alphaShape(verticesA,inf);
    Sb = alphaShape(verticesB,inf);
    Sc = alphaShape(verticesC,inf);
%     figure;
%     plot(Sa);
%     hold on;
%     plot(Sb);
%     plot(Sc);
end

Fz = real(Fz)

% figure;
% plot(zz,Fz);
% legend(num2str(myytheta'));
% grid on;

K = zeros(size(Fz));
K(2:end-1,:) = (Fz(3:end,:)-Fz(1:end-2,:))/(2*dz);
K(1,:) = (Fz(2,:)-Fz(1,:))/dz;
K(end,:) = (Fz(end,:)-Fz(end-1,:))/dz;
% figure;
% plot(zz,K);
% legend(num2str(myytheta'));
% grid on;

% Knorm = K./Fz;
% figure;
% plot(zz,Knorm);
% grid on;
% legend(num2str(myytheta'));

figure;
subplot(2,2,[1,3]);
plot(Sc);
title([num2str(scaleA),', ',num2str(scaleC)]);
hold on;
plot(Sa,'FaceColor','k');
plot(Sb);
subplot(2,2,2);
plot(zz,Fz);
grid on;
subplot(2,2,4);
plot(zz,K);
grid on;

