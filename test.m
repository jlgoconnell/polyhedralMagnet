
close all
clear
clc

tic

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

plot(Sa);
hold on;
plot(Sb);
grid on;

tic;
for i = 1:6
Ffft(i,1:3) = polyhedronForceFFT(verticesA,verticesB,magA,magB,i)';
% t = 23.66
tfft(i,1) = toc;
end

[Flinear,~,~,tlinear] = polyhedronForce(verticesA,verticesB,magA,magB,mean(verticesB),1e-30,tfft(end));
Flinear(1,:) = [];
% polyhedronForceQuadric(verticesA,verticesB,magA,magB,16)

x = 1:length(Flinear)
tquadric = zeros(length(x),1);
Fquadric = zeros(length(x),3);
tic;
for i = x
    FF = polyhedronForceQuadric(verticesA,verticesB,magA,magB,i);
    Fquadric(i,1:3) = FF;
    tquadric(i,1) = toc;
end

errorlinear = abs(Flinear-Flinear(end,:));
errorfft = abs(Ffft-Ffft(end,:));
errorquadric = abs(Fquadric-Fquadric(end,:));

figure;
plot(x,errorlinear(:,3),'r.-',x,errorquadric(:,3),'k.-');
grid on;
legend('Dumb','Quadric');
title('Error');
xlabel('Number of iterations');
ylabel('Error (N)');

figure;
semilogx(tlinear,Flinear(:,3),'r.-',tfft(3:end),Ffft(3:end,3),'b.-',tquadric,Fquadric(:,3),'k.-');
grid on;
legend('Dumb','FFT','Quadric');
title('z-Force');
xlabel('Time (s)');
ylabel('Force (N)');

figure;
loglog(tlinear,errorlinear(:,3),'r.-',tfft(3:end),errorfft(3:end,3),'b.-',tquadric,errorquadric(:,3),'k.-');
grid on;
legend('Dumb','FFT','Quadric');
title('Error');
xlabel('Time (s)');
ylabel('Error (N)');


toc


% figure;
% semilogx(tlinear,Flinear(:,3),'r.-',tquadric,Fquadric(:,3),'b.-');
% grid on;
% legend('Quadric','Dumb');
% title('z-Force between two magnets');
% xlabel('Time taken (seconds)');
% ylabel('Force (N)');
% 
% errorquadric = abs(Fquadric-Fquadric(end,:));
% errorlinear = abs(Flinear-Flinear(end,:));
% figure;
% semilogx(tlinear(3:end),errorlinear(3:end,3),'r.-',tquadric,errorquadric(:,3),'b.-');
% grid on;
% legend('Dumb','Quadric');
% title('Error in force calculation');
% xlabel('Time taken (seconds)');
% ylabel('Error (N)');
% 
% % 
% % Fquadric
% figure;
% plot(x,errorlinear(:,3),'r.-',x,errorquadric(:,3),'b.-');
% title('Error');
% legend('Dumb','Quadric');
% grid on;