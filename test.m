
clc
clear

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
verticesB = verticesA + repmat([0.002,0.003,0],length(verticesA),1);
verticesB = verticesB + repmat([0,0,0.001-min(verticesB(:,3))+max(verticesA(:,3))],length(verticesB),1);
Sb = alphaShape(verticesB,Inf);
magB = [0,0,-1];

plot(Sa);
hold on;
plot(Sb);
grid on;

[Flinear,~,~,tlinear] = polyhedronForce(verticesA,verticesB,magA,magB,mean(verticesB),1e-30,20);
Flinear(1,:) = [];
% polyhedronForceQuadric(verticesA,verticesB,magA,magB,16)

x = 1:length(Flinear);
tquadric = zeros(length(x),1);
Fquadric = zeros(length(x),3);
tic;
for i = x
    FF = polyhedronForceQuadric(verticesA,verticesB,magA,magB,i);
    Fquadric(i,1:3) = FF;
    tquadric(i,1) = toc;
end

% 
% Fquadric
% figure;
% plot(x,Fquadric);
% title('Quadric method');
% legend('x','y','z');
% grid on;