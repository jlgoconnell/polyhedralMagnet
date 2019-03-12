
% Script to mesh magnet B and use quadric elements to calculate the force
% on it.
%
% James O'Connell 5th March 2019

close all;
clear all;
clc;

phi = (1+sqrt(5))/2;
AA = 0.01*phi;
verticesA = AA*[1,1,1;1,1,-1;1,-1,1;1,-1,-1;-1,1,1;-1,1,-1;-1,-1,1;-1,-1,-1;...
    0,phi,1/phi;0,phi,-1/phi;0,-phi,1/phi;0,-phi,-1/phi;...
    1/phi,0,phi;1/phi,0,-phi;-1/phi,0,phi;-1/phi,0,-phi;...
    phi,1/phi,0;phi,-1/phi,0;-phi,1/phi,0;-phi,-1/phi,0];
thetax = 31.715*pi/180;
R_x = [1,0,0;0,cos(thetax),-sin(thetax);0,sin(thetax),cos(thetax)];
verticesA = verticesA*R_x;
magA = [0,0,1.3];
Sa = alphaShape(verticesA,Inf);

verticesB = verticesA + repmat([0,0,max(verticesA(:,3))-min(verticesA(:,3))+0.01],length(verticesA),1);
magB = [0,0,-1.3];
Sb = alphaShape(verticesB,Inf);

% Decompose into polygon faces
subdivideNumber = 8;
Fac = minConvexHull(verticesB);
[Ver,~] = surfToMesh(verticesB(:,1),verticesB(:,2),verticesB(:,3));
norms = meshFaceNormals(Ver,Fac);

fac = triangulateFaces(Fac);
% [Ver,fac] = subdivideMesh(Ver,fac,subdivideNumber);

edges = fac(:,[1,2,2,3,3,1])';
edges = reshape(edges,2,numel(edges)/2)';
facc = fac(:);

for i = 1:size(edges,1)
    Ver(end+1,:) = 0.5*(Ver(edges(i,1),:)+Ver(edges(i,2),:));
    facc(end+1) = size(Ver,1);
%     [ver,ia,ic] = unique(Ver,'rows');
    ver = Ver;
end
% a = 1:size(Ver,1);
% facc = a(facc)';
facc = reshape(facc,numel(facc)/6,6);

B = polyhedronField(verticesA,magA,ver);

vv = Ver(facc,:);
plot3(vv(:,1),vv(:,2),vv(:,3),'go');
hold on;
drawMesh(Ver,fac);

for i = 1:length(Fac)
    
    n = norms(i,:);
    thetay = atan2(n(1),-n(3));
    thetax = atan2(n(2),-sqrt(n(1)^2+n(3)^2));
    Ry = [cos(thetay),0,sin(thetay);0,1,0;-sin(thetay),0,cos(thetay)];
    Rx = [1,0,0;0,cos(thetax),-sin(thetax);0,sin(thetax),cos(thetax)];
    Rn = Rx*Ry;
    R = Rn^-1;
    
    
    
end








