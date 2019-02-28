
% Script to attempt a force calculation using vectorisation on a quadric
% surface element
%
% James O'Connell 28th Feb 2019

close all;
clear all;
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
magA = [0,0,1.3];
Sa = alphaShape(verticesA,Inf);

verticesB = verticesA + repmat([0,0,max(verticesA(:,3))-min(verticesA(:,3))+0.01],length(verticesA),1);
magB = [0,0,-1.3];
Sb = alphaShape(verticesB,Inf);

[Factual,Tactual] = polyhedronForce(verticesA,verticesB,magA,magB,mean(verticesB),1e-30,1);
Factual

% Define the facets
Fac = minConvexHull(verticesB);
[Ver,~] = surfToMesh(verticesB(:,1),verticesB(:,2),verticesB(:,3));
norms = meshFaceNormals(Ver,Fac);

subdivisionsize = 2;
for i = 1:size(norms,1)
    n = norms(i,:);
    MdotN = myDot(n,magB);
    thetay = atan2(n(1),-n(3));
    thetax = atan2(n(2),-sqrt(n(1)^2+n(3)^2));
    Ry = [cos(thetay),0,sin(thetay);0,1,0;-sin(thetay),0,cos(thetay)];
    Rx = [1,0,0;0,cos(thetax),-sin(thetax);0,sin(thetax),cos(thetax)];
    Rn = Rx*Ry;
    R = Rn^-1;
    
    [fac,~] = triangulateFaces(Fac{i});
    [ver,fac] = subdivideMesh(Ver,fac,2);
    
end