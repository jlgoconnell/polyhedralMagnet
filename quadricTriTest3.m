
% Script to mesh magnet B and use quadric elements to calculate the force
% on it.
%
% James O'Connell 4th March 2019

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
subdivideNumber = 2;
Fac = minConvexHull(verticesB);
[Ver,~] = surfToMesh(verticesB(:,1),verticesB(:,2),verticesB(:,3));
norms = meshFaceNormals(Ver,Fac);

factri = triangulateFaces(Fac);
[ver,fac] = subdivideMesh(Ver,factri,subdivideNumber);
[verr,facc] = subdivideMesh(ver,fac,2);
B = polyhedronField(verticesA,magA,verr);

for i = 1%:length(Fac)
    
    n = norms(i,:);
    thetay = atan2(n(1),-n(3));
    thetax = atan2(n(2),-sqrt(n(1)^2+n(3)^2));
    Ry = [cos(thetay),0,sin(thetay);0,1,0;-sin(thetay),0,cos(thetay)];
    Rx = [1,0,0;0,cos(thetax),-sin(thetax);0,sin(thetax),cos(thetax)];
    Rn = Rx*Ry;
    R = Rn^-1;
    
    trifaceinfo = reshape(facc',12,numel(facc)/12)';
    faceinfo = zeros(size(trifaceinfo,1),6);
    for j = 1:size(trifaceinfo,1)
        faceinfo(j,1:6) = unique(trifaceinfo(j,:));
    end
    
    myver = verr*R;
%     xmat = [verr(:,1).^2,verr(:,1),verr(:,1).*verr(:,2),verr(:,2).^2
    
end

% for i = 1:length(Fac)
%     
%     n = norms(i,:)/norm(norms(i,:));
%     thetay = atan2(n(1),-n(3));
%     thetax = atan2(n(2),-sqrt(n(1)^2+n(3)^2));
%     Ry = [cos(thetay),0,sin(thetay);0,1,0;-sin(thetay),0,cos(thetay)];
%     Rx = [1,0,0;0,cos(thetax),-sin(thetax);0,sin(thetax),cos(thetax)];
%     Rn = Rx*Ry;
%     R = Rn^-1;
%     
%     [fac,~] = triangulateFaces(Fac{i});
%     [ver,fac2] = subdivideMesh(Ver,fac,subdivideNumber);
%     AA = reshape(fac2',12,numel(fac2)/12)';
%     for j = 1:size(AA,1)
%         BB(j,:) = unique(AA(j,:));
%     end
%     % B is a matrix containing information about all the field points
%     
%     thesepts = ver(BB,:);
%     B = polyhedronField(verticesA,magA,thesepts);
%     Brot = B*R;
%     
% end








