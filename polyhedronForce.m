
% 

function [F,T] = polyhedronForce(verticesA,verticesB,magA,magB,torquepoint,tolerance,timeout,varargin)

% If missing arguments:
if nargin < 7
    timeout = Inf;
end
if nargin < 6
    tolerance = 1e-2;
end
if nargin < 5
    torquepoint = mean(verticesB);
end

tic;

% Initial conditions
F = [0,0,0];
T = [0,0,0];
Fold = 2*[tolerance,tolerance,tolerance];
Told = Fold;

% Choose a suitable starting mesh parameter
i = 1;

while max(abs([F,T]-[Fold,Told])) > tolerance && toc < timeout
    
    % Temp values to compare error for convergence
    Fold = F;
    Told = T;
    
    % Set up and subdivide mesh
    Fac = minConvexHull(verticesB);
    [Ver,~] = surfToMesh(verticesB(:,1),verticesB(:,2),verticesB(:,3));
    n = meshFaceNormals(Ver,Fac);
    MdotN = 1/(pi*4e-7)*dot(repmat(magB,size(n,1),1)',n')';
    Fac = Fac(abs(MdotN)>eps);
    [Fac,~] = triangulateFaces(Fac);
    [Ver,Fac] = subdivideMesh(Ver,Fac,2+i);
    n = meshFaceNormals(Ver,Fac);
    
    % Calculate field
    MdotN = 1/(pi*4e-7)*dot(repmat(magB,size(n,1),1)',n')';
    Fac = Fac(abs(MdotN)>1e-8,:);
    dA = meshFaceAreas(Ver,Fac);
    obspt = meshFaceCentroids(Ver,Fac);
    MdotN = MdotN(abs(MdotN)>1e-8,:);
    B = polyhedronField(verticesA,magA,obspt);
    
    % Calculate force and torque
    F = sum(B.*MdotN.*dA,1);
    myleverpoint = obspt-torquepoint;
    T = sum(MdotN.*cross(myleverpoint,B).*dA);
    
    i = i + 1
    
end

toc

end
