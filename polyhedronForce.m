
% Function to estimate the force and torque between two polyhedral
% permanent magnets, both with constant uniform magnetisation. Let one
% magnet be magnet A and the other be magnet B. This function calculates
% the force and torque on magnet B.
%
% Inputs:
% verticesA: The vertices of one magnet in an (n x 3) matrix.
% verticesB: The vertices of the magnet we are calculating the force and
% torque on in an (n x 3) matrix.
% magA: The magnetisation vector of magnet A in Teslas.
% magB: The magnetisation vector of magnet B in Teslas.
% torquepoint: The point about which the torque on magnet B is calculated.
% tolerance: The maximum difference between the torque and force before
% convergence is reached.
% timeout: The maximum time this function is allowed to run for. After
% this, the function will finish the current iteration then return the
% results for the current iteration.
%
% Outputs:
% F: A vector representing the x, y, and z forces on magnet B.
% T: A vector representing the x, y, and z torques on magnet B.

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
i = 10;
ctr = 0;
oldC = Inf;
C = 0;

while max(abs(C-oldC)) > tolerance && toc < timeout
    
    % Temp values to compare error for convergence
    Fold = F;
    Told = T;
    oldC = C;
    
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
    
    i = i + 1;
    ctr = ctr + 1;
    
    C = [(exp(1)*F-Fold)/(exp(1)-1),(exp(1)*T-Told)/(exp(1)-1)];
    
%     if i > 11
%     figure(1);
%     plot(i-1:i,[Fold(1),F(1)],'r');
%     hold on;
%     figure(2);
%     plot(i-1:i,[Fold(2),F(2)],'k');
%     hold on;
%     figure(3);
%     plot(i-1:i,[Fold(3),F(3)],'b');
%     hold on;
%     end
    
end

F = C(1:3);
T = C(4:6);
ctr

toc

end
