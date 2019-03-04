
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

function [F,T,C,t] = polyhedronForce(verticesA,verticesB,magA,magB,torquepoint,tolerance,timeout,varargin)

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
t = [];

% Choose a suitable starting mesh parameter
i = 1;
ctr = 0;
oldC = 0;
C = [0,0,0,0,0,0;Inf,Inf,Inf,Inf,Inf,Inf];

% figure;

while max(abs(C(end,:)-C(end-1,:))) > tolerance && toc < timeout
%     max(abs(C(end,:)-[F(end,:),T(end,:)]))
    % Temp values to compare error for convergence
%     Fold = F(end,:);
%     Told = T(end,:);
%     oldC = C(end,:);
    
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
    F(end+1,:) = sum(B.*MdotN.*dA,1);
    myleverpoint = obspt-torquepoint;
    T(end+1,:) = sum(MdotN.*cross(myleverpoint,B).*dA);
    
    i = i + 1;
    ctr = ctr + 1;
    
    if ctr >= 2
        C(end+1,:) = [(F(end-2,:).*F(end,:)-F(end-1,:).^2)./(F(end-2,:)-2*F(end-1,:)+F(end,:)),(T(end-2,:).*T(end,:)-T(end-1,:).^2)./(T(end-2,:)-2*T(end-1,:)+T(end,:))];
    else
        C(end+1,1:6) = zeros(1,6);
    end
    
    
    
    t(end+1,1) = toc;
    
%     if ctr > 1
%         subplot(2,3,1);
%         plot(ctr-1:ctr,[F(end-1,1),F(end,1)],'r');
%         hold on;
%         plot(ctr,D(end,1),'r.');
%         grid on;
%         ylabel('Fx');
%         subplot(2,3,2);
%         plot(ctr-1:ctr,[F(end-1,2),F(end,2)],'k');
%         hold on;
%         plot(ctr,D(end,2),'k.');
%         grid on;
%         ylabel('Fy');
%         subplot(2,3,3);
%         plot(ctr-1:ctr,[F(end-1,3),F(end,3)],'b');
%         hold on;
%         plot(ctr,D(end,3),'b.');
%         grid on;
%         ylabel('Fz');
%         subplot(2,3,4);
%         plot(ctr-1:ctr,[T(end-1,1),T(end,1)],'r');
%         hold on;
%         plot(ctr,D(end,4),'r.');
%         grid on;
%         ylabel('Tx');
%         subplot(2,3,5);
%         plot(ctr-1:ctr,[T(end-1,2),T(end,2)],'k');
%         hold on;
%         plot(ctr,D(end,5),'k.');
%         grid on;
%         ylabel('Ty');
%         subplot(2,3,6);
%         plot(ctr-1:ctr,[T(end-1,3),T(end,3)],'b');
%         hold on;
%         plot(ctr,D(end,6),'b.');
%         grid on;
%         ylabel('Tz');
%     end
% max(abs(C(end,:)-C(end-1,:)))
    
end

C = C(3:end,:);
% 
figure;
drawMesh(Ver,Fac,'white')

% max(abs([F(end,:),T(end,:)]-[Fold,Told]))

% F = F(2:end,:);
% T = T(2:end,:);

% F = C(end,1:3);
% T = C(end,4:6);

% F = F(2:end,:);
% T = T(2:end,:);
% ctr

% ctr

toc

% max(abs(C(end,:)-[F(end,:),T(end,:)]))

end
