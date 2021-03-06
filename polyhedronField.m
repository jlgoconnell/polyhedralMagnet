
% Function to calculate the magnetic field at the point obspt due to a
% polyhedral permanent magnet with a list of vertices and a given
% magnetisation vector mag.
%
% Inputs:
% vertices: An (n x 3) matrix of coordinates of the polyhedron vertices.
% mag: A (1 x 3) magnetisation vector in A/m.
% obspt: An (n x 3) matrix of coordinates for which the field is to be
% calculated at.
%
% Output:
% B: An (n x 3) matrix of the magnetic field for each obspt.
%
% James O'Connell 20th Feb 2019.

function B = polyhedronField(vertices, mag, obspt, Fac, varargin)

B = zeros(size(obspt));

% Get face and vertex information of the polyhedron
if nargin < 4
    Fac = minConvexHull(vertices);
end
[Ver,~] = surfToMesh(vertices(:,1),vertices(:,2),vertices(:,3));
norms = meshFaceNormals(Ver,Fac);
MdotN = dot(repmat(mag,size(norms,1),1)',norms')';
Fac = Fac(abs(MdotN)>eps);
norms = meshFaceNormals(Ver,Fac);

% Rotate each face and split into trapezia to calculate field
for i = 1:length(Fac)
    
    facepts = Ver(Fac{i},:);
    facepts = round(facepts,8);
    
    % Rotate the face to be parallel to the XY plane
    n = norms(i,:)/norm(norms(i,:));
    m = sqrt(n(2)^2+n(3)^2);
    sg = 1;%(n(3)>0)*2-1; % sign() function such that sign(0) = 1
    if m < eps % Special case when a y-rotation of 90deg is needed
        R = [0,0,1;0,1,0;-1,0,0];
    else
        R = [m,            0,          sg*n(1);...
             -n(1)*n(2)/m, sg*n(3)/m,  sg*n(2);...
             -n(1)*n(3)/m, -sg*n(2)/m, sg*n(3)];
    end
    Rn = R';
    facepts = facepts*R;
    thisz = facepts(1,3);
    
    % Decompose into trapezia
    trappts = trapDecomp(facepts(:,1:2));
    trappts = [trappts, facepts(1,3)*ones(length(trappts),1)];
    
    % Now calculate field of each trapezium
    Bxpoly = zeros(size(obspt,1),1);
    Bypoly = Bxpoly;
    Bzpoly = Bxpoly;
    MdotN = myDot(mag,n);
    
    for j = 1:length(trappts)/2-1
        thisobspt = obspt*R;
        thesepts = trappts(2*j-1:2*j+2,:);
        Btrap = trapField(thesepts,MdotN,thisobspt);
        Bxpoly = Bxpoly + Btrap(:,1);
        Bypoly = Bypoly + Btrap(:,2);
        Bzpoly = Bzpoly + Btrap(:,3);
    end
    
    B = B + [Bxpoly,Bypoly,Bzpoly]*Rn;
    
%     assert(sum(isnan(B(:)))==0);
    
end

end










