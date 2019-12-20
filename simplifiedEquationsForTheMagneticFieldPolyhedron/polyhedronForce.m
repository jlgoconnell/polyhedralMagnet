
% Function to calculate the magnetic force and torque between two magnets,
% magnet A and magnet B with magnetisation strengths magA and magB.
%
% Inputs:
% verticesA: An (n x 3) matrix of the vertices of magnet A.
% verticesB: An (n x 3) matrix of the vertices of magnet B.
% magA: A (1 x 3) vector describing the magnetisation of magnet A in A/m.
% magB: A (1 x 3) vector describing the magnetisation of magnet B in A/m.
% meshnum: A scalar describing how fine the mesh should be on magnet B. A
% larger number will give a finer mesh, leading to a more accurate but
% slower calculation. A mesh number between 5 and 10 will usually give a
% fast and fairly accurate calculation.
% torquept: A (1 x 3) vector describing the point about which the torque on
% magnet B should be calculated.
% d: An (n x 3) matrix of displacement vectors over which the force and
% torque should be calculated.
%
% Output:
% F: An (n x 3) matrix of the magnetic force on magnet B at each
% displacement.
% T: An (n x 3) matrix of the magentic torque on magnet B at each
% displacement.
%
% James O'Connell 20th Feb 2019.


function [F,T,t] = polyhedronForce(verticesA,verticesB,magA,magB,meshnum,torquept,d,varargin)

tic;

if nargin < 7
    d = [0,0,0];
end

if nargin < 6
    torquept = mean(verticesB);
end

f = minConvexHull(verticesB);
[v,~] = surfToMesh(verticesB(:,1),verticesB(:,2),verticesB(:,3));
[f,~] = triangulateFaces(f);
MdotN = meshFaceNormals(v,f)*magB';
f = f(abs(MdotN)>eps,:);
if meshnum > 1
    [v,f] = subdivideMesh(v,f,meshnum);
end

e = zeros(numel(f), 2);
for i = 1:length(f)
    ff = f(i, :);
    e(3*i-2:3*i, :) = [ff' ff([2:end 1])'];
end
[e,~,IE] = unique(sort(e, 2),'rows');

MdotN = meshFaceNormals(v,f)*magB';
areas = meshFaceAreas(v,f);

midpts = 0.5*(v(e(:,1),:)+v(e(:,2),:));
fmids = reshape(IE,3,numel(IE)/3)';

% Work out the convex hull of magnet A:
if iscell(verticesA)
    for i = 1:length(verticesA)
        Fac{i} = minConvexHull(verticesA{i});
    end
else
    Fac = minConvexHull(verticesA);
end

% Loop through all displacements and calculate field, force, and torque:
for i = 1:size(d,1)
    midptsd = midpts + repmat(d(i,:),size(midpts,1),1);
    
    if iscell(verticesA)
        Bfield = [0,0,0];
        for j = 1:length(verticesA)
            Bfieldtemp = polyhedronField(verticesA{j},magA{j},midptsd,Fac{j});
            Bfield = Bfield + Bfieldtemp;
        end
    else
        Bfield = polyhedronField(verticesA,magA,midptsd,Fac);
    end

    F(i,:) = sum(Bfield(fmids',:).*repelem(MdotN.*areas,3,1))/3;
    myleverpoint = midptsd(fmids',:)-torquept;
    T(i,:) = sum(cross(myleverpoint,Bfield(fmids',:)).*repelem(MdotN.*areas,3,1))/3;
    t(i) = toc;
end

end