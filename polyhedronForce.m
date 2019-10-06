


function [F,T,t] = polyhedronForce(verticesA,verticesB,magA,magB,meshnum,torquept,d,varargin)

tic;

if nargin < 7
    d = [0,0,0];
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

    F(i,:) = sum(Bfield(fmids',:).*repelem(MdotN.*areas,3,1))/(12*pi*10^-7);
    myleverpoint = midptsd(fmids',:)-torquept;
    T(i,:) = sum(cross(myleverpoint,Bfield(fmids',:)).*repelem(MdotN.*areas,3,1))/(12*pi*10^-7);
    t(i) = toc;
end

end