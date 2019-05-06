


function [F,T,t] = polyhedronForce(verticesA,verticesB,magA,magB,meshnum,torquept)

tic;

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

Bfield = polyhedronField(verticesA,magA,midpts);

F = sum(Bfield(fmids',:).*repelem(MdotN.*areas,3,1))/(12*pi*10^-7);
myleverpoint = midpts(fmids',:)-torquept;
T = sum(cross(myleverpoint,Bfield(fmids',:)).*repelem(MdotN.*areas,3,1))/(12*pi*10^-7);
t = toc;

end