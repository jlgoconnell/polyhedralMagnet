


function [F, t] = polyhedronForceGauss(verticesA,verticesB,magA,magB,meshnum)

tic;

f = minConvexHull(verticesB);
[v,~] = surfToMesh(verticesB(:,1),verticesB(:,2),verticesB(:,3));
[f,~] = triangulateFaces(f);
if meshnum > 1
    [v,f] = subdivideMesh(v,f,meshnum);
end
e = meshEdges(f);

% [v2,f2] = subdivideMesh(v,f,2);
MdotN = meshFaceNormals(v,f)*magB';
areas = meshFaceAreas(v,f);

% Try 2:
% m = max(f(:));
midpts = 0.5*(v(e(:,1),:)+v(e(:,2),:));
fmids = cell2mat(meshFaceEdges(v,e,f));

% a = f2';
% a = a(a > m);
% fmids = sort(reshape(a,9,numel(a)/9)',2);
% fmids = fmids(:,[1,4,7])-m;

Bfield = polyhedronField(verticesA,magA,midpts);

if sum(isnan(Bfield))>0
    fprintf('error');
end

F = sum(Bfield(fmids',:).*repelem(MdotN.*areas,3,1))/(12*pi*10^-7);

t = toc;

end