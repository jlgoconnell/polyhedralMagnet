

function [polygons,M] = addPolyhedron(points,magnetisation)

xc = mean(points);
V = [];
M = {};

% This code calculates all the planes which form exterior faces of the
% polyhedron:
S = alphaShape(points,Inf);
faces = boundaryFacets(S);
t1 = points(faces(:,1),:); % Vectors of points which form
t2 = points(faces(:,2),:); % the triangular faces of the
t3 = points(faces(:,3),:); % exterior faces
V = cross(t3-t1,t2-t1);
V = V./sqrt(V(:,1).^2+V(:,2).^2+V(:,3).^2);
V(:,4) = sum(V.*t1,2);
V = uniquetol(V,'ByRows',true); % A matrix representation of the plane equations


% This code sorts the points into each polygon face:
for plane = 1:length(V)
    polygons{plane}(1,:) = [1,1,1]; % Need to have one entry to make it work
    for p = 1:length(points)
        % If a point lies on a plane, add that point to a cell array
        if abs(myDot(V(plane,1:3),points(p,:)) - V(plane,4)) < 1e-6
            polygons{plane}(end+1,:) = points(p,:);
        end
    end
    polygons{plane} = polygons{plane}(2:end,:); % Delete that one entry
    
    % This code will sort the points into anticlockwise order:
    xc = mean(polygons{plane});
    ref = polygons{plane}(1,:)-xc;
    meas = polygons{plane}(:,:)-xc;
    theta = [];
    for a = 1:length(polygons{plane})
        theta(a) = acos(myDot(ref,meas(a,:))/(norm(ref)*norm(meas(a,:))));
        if sign(myDot(V(plane,1:3),cross(ref,meas(a,:)))) < 0
            theta(a) = -theta(a);
        end
    end
    theta = real(theta);
    [~,I] = sort(theta,'ascend');
    polygons{plane} = polygons{plane}(I,:)';
    M{plane} = -magnetisation;
    
end

end