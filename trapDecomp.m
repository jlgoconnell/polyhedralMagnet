
% Function to decompose a polygon in the XY plane into a set of trapezia
% with the parallel sides parallel to the Y-axis.
%
% James O'Connell 21st March 2019


function trapVertices = trapDecomp(vertices)

vertices = unique(sortrows(vertices),'rows');
trapVertices = vertices;

conlist = convhull(vertices);
xv = vertices(conlist,1);
yv = vertices(conlist,2);

for j = 1:length(xv)-1
    tempvec = (xv(1:end-1)-xv(j)).*(xv(2:end)-xv(j));
    inds = tempvec<0;
    
    if sum(inds) > 0
        intersecverts = [xv(tempvec<0),yv(tempvec<0);xv([false;tempvec<0]),yv([false;tempvec<0])];
        m = (intersecverts(2,2)-intersecverts(1,2))/(intersecverts(2,1)-intersecverts(1,1));
        c = intersecverts(1,2)-m*intersecverts(1,1);
        ypt = m*xv(j)+c;
        trapVertices = [trapVertices; xv(j), ypt];
    else if sum(tempvec==0) == 2
            trapVertices = [trapVertices; xv(j),yv(j)];
        end
    end
end

trapVertices = sortrows(trapVertices);

end