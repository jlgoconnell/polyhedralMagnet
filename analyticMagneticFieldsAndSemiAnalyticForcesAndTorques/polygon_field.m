% This function takes in a set of points in a plane initpoints, a point of
% interest ptinterest, and a magnetisation vector M, and calculates the
% field of the polygon at the point of interest created by the
% anticlockwise-defined points.



function [Bx,By,Bz] = polygon_field(initpoints,ptinterest,M)

Bx = 0;
By = 0;
Bz = 0;

mu0 = 4*pi*10^(-7);

p = initpoints - ptinterest;
n = cross((p(:,2)-p(:,1)),(p(:,3)-p(:,1)));
n = n/norm(n);

% Calculate the necessary angles and rotational matrices
thetay = atan2(n(1),-n(3));
thetax = atan2(n(2),-sqrt(n(1)^2+n(3)^2));
Ry = [cos(thetay),0,sin(thetay);0,1,0;-sin(thetay),0,cos(thetay)];
Rx = [1,0,0;0,cos(thetax),-sin(thetax);0,sin(thetax),cos(thetax)];
R = Rx*Ry;
% Calculate new points (parallel to XY plane)
pp = R*p;
z_off = pp(3,1);

% Set up the vertices of the XY polygon:
myvertices = pp(1:2,:)';
numTri = size(myvertices,1);

for tri = 1:numTri
    x1 = myvertices(tri,1);
    x2 = myvertices(mod(tri,numTri)+1,1);
    y1 = myvertices(tri,2);
    y2 = myvertices(mod(tri,numTri)+1,2);

    theta = atan2((x1-x2),(y2-y1));
    X = [x1,y1;x2,y2]*[cos(theta),-sin(theta);sin(theta),cos(theta)];
    x1 = X(1,1);
    y1 = X(1,2);
    y2 = X(2,2);
    
    MdotN = myDot(M,n);
    if abs(MdotN) > 0
        B = mu0*myDot(M,n)/(4*pi)*b_field_tri(x1,0,z_off,y1,y2)*[cos(-theta),...
            -sin(-theta),0;sin(-theta),cos(-theta),0;0,0,1];
    else
       B = [0;0;0];
    end
    Bx = Bx + B(1);
    By = By + B(2);
    Bz = Bz + B(3);
    
end
    
B = [Bx,By,Bz]';
B = R^(-1)*B;

Bx = B(1);
By = B(2);
Bz = B(3);

end