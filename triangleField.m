
% Function to calculate the magnetic field of a triangular surface with
% magnetic surface charge sigma at a given set of points.

function B = triangleField(vertices,sigma,obspt)

B = zeros(size(obspt));

% Calculate rotation matrix
pt1 = vertices(1,:);
pt2 = vertices(2,:);
pt3 = vertices(3,:);
z = cross(pt2-pt1,pt3-pt1);
z = z/norm(z);
AB = pt2-pt1;
AC = pt3-pt1;
if det([z;AC;AB]) > 0
    y = AB;
else
    y = AC;
    vertices = [vertices(1,:);vertices(3,:);vertices(2,:)];
end
y = y/norm(y);
x = cross(y,z);
R = [x',y',z'];

vertices = [vertices;vertices(3,:)];
rotverts = vertices*R;
rotobs = obspt*R;

B = trapField(rotverts,sigma,rotobs);
B = B*R';

end