
% Function to calculate the magnetic field of a tetrahedral permanent
% magnet at a given point.
%
% Inputs:
% vertices: A (4 x 3) matrix of coordinates of the tetrahedron vertices.
% mag: A (1 x 3) magnetisation vector in A/m.
% obspt: An (n x 3) matrix of coordinates for which the field is to be
% calculated at.
%
% Output:
% B: An (n x 3) matrix of the magnetic field for each obspt.
%
% James O'Connell 7th November 2019.


function B = tetField(vertices,mag,obspt)

% Set up variables
B = zeros(size(obspt));
vertices = vertices([1:4,1:3],:);

% Perform calcs for each face of the tet
for fac = 1:4
    facverts = vertices(fac:fac+2,:);
    
    % Calculate necessary rotation matrix
    z = cross(facverts(2,:)-facverts(1,:),facverts(3,:)-facverts(1,:));
    tempvec = vertices(fac+3,:)-vertices(fac+2,:);
    z = -sign(myDot(z,tempvec))*z/norm(z); % Make sure it's outward-facing
    MdotN = myDot(mag,z);
    if abs(MdotN) > eps
        AB = facverts(2,:)-facverts(1,:);
        AC = facverts(3,:)-facverts(1,:);
        thing = det([z;AC;AB]);
        if thing > 0
            y = AB;
        else
            y = AC;
            facverts = [facverts(1,:);facverts(3,:);facverts(2,:)];
        end
        y = y/norm(y);
        x = cross(y,z);
        R = [x',y',z'];

        % Solve field for each face
        trappts = round(facverts([1:3,3],:)*R,10);
        obsptr = obspt*R;
        Btemp = trapField(trappts,MdotN,obsptr);
        B = B + Btemp*R';
    end
end




end