
% Function to calculate the magnetic field due to a charged trapezium plate
% assuming the plate is parallel to the XY axis.
%
% Inputs:
% vertices: A (4 x 3) matrix of vertices. Note that if a triangle needs to
% be calculated, this can have one vertex repeated. The vertices must be
% ordered from smaller x value to larger x value, then from smaller y to
% larger y when x values are equal.
% MdotN: A scalar representing the magnetic charge density in A/m.
% obspt: An (n x 3) matrix of coordinates which the magnetic field will be
% calculated at.
% 
% Output:
% B: An (n x 3) matrix of the magnetic field at each obspt.
% 
% James O'Connell on 20th Feb 2019.

function B = trapField(vertices,MdotN,obspt)

mu0MdotN = 4*pi*10^(-7)*MdotN;
myeps = sqrt(eps);

% Get point and line information
x1 = vertices(1,1);
x2 = vertices(3,1);
m = [(vertices(3,2)-vertices(1,2))/(vertices(3,1)-vertices(1,1)),...
    (vertices(4,2)-vertices(2,2))/(vertices(4,1)-vertices(2,1))];
c = [vertices(1,2)-m(1)*vertices(1,1),vertices(2,2)-m(2)*vertices(2,1)];

x = obspt(:,1);
y = obspt(:,2);
z = obspt(:,3);
zd = vertices(1,3);

% Set up some indexing stuff
p = [1,2,1,2];
q = [1,1,2,2];
m = repmat(m,1,2);
c = repmat(c,1,2);
xq = [x1,x1,x2,x2];
sqrt1plusm = sqrt(1+m.^2);

% Set up intermediate variables
X = xq-x;
Y = c+m.*xq-y;
Z = zd-z;
R = sqrt(X.^2+Y.^2+Z.^2);
S = X+m.*Y+sqrt1plusm.*R;
T = R+Y;
U = (m.*(X.^2+Z.^2)-X.*Y)./(Z.*R);

% Singularity treatment:
% If S = 0:
if any(abs(Z)<myeps) && any(any(abs(Y-m.*X)<myeps)) && any(any(X<0))
    indS = (abs(Z)<myeps) & (abs(Y-m.*X)<myeps) & (X<0);
    S(indS) = 1./R(indS);
end

% If T = 0:
if any(abs(Z)<myeps) && any(any(abs(X)<myeps)) && any(any(Y<0))
    indT = (abs(Z)<myeps) & (abs(X)<myeps) & (Y<0);
    T(indT) = 1./R(indT);
end

% If U = 0/0:
if any(abs(Z)<myeps) && (any(any(abs(X)<myeps)) || any(any(abs(Y-m.*X)<myeps)))
    indU = (abs(Z)<myeps) & ((abs(X)<myeps) | (abs(Y-m.*X)<myeps));
    U(indU) = 0;
end

% Calculate the field:
logS = log(S);
logT = log(T);
atanU = atan(U);
myBx = (-1).^(p+q).*(logT-m./sqrt1plusm.*logS);
myBy = (-1).^(p+q)./sqrt1plusm.*logS;
myBz = (-1).^(p+q).*atanU;

Bx = mu0MdotN/(4*pi)*sum(myBx,2);
By = mu0MdotN/(4*pi)*sum(myBy,2);
Bz = mu0MdotN/(4*pi)*sum(myBz,2);

% This makes the B-field correct on the surface of a magnet:
% Bz(abs(Z)<1000*eps) = abs(Bz(abs(Z)<1000*eps));

B = [Bx,By,Bz];

end














