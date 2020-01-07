
% Function to calculate the magnetic field across a square grid over a
% vertically-magnetised frustum magnet.

% James O'Connell 7th January 2020

function B = frustumField(l,L,h,magn,XY,Z,n)

vertices = [l/2,l/2,0;l/2,-l/2,0;-l/2,l/2,0;-l/2,-l/2,0;...
    L/2,L/2,-h;L/2,-L/2,-h;-L/2,L/2,-h;-L/2,-L/2,-h];

slanteddensity = (L-l)/2/sqrt(((L-l)/2)^2+h^2)*magn;

x = linspace(-XY/2,XY/2,n);
y = x';
[x,y] = meshgrid(x,y);
X = x(:);
Y = y(:);
obspt = [X,Y,repmat(Z,size(X))];

% Calculate the contributions from the surfaces:
theta = atan(2*h/(L-l));
R = [cos(theta),0,-sin(theta);0,1,0;sin(theta),0,cos(theta)];
slapts = [-L/2,-L/2,-h;-L/2,L/2,-h;-l/2,-l/2,0;-l/2,l/2,0]*R;
slaobs = obspt*R;

Btop = trapField(vertices([4,3,2,1],:),magn,obspt);
Bbot = trapField(vertices([8,7,6,5],:),-magn,obspt);
Bsla = trapField(slapts,slanteddensity,slaobs)*R^-1;
Bxsla = reshape(Bsla(:,1),size(x));
Bysla = reshape(Bsla(:,2),size(x));
Bzsla = reshape(Bsla(:,3),size(x));
% Rotate the coordinate systems:
Bxslaf = Bxsla + rot90(Bysla,1) - rot90(Bxsla,2) - rot90(Bysla,3);
Byslaf = Bysla - rot90(Bxsla,1) - rot90(Bysla,2) + rot90(Bxsla,3);
Bzslaf = Bzsla + rot90(Bzsla,1) + rot90(Bzsla,2) + rot90(Bzsla,3);
Bsla = [Bxslaf(:),Byslaf(:),Bzslaf(:)];


B = Btop + Bbot + Bsla;




end