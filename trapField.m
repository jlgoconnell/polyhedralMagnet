
% Function to calculate the magnetic field due to a charged trapezium plate
% assuming the plate is parallel to the XY axis.
%
% Inputs:
% vertices: A (4 x 3) matrix of vertices. Note that if a triangle needs to
% be calculated, this can have one vertex repeated. The vertices must be
% ordered from smaller x value to larger x value, then from smaller y to
% larger y when x values are equal.
% MdotN: A scalar representing the magnetic charge density in Teslas.
% obspt: An (n x 3) matrix of coordinates which the magnetic field will be
% calculated at.
% 
% Output:
% B: An (n x 3) matrix of the magnetic field at each obspt.
% 
% James O'Connell on 20th Feb 2019.

function B = trapField(vertices,MdotN,obspt)

% Get point and line information
xq = [vertices(1,1),vertices(1,1),vertices(3,1),vertices(3,1)];
m = [(vertices(3,2)-vertices(1,2))/(vertices(3,1)-vertices(1,1)),...
    (vertices(4,2)-vertices(2,2))/(vertices(4,1)-vertices(2,1))];
m = [m,m];
c = [vertices(1,2)-m(1)*vertices(1,1),vertices(2,2)-m(2)*vertices(2,1)];
c = [c,c];

% Run the equations
x = obspt(:,1);
y = obspt(:,2);
z = obspt(:,3);
zd = vertices(1,3);

[x,m] = meshgrid(x,m);
[y,c] = meshgrid(y,c);
[z,xq] = meshgrid(z,xq);
x = x';
y = y';
z = z';
m = m';
c = c';
xq = xq';
sqrtm = sqrt(1+m.^2);

edgecasecheck = (c==y) & (m==0); % An edge case that occurs sometimes
S = sqrt((x-xq).^2+(y-m.*xq-c).^2+(z-zd).^2);
negS = S.*[1,-1,1,-1];
M = (z-zd).^2.*(m.*(m.*xq+c-y)+xq-x+m.*negS);
N = -(z-zd).*((c+m.*x-y).*(c+negS+m.*xq-y)+(z-zd).^2);
C = (c+m.*x-y).^2.*(xq-x)+(z-zd).^2.*(m.^2.*(x-xq)+2*m.*(c-y+m.*x));
D = (z-zd).*(-(c-y).^2+2*m.*(xq-2*x).*(c-y)+m.^2.*x.*(-3*x+2*xq)+m.^2.*(z-zd).^2);
R = c.*m-x+xq+m.^2.*xq-m.*y+sqrtm.*S;
myBx = (2*[-1,1,1,-1].*m./sqrtm.*log(R)-[1,1,-1,-1].*log((C.^2+D.^2)./(M.^2+N.^2)));
% myBx = 2*[m(:,1)/sqrtm(:,1)*log(R(:,3)./R(:,1)),m(:,3)/sqrtm(:,3)*log(R(:,2)./R(:,4))]+log(((M(:,1:2).^2+N(:,1:2).^2).*(C(:,3:4).^2+D(:,3:4).^2))./((C(:,1:2).^2+D(:,1:2).^2).*(M(:,3:4).^2+N(:,3:4).^2)));
myBx(isinf(myBx)) = 0;
Bx = MdotN/(8*pi)*sum(myBx,2);
myBy = -[-1,1,1,-1]./sqrtm.*log(-x+xq+m.*(c+m.*xq-y)+sqrtm.*S);
By = MdotN/(4*pi)*sum(myBy,2);
rr = M.*C+N.*D;
x1 = rr(:,1:2);
x2 = rr(:,3:4);
rr = N.*C-M.*D;
y1 = rr(:,1:2);
y2 = rr(:,3:4);
myBz = -atan2(y1.*x2-y2.*x1,x1.*x2+y1.*y2);
Bz = MdotN/(4*pi)*(sum(myBz,2));

B = [Bx,By,Bz];

end














