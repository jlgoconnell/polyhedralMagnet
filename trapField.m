
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
x1 = vertices(1,1);
x2 = vertices(3,1);
m = [(vertices(3,2)-vertices(1,2))/(vertices(3,1)-vertices(1,1)),...
    (vertices(4,2)-vertices(2,2))/(vertices(4,1)-vertices(2,1))];
c = [vertices(1,2)-m(1)*vertices(1,1),vertices(2,2)-m(2)*vertices(2,1)];

x = obspt(:,1);
y = obspt(:,2);
z = obspt(:,3);
zd = vertices(1,3);


tic;
p = [1,2];

mn0 = m(m~=0);
cn0 = c(m~=0);
p = p(m~=0);
sqrtm = sqrt(1+mn0.^2);

S1 = sqrt((x-x1).^2+(y-mn0.*x1-cn0).^2+(z-zd).^2);
R1 = cn0.*mn0-x+x1+mn0.^2.*x1-mn0.*y+sqrtm.*S1;
I = [1,-1];
negS1 = S1.*I(m~=0);
M1 = (z-zd).^2.*(mn0.*(mn0.*x1+cn0-y)+x1-x+mn0.*negS1);
N1 = -(z-zd).*((cn0+mn0.*x-y).*(cn0+negS1+mn0.*x1-y)+(z-zd).^2);
C1 = (cn0+mn0.*x-y).^2.*(x1-x)+(z-zd).^2.*(mn0.^2.*(x-x1)+2*mn0.*(cn0-y+mn0.*x));
D1 = (z-zd).*(-(cn0-y).^2+2*mn0.*(x1-2*x).*(cn0-y)+mn0.^2.*x.*(-3*x+2*x1)+mn0.^2.*(z-zd).^2);
S2 = sqrt((x-x2).^2+(y-mn0.*x2-cn0).^2+(z-zd).^2);
R2 = cn0.*mn0-x+x2+mn0.^2.*x2-mn0.*y+sqrtm.*S2;
negS2 = S2.*I(m~=0);
M2 = (z-zd).^2.*(mn0.*(mn0.*x2+cn0-y)+x2-x+mn0.*negS2);
N2 = -(z-zd).*((cn0+mn0.*x-y).*(cn0+negS2+mn0.*x2-y)+(z-zd).^2);
C2 = (cn0+mn0.*x-y).^2.*(x2-x)+(z-zd).^2.*(mn0.^2.*(x-x2)+2*mn0.*(cn0-y+mn0.*x));
D2 = (z-zd).*(-(cn0-y).^2+2*mn0.*(x2-2*x).*(cn0-y)+mn0.^2.*x.*(-3*x+2*x2)+mn0.^2.*(z-zd).^2);

myBx = zeros(length(x),2);
myBx(:,m~=0) = 2*(-1).^p.*mn0./sqrtm.*log(R1./R2) + log((M1.^2+N1.^2).*(C2.^2+D2.^2)./((C1.^2+D1.^2).*(M2.^2+N2.^2)));
toc;





xq = [vertices(1,1),vertices(1,1),vertices(3,1),vertices(3,1)];
m = [m,m];
c = [c,c];


[x,m] = meshgrid(x,m);
[y,c] = meshgrid(y,c);
[z,xq] = meshgrid(z,xq);
x = x';
y = y';
z = z';
m = m';
c = c';
xq = xq';

p = repmat([1,2],size(m,1),2);
q = repelem([1,2],size(m,1),2);
mn0 = m~=0;
m0 = m==0;
pmn0 = reshape(p(mn0),size(m,1),[]);
qmn0 = reshape(q(mn0),size(m,1),[]);
pm0 = reshape(p(m0),size(m,1),[]);
qm0 = reshape(q(m0),size(m,1),[]);

% ymxc = y==m.*xq+c; % An edge case that occurs sometimes
S = sqrt((x-xq).^2+(y-m.*xq-c).^2+(z-zd).^2);
negS = S.*[1,-1,1,-1];
M = (z-zd).^2.*(m.*(m.*xq+c-y)+xq-x+m.*negS);
N = -(z-zd).*((c+m.*x-y).*(c+negS+m.*xq-y)+(z-zd).^2);
C = (c+m.*x-y).^2.*(xq-x)+(z-zd).^2.*(m.^2.*(x-xq)+2*m.*(c-y+m.*x));
D = (z-zd).*(-(c-y).^2+2*m.*(xq-2*x).*(c-y)+m.^2.*x.*(-3*x+2*xq)+m.^2.*(z-zd).^2);
R = c.*m-x+xq+m.^2.*xq-m.*y+sqrtm.*S;
myBx = (2*(-1).^(pmn0+qmn0+1).*reshape(m(mn0)./sqrtm(mn0).*log(R(mn0)),size(m,1),[])+(-1).^(qmn0).*reshape(log((C(mn0).^2+D(mn0).^2)./(M(mn0).^2+N(mn0).^2)),size(m,1),[])); % m ~= 0
% myBx = 2*[m(:,1)/sqrtm(:,1)*log(R(:,3)./R(:,1)),m(:,3)/sqrtm(:,3)*log(R(:,2)./R(:,4))]+log(((M(:,1:2).^2+N(:,1:2).^2).*(C(:,3:4).^2+D(:,3:4).^2))./((C(:,1:2).^2+D(:,1:2).^2).*(M(:,3:4).^2+N(:,3:4).^2)));
% myBx(ymxc) = 0;
myBx = [myBx,(-1).^(pm0+qm0).*reshape(log((c(m0)-y(m0)+S(m0))./(y(m0)-c(m0)+S(m0))),size(m,1),[])]; % m == 0
% myBx(isinf(myBx)) = 0;
myBx(x==xq) = 0;
assert(sum(sum(isinf(myBx))) == 0)
assert(sum(sum(isnan(myBx))) == 0)
Bx = MdotN/(8*pi)*sum(myBx,2);
myBy = -[-1,1,1,-1]./sqrtm.*log(-x+xq+m.*(c+m.*xq-y)+sqrtm.*S);
myBy(y==m.*xq+c) = 0;
assert(sum(sum(isinf(myBy))) == 0)
assert(sum(sum(isnan(myBy))) == 0)
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














