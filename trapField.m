
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

% Set up some indexing stuff
p = [1,2,1,2];
q = [1,1,2,2];
m = repmat(m,1,2);
c = repmat(c,1,2);
xq = [x1,x1,x2,x2];

% Set up intermediate variables
S = sqrt((x-xq).^2+(y-m.*xq-c).^2+(z-zd).^2);
C = (xq-x).*(c+m.*x-y).^2+2*m.*(z-zd).^2.*(c+m.*x-y)-m.^2.*(z-zd).^2.*(xq-x);
D = -(z-zd).*((c+m.*x-y).*(c+m.*x-y-2*m.*(xq-x))-m.^2.*(z-zd).^2);
M = (c+m.*x-y).*(c-(-1).^p.*S+m.*xq-y)+(z-zd).^2;
N = (xq-x+m.*(c-(-1).^p.*S+m.*xq-y)).*(z-zd);
R = xq-x+m.*(c+m.*xq-y)+sqrt(1+m.^2).*S;

% Modify intermediate variables for singular cases
% indm = repmat(m,size(x))==0;
% mm = m.*indm;
% cc = c.*indm;
% yy = y.*indm;
% xx = x.*indm;
% zz = z.*indm;
% xqq = xq.*indm;
% R(indm) = (S(indm)+cc(indm)-yy(indm))./sqrt((xqq(indm)-xx(indm)).^2+(zz(indm)-zd).^2);
% M(indm) = 1;
% N(indm) = 0;
% C(indm) = 1;
% D(indm) = 0;

indxz = (abs(x-xq)<eps)&(abs(z-zd)<eps);
NUM = [(c(2)+m(2)*x-y-abs(c(2)+m(2)*x-y))./abs(c(1)+m(1)*x-y)-(c(1)+m(1)*x-y+abs(m(1)*x+c(1)-y))./abs(m(2)*x+c(2)-y),ones(size(x))];
DEN = [2*(c(1)+m(1)*x-y).*(c(2)+m(2)*x-y),ones(size(x))];
NUM = [NUM,NUM];
DEN = [DEN,DEN];
M(indxz) = NUM(indxz);
C(indxz) = DEN(indxz);
YY = m.*x+c-y;
num = 2+zeros(size(indxz));%2*sign(YY).^p;
den = (-1).^(p+1).*YY-abs(YY)+fliplr(YY)+(-1).^(p+1).*abs(fliplr(YY));
MM = M;
CC = C;
M(indxz) = num(indxz);
C(indxz) = den(indxz);

indyz = (abs(y-m.*x-c)<eps)&(abs(z-zd)<eps);
NUM = ones(size(indyz));
DEN = ones(size(indyz));
M(indyz) = NUM(indyz);
C(indyz) = DEN(indyz);

indR = abs(R)<eps;
xx = x.*indR;
mm = m.*indR;
xqq = xq.*indR;
R(indR) = (sqrt(1+mm(indR).^2)-mm(indR).^2)./((xx(indR)-xqq(indR)).*sqrt(1+mm(indR).^2));

indmy = (abs(y-c)<eps)&(abs(m)<eps);
zz = repmat(z,1,4);
% C(indmy) = (xqq(indmy)-xx(indmy)).^2+(zz(indmy)-zd).^2;
C(indmy) = S(indmy);
D(indmy) = 0;

% Solve the equations
myBx = 2*(-1).^(p+q).*m./sqrt(1+m.^2).*log(R)+(-1).^q.*log((M.^2+N.^2)./(C.^2+D.^2));
myBy = (-1).^(p+q)./sqrt(1+m.^2).*log(R);
myBz = (-1).^(q).*atan2((M.*C+N.*D),(N.*C-M.*D));

Bx = -MdotN/(8*pi)*sum(myBx,2);
By = MdotN/(4*pi)*sum(myBy,2);
Bz = -MdotN/(4*pi)*sum(wrapToPi(myBz(:,1:2)+myBz(:,3:4)),2);

B = [Bx,By,Bz];

end














