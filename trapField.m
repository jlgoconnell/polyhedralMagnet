
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

p = [1,2,1,2];
q = [1,1,2,2];
m = repmat(m,1,2);
c = repmat(c,1,2);
xq = [x1,x1,x2,x2];

S = sqrt((x-xq).^2+(y-m.*xq-c).^2+(z-zd).^2);
C = (xq-x).*(c+m.*x-y).^2+2*m.*(z-zd).^2.*(c+m.*x-y)-m.^2.*(z-zd).^2.*(xq-x);
D = -(z-zd).*((c+m.*x-y).*(c+m.*x-y-2*m.*(xq-x))-m.^2.*(z-zd).^2);
M = (c+m.*x-y).*(c-(-1).^p.*S+m.*xq-y)+(z-zd).^2;
N = (xq-x+m.*(c-(-1).^p.*S+m.*xq-y)).*(z-zd);
R = xq-x+m.*(c+m.*xq-y)+sqrt(1+m.^2).*S;

indxz = (x==xq)&(z==zd);
indxz = (abs(x-xq)<0.0001)&(abs(z-zd)<0.0001);
NUM = fliplr(c+m.*xq-y-abs(c+m.*xq-y))./abs(c+m.*xq-y)-(c+m.*xq-y+abs(m.*xq+c-y))./fliplr(abs(m.*xq+c-y));
DEN = 2*(c+m.*xq-y).*fliplr(c+m.*xq-y);
% NUM = (c(2)+m(2)*xq-y-abs(c(2)+m(2)*xq-y))./abs(c(1)+m(1)*xq-y)-(c(1)+m(1)*xq-y+abs(m(1)*xq+c(1)-y))./abs(m(2)*xq+c(2)-y);
% DEN = 2*(c(1)+m(1)*xq-y).*(c(2)+m(2)*xq-y);
M(indxz) = [sum(NUM(indxz),2),ones(length(NUM(indxz)),1)];
C(indxz) = [sum(DEN(indxz),2),ones(length(NUM(indxz)),1)];
% M(indxz) = [NUM(sum(indxz,2)~=0),ones(size(NUM(sum(indxz,2)~=0)))];
% C(indxz) = [DEN(sum(indxz,2)~=0),ones(size(DEN(sum(indxz,2)~=0)))];


% myBx = 2*(-1).^(p+q).*m./sqrt(1+m.^2).*log(R)+(-1).^q.*log((M.^2+N.^2)./((xq-x).^2+(z-zd).^2));
myBx = 2*(-1).^(p+q).*m./sqrt(1+m.^2).*log(R)+(-1).^q.*log((M.^2+N.^2)./(C.^2+D.^2));
myBy = (-1).^(p+q)./sqrt(1+m.^2).*log(R);
myBz = (-1).^(q).*atan2((M.*C+N.*D),(N.*C-M.*D));

Bx = -MdotN/(8*pi)*sum(myBx,2);
By = MdotN/(4*pi)*sum(myBy,2);
Bz = -MdotN/(4*pi)*sum(wrapToPi(myBz(:,1:2)+myBz(:,3:4)),2);

B = [Bx,By,Bz];

% Some bugtesting variables:
Y = m.*x+c-y;
Ys = m.*xq+c-y-(-1).^p.*S;
Z = z-zd;
X = xq-x;
A = (-1).^q.*log((M.^2+N.^2)./(C.^2+D.^2));

% assert(sum(isnan(B(:)))==0)
% B(isnan(B)) = 0;

Y1 = Y(:,1);
Y2 = Y(:,2);
Ys11 = Ys(:,1);
Ys21 = Ys(:,2);
X1 = X(:,1);

% NUM = (c(2)+m(2)*x-y-abs(c(2)+m(2)*x-y))./abs(c(1)+m(1)*x-y)-(c(1)+m(1)*x-y+abs(m(1)*x+c(1)-y))./abs(m(2)*x+c(2)-y);
% DEN = 2*(c(1)+m(1)*x-y).*(c(2)+m(2)*x-y);

end














