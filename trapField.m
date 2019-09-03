
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

% % Set up intermediate variables
% S = sqrt((x-xq).^2+(y-m.*xq-c).^2+(z-zd).^2);
% C = (xq-x).*(c+m.*x-y).^2+2*m.*(z-zd).^2.*(c+m.*x-y)-m.^2.*(z-zd).^2.*(xq-x);
% D = -(z-zd).*((c+m.*x-y).*(c+m.*x-y-2*m.*(xq-x))-m.^2.*(z-zd).^2);
% M = (c+m.*x-y).*(c-(-1).^p.*S+m.*xq-y)+(z-zd).^2; % NI
% N = (xq-x+m.*(c-(-1).^p.*S+m.*xq-y)).*(z-zd); % NR
% R = xq-x+m.*(c+m.*xq-y)+sqrt(1+m.^2).*S;
% 
% % % Solve singularities:
% % indxz = (abs(x-xq)<eps)&(abs(z-zd)<eps);
% % YY = m.*x+c-y;
% % num = 2+zeros(size(indxz));
% % den = (-1).^(p+1).*YY-abs(YY)+fliplr(YY)+(-1).^(p+1).*abs(fliplr(YY));
% % M(indxz) = num(indxz);
% % N(indxz) = 0;
% % C(indxz) = den(indxz);
% % D(indxz) = 0;
% % 
% % indyz = (abs(y-m.*x-c)<eps)&(abs(z-zd)<eps);
% % NUM = ones(size(indyz));
% % DEN = ones(size(indyz));
% % M(indyz) = NUM(indyz);
% % N(indyz) = 0;
% % C(indyz) = 0;
% % D(indyz) = DEN(indyz);
% % 
% indR = abs(R)<eps;
% indRc = sum(indR,2)~=0;
% indRr = sum(indR,1)~=0;
% xx = x.*indR;
% mm = m.*indR;
% xqq = xq.*indR;
% R(indR) = (sqrt(1+mm(indR).^2)-mm(indR).^2)./((xx(indR)-xqq(indR)).*sqrt(1+mm(indR).^2));
% % 
% % indmy = (abs(y-c)<eps)&(abs(m)<eps);
% % zz = repmat(z,1,4);
% % C(indmy) = S(indmy);
% % D(indmy) = 0;
% 
% % Debugging variables:
% Y = c+m.*x-y;
% Yq = c+m.*xq-y;
% X = xq-x;
% Z = zd-z;
% Sp = (-1).^p.*S;
% ni = -(Y.*(S+Yq)+Z.^2).*Z;
% nr = -(X+m.*(S+Yq)).*Z.^2;
% di = -Z.*(Y.*(Y-2*m.*X)-m.^2.*Z.^2);
% dr = X.*Y.^2+2*m.*Z.^2.*Y-m.^2.*X.*Z.^2;
% 
% % Solve the equations
% myBx = 2*(-1).^(p+q).*m./sqrt(1+m.^2).*log(R)+(-1).^q.*log((M.^2+N.^2)./(C.^2+D.^2));
% myBy = (-1).^(p+q)./sqrt(1+m.^2).*log(R);
% myBz = (-1).^(q).*atan2((z-zd).*(M.*C+N.*D),(z-zd).*(N.*C-M.*D));
% 
% % Massive simplification to the x-field equation:
% % Note that this breaks if z-zd = 0 and x-xq = 0.
% myBx = (-1).^(p+q).*(m./sqrt(1+m.^2).*log(R)+log(S-c-m.*xq+y));
% 
% Bx = -MdotN/(4*pi)*sum(myBx,2);
% By = MdotN/(4*pi)*sum(myBy,2);
% Bz = -MdotN/(4*pi)*sum(wrapToPi(myBz(:,1:2)+myBz(:,3:4)),2);
% % Bz = -MdotN/(4*pi)*sum(mod(myBz(:,1:2)+myBz(:,3:4)+pi,2*pi)-pi,2);
% Bz = -MdotN/(4*pi)*sum(myBz,2);
% 
% B = [Bx,By,Bz];
% 
% % More debugging variables:
% MM = Z.^2.*(X+m.*(S+Yq));
% NN = -Z.*(Y.*(S+Yq)+Z.^2);
% CC = X.*Y.^2+2*m.*Y.*Z.^2-m.^2.*X.*Z.^2;
% DD = Z.*(m.^2.*Z.^2-Y.^2+2*m.*X.*Y);
% 
% newBz = (-1).^(p+q).*atan2(Z.*(X.*Y-m.*Z.^2),Z.*(Z.*S))-2*(-1).^q.*atan2(m.*Z,Y);

X = xq-x;
Y = c+m.*x-y;
Z = zd-z;
S = sqrt(X.^2+(Y+m.*X).^2+Z.^2);
R = X.*(1+m.^2)+m.*Y+sqrt(1+m.^2).*S;

myBx = -(-1).^(p+q).*(m./sqrt(1+m.^2).*log(R)+log(S-Y-m.*X));
myBy = (-1).^(p+q)./sqrt(1+m.^2).*log(R);
myBz = -(-1).^(p+q).*atan2(Z.*(X.*Y-m.*Z.^2),Z.*(Z.*S))+2*(-1).^q.*atan2(m.*Z,Y);

Bx = MdotN/(4*pi)*sum(myBx,2);
By = MdotN/(4*pi)*sum(myBy,2);
Bz = MdotN/(4*pi)*sum(myBz,2);

B = [Bx,By,Bz];

end














