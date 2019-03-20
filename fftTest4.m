
% Still testing FFT stuff
%
% James O'Connell 18th March 2019

% clear;
close all;
clc;

% Define the first magnet
phi = (1+sqrt(5))/2;
A = 0.01*phi;
verticesA = A*[1,1,1;1,1,-1;1,-1,1;1,-1,-1;-1,1,1;-1,1,-1;-1,-1,1;-1,-1,-1;...
    0,phi,1/phi;0,phi,-1/phi;0,-phi,1/phi;0,-phi,-1/phi;...
    1/phi,0,phi;1/phi,0,-phi;-1/phi,0,phi;-1/phi,0,-phi;...
    phi,1/phi,0;phi,-1/phi,0;-phi,1/phi,0;-phi,-1/phi,0];
thetax = 31.715*pi/180;
R_x = [1,0,0;0,cos(thetax),-sin(thetax);0,sin(thetax),cos(thetax)];
verticesA = verticesA*R_x;
Sa = alphaShape(verticesA,Inf);
magA = [0,0,1];

z = 0.02235;

% Define geometry area
x1 = 0.01;
x2 = 0.02;
y11 = -0.02;
y12 = -0.015;
y21 = 0.02;
y22 = 0.01;
ys = [y11,y12,y21,y22];
ymin = min(ys);

% Calculate field at all points
n = 64;
x = linspace(x1,x2,n);
y = linspace(min(ys),max(ys),n);
[X,Y] = meshgrid(x,y);
B = polyhedronField(verticesA,magA,[X(:),Y(:),repmat(z,size(X(:)))]);
Bz = reshape(B(:,3),size(X));
f = Bz;

% Calculate FFT
ff = fft2(f);

% Shift frequencies appropriately
freqsx = 1/(x(2)-x(1))*(0:(n-1))/n;
freqsy = 1/(y(2)-y(1))*(0:(n-1))/n;
freqsx = fftshift(freqsx);
freqsy = fftshift(freqsy);
freqsx = freqsx - freqsx(1);
freqsy = freqsy - freqsy(1);

[Freqsx,Freqsy] = meshgrid(freqsx,freqsy);
A = 2*pi*1i*Freqsx;
B = 2*pi*1i*Freqsy;

m1 = (y12-y11)/(x2-x1);
m2 = (y22-y21)/(x2-x1);
if m1 == 0
    m1 = 1*eps;
end
if m2 == 0
    m2 = 1*eps;
end
c1 = y11-m1*x1;
c2 = y21-m2*x1;

% Define grid
nn = 201;
x = repmat(linspace(x1,x2,nn),nn,1);
ms = repmat(linspace(m1,m2,nn)',1,nn);
cs = repmat(linspace(c1,c2,nn)',1,nn);
y = ms.*x+cs;

% Solve integral equations
F = zeros(size(A));
F(1,1) = ff(1,1)*((m2-m1)/2*(x2^2-x1^2)+(c2-c1)*(x2-x1));
F(1,2:end) = ff(1,2:end)./A(1,2:end).*(exp(A(1,2:end)*(x2-x1)).* ...
    (y22-y12-(m2-m1)./A(1,2:end))-y21+y11+(m2-m1)./A(1,2:end));
F(2:end,1) = ff(2:end,1)./(B(2:end,1).^2.*exp(B(2:end,1)*ymin)).*((exp(B(2:end,1)* ...
    y22)-exp(B(2:end,1)*y21))/m2 - (exp(B(2:end,1)*y12)-exp(B(2:end,1) ...
    *y11))/m1);
F(2:end,2:end) = ff(2:end,2:end)./(B(2:end,2:end).* ...
    exp(A(2:end,2:end)*x1+B(2:end,2:end)*ymin)).* ...
    ((exp(A(2:end,2:end)*x2+B(2:end,2:end)*y22)-exp(A(2:end,2:end)*x1+...
    B(2:end,2:end)*y21))./(A(2:end,2:end)+B(2:end,2:end)*m2) - ...
    (exp(A(2:end,2:end)*x2+B(2:end,2:end)*y12)-exp(A(2:end,2:end)*x1+...
    B(2:end,2:end)*y11))./(A(2:end,2:end)+B(2:end,2:end)*m1));
F(isnan(F)) = 0;
inttotal = sum(F(:))/n^2;

% Estimate actual integral
midptsx = 0.25*(x(1:end-1,1:end-1)+x(1:end-1,2:end)+x(2:end,1:end-1)+x(2:end,2:end));
midptsy = 0.25*(y(1:end-1,1:end-1)+y(1:end-1,2:end)+y(2:end,1:end-1)+y(2:end,2:end));
BB = polyhedronField(verticesA,magA,[midptsx(:),midptsy(:),repmat(z,size(midptsx(:)))]);
Bz = reshape(BB(:,3),size(midptsx));
dx = x(1:end-1,2:end)-x(1:end-1,1:end-1);
dy = y(2:end,:)-y(1:end-1,:);
area = 0.5*(dy(:,2:end)+dy(:,1:end-1)).*dx;

actual = sum(sum(Bz.*area))
estimate = real(inttotal)
pcerror = (estimate-actual)/actual*100





