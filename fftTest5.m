
% Still testing FFT stuff
%
% James O'Connell 18th March 2019

clear;
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

z = 0.0225;

% Define geometry area
x1 = 0.01;
x2 = 0.02;
y11 = -0.02;
y12 = -0.015;
y21 = 0.02;
y22 = 0.01;
% y22 = y21;
% y12 = y11;
ys = [y11,y12,y21,y22];
ymin = min(ys);
ymax = max(ys);

% Calculate field at all points
n = 64;
deltad = 0.1;
xstart = x1-deltad*(x2-x1);
xend = x2+deltad*(x2-x1);
ystart = ymin-deltad*(ymax-ymin);
yend = ymax+deltad*(ymax-ymin);
x = linspace(xstart,xend,n);
y = linspace(ystart,yend,n);
[X,Y] = meshgrid(x,y);
B = polyhedronField(verticesA,magA,[X(:),Y(:),repmat(z,size(X(:)))]);
Bz = reshape(B(:,3),size(X));
f = Bz;

figure;
surf(X,Y,Bz);

% Calculate FFT
tic;
ff = fft2(f);
toc;

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
c1 = y11-m1*x1;
c2 = y21-m2*x1;

% Define grid
nn = 51;
x = repmat(linspace(x1,x2,nn),nn,1);
ms = repmat(linspace(m1,m2,nn)',1,nn);
cs = repmat(linspace(c1,c2,nn)',1,nn);
y = ms.*x+cs;

tic;
% Solve integral equations using left and right parts of the expression:
Fl = zeros(size(A));
Fr = zeros(size(A));
% A ~= 0, B ~= 0:
[IJ1] = find(A+B*m1~=0);
[IJ2] = find(A+B*m2~=0);
Fl(IJ2) = ff(IJ2)./B(IJ2).* ...
    ((exp(A(IJ2)*(x2-xstart)+B(IJ2)*(y22-ystart))-exp(A(IJ2)*(x1-xstart)+...
    B(IJ2)*(y21-ystart)))./(A(IJ2)+B(IJ2)*m2));
Fr(IJ1) = -ff(IJ1)./B(IJ1).* ...
    ((exp(A(IJ1)*(x2-xstart)+B(IJ1)*(y12-ystart))-exp(A(IJ1)*(x1-xstart)+...
    B(IJ1)*(y11-ystart)))./(A(IJ1)+B(IJ1)*m1));
[IJ1] = find(A+B*m1==0);
[IJ2] = find(A+B*m2==0);
Fl(IJ2) = ff(IJ2)./(B(IJ2).* ...
    exp(A(IJ2)*xstart+B(IJ2)*ystart)).*(x2-x1).*exp(B(IJ2)*c2);
Fr(IJ1) = ff(IJ1)./(B(IJ1).* ...
    exp(A(IJ1)*xstart+B(IJ1)*ystart)).*(x2-x1).*exp(B(IJ1)*c1);
% A = 0, B ~= 0:
if m1 == 0
    Fr(2:end,1) = -ff(2:end,1)./(B(2:end,1).^2.*exp(B(2:end,1)*ystart)).* ... % I THINK THIS IS BROKEN
        (x2-x1).*exp(B(2:end,1)*c1);
else
    Fr(2:end,1) = -ff(2:end,1)./(B(2:end,1).^2.*exp(B(2:end,1)*ystart)).* ...
        ((exp(B(2:end,1)*y12)-exp(B(2:end,1)*y11))/m1);
end
if m2 == 0
    Fl(2:end,1) = ff(2:end,1)./(B(2:end,1).^2.*exp(B(2:end,1)*ystart)).* ... % I THINK THIS IS BROKEN
        (x2-x1).*exp(B(2:end,1)*c2);
else
    Fl(2:end,1) = ff(2:end,1)./(B(2:end,1).^2.*exp(B(2:end,1)*ystart)).* ...
        ((exp(B(2:end,1)*y22)-exp(B(2:end,1)*y21))/m2);
end
% A ~= 0, B = 0:
Fl(1,2:end) = ff(1,2:end)./A(1,2:end).*(exp(A(1,2:end)*(x2-xstart)).*(y22-y12-(m2-m1)./A(1,2:end))-exp(A(1,2:end)*(x1-xstart)).*(y21-y11-(m2-m1)./A(1,2:end)));
Fr(1,2:end) = 0;
% A = 0, B = 0:
Fl(1,1) = ff(1,1)*((m2-m1)/2*(x2^2-x1^2)+(c2-c1)*(x2-x1));
Fr(1,1) = 0;

FF = Fl+Fr;
inttotal = real(sum(FF(:))/n^2);
toc;

% Estimate actual integral
midptsx = 0.25*(x(1:end-1,1:end-1)+x(1:end-1,2:end)+x(2:end,1:end-1)+x(2:end,2:end));
midptsy = 0.25*(y(1:end-1,1:end-1)+y(1:end-1,2:end)+y(2:end,1:end-1)+y(2:end,2:end));
for j = 1:length(midptsx)
    j
    for k = 1:length(midptsy)
        
        xx = midptsx(k,j);
        yy = midptsy(k,j);
        
        data(k,j) = 1/n^2*sum(sum((ff).*exp(2*pi*1i*Freqsx*(xx-xstart)+2*pi*1i*Freqsy*(yy-ystart))));
        
    end
end
data = real(data);
BB = polyhedronField(verticesA,magA,[midptsx(:),midptsy(:),repmat(z,size(midptsx(:)))]);
Bz = reshape(BB(:,3),size(midptsx));
dx = x(1:end-1,2:end)-x(1:end-1,1:end-1);
dy = y(2:end,:)-y(1:end-1,:);
area = 0.5*(dy(:,2:end)+dy(:,1:end-1)).*dx;

actual = sum(sum(Bz.*area))
estimate = real(inttotal)
pcerror = (estimate-actual)/actual*100

figure;
surf(midptsx,midptsy,Bz);

figure;
surf(midptsx,midptsy,data);

figure;
surf(midptsx,midptsy,data-Bz);

figure;
subplot(1,3,2)
surf(midptsx,midptsy,data)
title('FFT approximation using 256x256 grid')
subplot(1,3,1)
surf(midptsx,midptsy,Bz)
title('Field data')
subplot(1,3,3)
surf(midptsx,midptsy,data-Bz);
title('Error')



