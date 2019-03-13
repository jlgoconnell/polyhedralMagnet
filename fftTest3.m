
% Learning how to use FFT
%
% James O'Connell 11th March 2019

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

n = 64;
x = linspace(-0.03,0.03,n);
y = linspace(-0.03,0.03,n);
[X,Y] = meshgrid(x,y);
B = polyhedronField(verticesA,magA,[X(:),Y(:),repmat(0.023,size(X(:)))]);
Bz = reshape(B(:,3),size(X));
f = Bz;
% f = cos(2*pi*1*X)+2*cos(2*pi*2*Y+2*pi*3*X)+3*cos(2*pi*3*Y);

figure;
plot3(X(:),Y(:),f(:),'r.');
hold on;

ff = fft2(f);
% ff = fftshift(ff);

freqsx = 1/(x(2)-x(1))*(0:(n-1))/n;
freqsy = 1/(y(2)-y(1))*(0:(n-1))/n;

freqsx(length(freqsx)/2:end) = freqsx(length(freqsx)/2:end)'-1/(x(2)-x(1));
freqsy(length(freqsy)/2:end) = freqsy(length(freqsy)/2:end)'-1/(x(2)-x(1));


% freqsx = (0:n-1);
% freqsy = freqsx;

% ff = ff/n^2;

% ff = ff(1:n/2+1,1:n/2+1);
% ff(2:end-1,2:end-1) = 2*ff(2:end-1,2:end-1);
% freqsx = freqsx(1:n/2+1);
% freqsy = freqsy(1:n/2+1);

n2 = n*4;
x2 = linspace(x(1),x(end),n2);
y2 = linspace(y(1),y(end),n2);
[X2,Y2] = meshgrid(x2,y2);
% f = cos(2*pi*1*X)+2*cos(2*pi*2*Y+2*pi*3*X);
data = zeros(length(x2),length(y2));

for j = 1:length(x2)
    j
    for k = 1:length(y2)
        
        [Freqsx,Freqsy] = meshgrid(freqsx,freqsy);
        
        xx = x2(j);
        yy = y2(k);
        
        data(k,j) = 1/n^2*sum(sum((ff).*exp(2*pi*1i*Freqsx*(xx-x(1))+2*pi*1i*Freqsy*(yy-y(1)))));
        
    end
end
data = real(data);

mydata = ifft2(ff);

% f2 = cos(2*pi*1*X2)+2*cos(2*pi*2*Y2+2*pi*3*X2)+3*cos(2*pi*3*Y2);
surf(X2,Y2,data);
% plot3(X2(:),Y2(:),f2(:),'ro');
% surf(X,Y,mydata);
grid on;

% error = f-data;
% figure;
% surf(X2,Y2,error);