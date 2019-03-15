
% Learning how to use FFT
%
% James O'Connell 10th March 2019

clear;
% close all;
clc;

n = 64;
<<<<<<< Updated upstream
x = linspace(0,2,n)-2;
y = linspace(0,2,n)-2;
=======
x = linspace(0,1,n);
y = linspace(0,1,n);
>>>>>>> Stashed changes
[X,Y] = meshgrid(x,y);
f = cos(2*pi*1*X)+2*cos(2*pi*2*Y+2*pi*3*X)+3*cos(2*pi*3*Y);

figure;
plot3(X(:),Y(:),f(:),'r.');
hold on;

ff = fft2(f);

xx = (x-x(1))/(x(end)-x(1));
yy = (y-y(1))/(y(end)-y(1));

dx = x(2)-x(1);
fs = 1/dx;
df = fs/length(x);

freqsx = 1/(x(2)-x(1))*(0:(n-1))/n;
freqsy = 1/(y(2)-y(1))*(0:(n-1))/n;

% freqsx = -fs/2:df:fs/2-df;
% freqsy = freqsx;

% freqsx = fftshift(freqsx);
% freqsy = fftshift(freqsy);

ff = ff/n^2;
% ff = fftshift(ff);

% ff = 2*ff(1:n/2+1,1:n/2+1);
% % ff(2:end-1,2:end-1) = 2*ff(2:end-1,2:end-1);
% freqsx = freqsx(1:n/2+1);
% freqsy = freqsy(1:n/2+1);

n2 = n*1;
x2 = linspace(x(1),x(end),n2);
y2 = linspace(y(1),y(end),n2);
[X2,Y2] = meshgrid(x2,y2);
% f = cos(2*pi*1*X)+2*cos(2*pi*2*Y+2*pi*3*X);
data = zeros(length(x2),length(y2));

a = 1

for j = 1:length(x2)
    
    for k = 1:length(y2)
        
        [Freqsx,Freqsy] = meshgrid(freqsx,freqsy);
        
        xx = x2(j);
        yy = y2(k);
        
        data(k,j) = sum(sum((ff).*exp(2*pi*1i*Freqsx*(xx-x(1))).*exp(2*pi*1i*Freqsy*(yy-y(1)))));
        
    end
end
data = real(data);

% mydata = ifft2(ff*n^2);

f2 = cos(2*pi*1*X2)+2*cos(2*pi*2*Y2+2*pi*3*X2)+3*cos(2*pi*3*Y2);
surf(X2,Y2,data);
plot3(X2(:),Y2(:),f2(:),'r.');
grid on;

% f = cos(2*pi*1*X2)+2*cos(2*pi*2*Y2+2*pi*3*X2)+3*cos(2*pi*3*Y2);
error = f2-data;
figure;
surf(X2,Y2,error);