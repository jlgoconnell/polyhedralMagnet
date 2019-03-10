
% Learning how to use FFT
%
% James O'Connell 10th March 2019

clear;
close all;
clc;

n = 16;
x = linspace(0,1,n);
y = linspace(0,1,n);
[X,Y] = meshgrid(x,y);
f = cos(2*pi*1*X)+2*cos(2*pi*2*Y+2*pi*3*X)+3*cos(2*pi*3*Y);

figure;
plot3(X(:),Y(:),f(:),'r.');
hold on;

ff = fft2(f);
ff = ff/n^2;

freqsx = 1/(x(2)-x(1))*(0:(n-1))/n;
freqsy = 1/(y(2)-y(1))*(0:(n-1))/n;


n2 = n*16;
x2 = linspace(x(1),x(end),n2);
y2 = linspace(y(1),y(end),n2);
[X2,Y2] = meshgrid(x2,y2);
% f = cos(2*pi*1*X)+2*cos(2*pi*2*Y+2*pi*3*X);
data = zeros(length(x2),length(y2));

a = 1

for j = 1:length(x2)
    j
    for k = 1:length(y2)
        
        [Freqsx,Freqsy] = meshgrid(freqsx,freqsy);
        
        xx = x2(j);
        yy = y2(k);
        
        data(k,j) = sum(sum((ff).*exp(2*pi*1i*Freqsx*xx).*exp(2*pi*1i*Freqsy*yy)));
        
    end
end
data = real(data);

f2 = cos(2*pi*1*X2)+2*cos(2*pi*2*Y2+2*pi*3*X2)+3*cos(2*pi*3*Y2);
surf(X2,Y2,data);
plot3(X2(:),Y2(:),f2(:),'ro');

% error = f-data;
% figure;
% surf(X2,Y2,error);