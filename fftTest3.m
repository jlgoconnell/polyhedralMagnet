
% Learning how to use FFT
%
% James O'Connell 10th March 2019

clear;
close all;
clc;

n = 64;
x = linspace(0,1,n);
y = linspace(0,1,n);
[X,Y] = meshgrid(x,y);
f = cos(2*pi*1*X)+2*cos(2*pi*2*Y+2*pi*3*X);

figure;
plot3(X(:),Y(:),f(:),'ro');
hold on;

ff = fft2(f);
ff = ff/n^2;

freqsx = 1/(x(2)-x(1))*(0:(n-1))/n;
freqsy = 1/(y(2)-y(1))*(0:(n-1))/n;

data = zeros(length(x),length(y));

n2 = n;
x2 = linspace(0,1,n2);
y2 = linspace(0,1,n2);
[X2,Y2] = meshgrid(x2,y2);
% f = cos(2*pi*1*X)+2*cos(2*pi*2*Y+2*pi*3*X);

a = 1

[Freqsx,Freqsy] = meshgrid(freqsx,freqsy);
Freqsx = reshape(Freqsx,[1,1,size(Freqsx)]);
Freqsy = reshape(Freqsy,[1,1,size(Freqsy)]);
% X = repmat(X,[1,1,size(Freqsx)]);
% Y = repmat(Y,[1,1,size(Freqsy)]);


data = sum(sum(ff.*(exp(2*pi*1i*Freqsx.*X).*exp(2*pi*1i*Freqsy.*Y)),3),4);
data = real(data);
% 
% 
% for j = 1:length(x)
%     for k = 1:length(y)
%         
%         [Freqsx,Freqsy] = meshgrid(freqsx,freqsy);
%         
%         xx = x(j);
%         yy = y(k);
%         
%         data(k,j) = sum(sum((ff).*exp(2*pi*1i*Freqsx*xx).*exp(2*pi*1i*Freqsy*yy)));
%         
%     end
% end
% data = real(data);
% 
surf(X,Y,data);

% error = f-data;
% figure;
% surf(X,Y,error);