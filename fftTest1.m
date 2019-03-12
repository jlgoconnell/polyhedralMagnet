
% Learning how to use FFT
%
% James O'Connell 9th March 2019

clear;
close all;
clc;

n = 16;
x = linspace(0,0.01,n);

f = cos(2*pi*10*x)+2*sin(2*pi*22*x);

ff = fft(f);
estf = ifft(ff);
ff = ff/n;
% ff = ff(1:n/2+1);
% ff(2:end-1) = 2*ff(2:end-1);

phase = angle(ff);

a = (n-1)/(max(x)-min(x));

k = 1:n;

% ff = fftshift(ff);

% freqs = a*x.*(k-1)/n;
freqs = 1/(x(2)-x(1))*(0:(n-1))/n;

freqs = freqs(1:length(ff));

% plot(freqs,abs(ff));
% xlim([0,10]);
% ffshift = [ff(1:n/2),fliplr(ff(2:(n/2+1)))];
% fff = (ff');
% freqss = freqs';
% phase = phase';
% fff = repmat(fff,1,length(x));
% freqss = repmat(freqss,1,length(x));
% phase = repmat(phase,1,length(x));
% xx = repmat(x,size(fff,1),1);
% data = fff.*(cos(2*pi*freqss.*xx-phase)+1i*sin(2*pi*freqss.*xx-phase));
% data = sum(data,1);

for j = 1:length(x)
    
    xx = x(j);
    
%     ff = ff.*(abs(ff)>0.1);
%     phase = phase.*(abs(ff)>0.1);
    
    data(j) = sum(abs(ff).*cos(2*pi*freqs*xx+phase));
    data(j) = sum((ff).*exp(2*pi*1i*freqs*xx));
    
end
data = real(data);



figure;
plot(x,f,'.--');
hold on
plot(x,data,'-');
% plot(x,estf,'--');
legend('Real','FFT');