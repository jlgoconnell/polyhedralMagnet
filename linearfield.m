
% Attempt to find the most linear magnetic field above a vertically
% magnetised frustum

% James O'Connell 11th September 2019

clear;
close all;
clc;

V = 8e-6;
h = 0.02;
mag = [0,0,1];

num = 1000;
x = 0*ones(num,1);
y = 0*ones(num,1);
z = linspace(0,0.01,num+1)';
z = z(2:end);
obspt = [x,y,z];

Theta = pi/180*linspace(60,120,7);

for i = 1:length(Theta)
    theta = Theta(i);
    b = (-6*h/tan(theta)+sqrt(36*V/h-12*h^2/(tan(theta))^2))/6;
    bb = b + 2*h/tan(theta);

    vertices = [b/2,b/2,0;b/2,-b/2,0;-b/2,b/2,0;-b/2,-b/2,0;...
        bb/2,bb/2,-h;bb/2,-bb/2,-h;-bb/2,bb/2,-h;-bb/2,-bb/2,-h];
%     S = alphaShape(vertices,inf); figure; plot(S);

    BB = polyhedronField(vertices,mag,obspt);
    B(:,i) = BB(:,3);
end

figure;
plot(1000*z,B);
legend(num2str(180/pi*Theta'))
grid on;
xlabel('z (mm)');
ylabel('Field strength (T)');

% D = diff(B);
% figure;
% plot(D)
% legend(num2str(180/pi*Theta'))