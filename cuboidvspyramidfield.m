
% Compare the magnetic field strength of a cuboid and a pyramid magnet

% James O'Connell 11th September 2019

clear;
close all;
clc;

V = 8e-6;
h = 0.02;
mag = [0,0,1.3];

num = 1000;
x = linspace(-0.02,0.02,num);
y = 0;
z = linspace(0.00,0.004,num+1)';
z = z(2:end);
[X,Y,Z] = meshgrid(x,y,z);
obspt = [X(:),Y(:),Z(:)];

% Theta = pi/180*linspace(atand(2/sqrt(3)),90,4);
Theta = pi/180*linspace(60,90,2);

for i = 1:length(Theta)
    theta = Theta(i);
    b = (-6*h/tan(theta)+sqrt(36*V/h-12*h^2/(tan(theta))^2))/6;
    bb = b + 2*h/tan(theta);

    vertices = [b/2,b/2,0;b/2,-b/2,0;-b/2,b/2,0;-b/2,-b/2,0;...
        bb/2,bb/2,-h;bb/2,-bb/2,-h;-bb/2,bb/2,-h;-bb/2,-bb/2,-h];
%     S = alphaShape(vertices,inf); figure; plot(S);

    BB = polyhedronField(vertices,mag,obspt);
    Bx(:,i) = BB(:,1);
    Bz(:,i) = BB(:,3);
end

xx = reshape(X,sqrt(numel(X(:))),sqrt(numel(X(:))));
zz = reshape(Z,sqrt(numel(Z(:))),sqrt(numel(Z(:))));

Btot = sqrt(Bx.^2+Bz.^2);
figure;
titles = {'Frustum magnet','Cube magnet'};
for i = 1:length(Theta)
    Btott = reshape(Btot(:,i),sqrt(numel(Btot(:,i))),sqrt(numel(Btot(:,i))));
    
    subplot(1,length(Theta),i);
    [C,h] = contourf(1000*xx,1000*zz,Btott);
    set(h,'LineColor','none');
%     title(['\theta = ',num2str(180/pi*Theta(i)),' degrees'], 'Interpreter','tex');
    title(titles{i});
    c = colorbar('southoutside');
    c.Label.String = 'Field strength (T)';
    caxis([0,0.8]);
    xlabel('x (mm)');
    ylabel('z (mm)');
    ylim([0,1000*max(z)]);
    
    
    
end

% figure;
% plot(1000*z,B);
% legend(num2str(180/pi*Theta'))
% grid on;
% xlabel('z (mm)');
% ylabel('Field strength (T)');