
% A halbach array using trapezial magnets to increase the field strength

% James O'Connell 11th September 2019

clear;
close all;
clc;

db = linspace(-0.005,0.005,11);

h = 0.01;

nummags = 16;
x = linspace(0.01*nummags*0.25,0.01*nummags*0.75,100);
z = linspace(0.004,0.006,11);
% z = z(2:end);
y = 0;
[X,Y,Z] = meshgrid(x,y,z);
% X = X(:);
% Y = Y(:);
% Z = Z(:);

Bx = zeros(numel(X),length(db));
By = zeros(numel(X),length(db));
Bz = zeros(numel(X),length(db));

for i = 1:length(db)
    d = db(i);
    
    b = 0.005+d;
    t = 0.005-d;
   
%     figure;
%     hold on;
    for j = 1:nummags
        x1 = 0.01*(j-1);
        x2 = 0.01*j;
        vertices = [x1,t,0;x1,-t,0;x1,b,-h;x1,-b,-h;x2,t,0;x2,-t,0;x2,b,-h;x2,-b,-h];
%         S = alphaShape(vertices,inf);
%         plot(S);
        mag = round([cos(pi*j/2),0,sin(pi*j/2)]);
        B = polyhedronField(vertices,mag,[X(:),Y(:),Z(:)]);
        Bx(:,i) = Bx(:,i) + B(:,1);
        Bz(:,i) = Bz(:,i) + B(:,3);
    end
    
%     figure;
%     quiver3(X(:),Y(:),Z(:),Bx(:,i),By(:,i),Bz(:,i));
%     view([0,1,0]);
end

Btot = sqrt(Bx.^2+Bz.^2);

plot(db,mean(Btot));

% X = reshape(X,size(X,2),size(X,3));
% Z = reshape(Z,size(Z,2),size(Z,3));
% figure;
% for i = 1:size(Btot,2)
%     Btott = reshape(Btot(:,i),size(X));
%     
%     subplot(size(Btot,2),1,i);
%     contourf(X,Z,Btott);
%     colorbar;
% end








