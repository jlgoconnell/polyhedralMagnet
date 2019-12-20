
% Function to calculate the magnetic field due to ferrous materials with
% non-unity permeability

% James O'Connell 17th December 2019

function [B,F,T,t] = calcMagField(magnets,obspt,meshDensity,torquepts,varargin)

tic;

if nargin < 4
    torquepts = [];
    for i = 1:length(magnets)
        torquepts(i,:) = mean(magnets{i}.vertices);
    end
end

mu0 = 4*pi*10^-7;

% Set up basic magnet mesh
numMags = length(magnets);
whichmag = [];
tris = [0,0,0];
pts = [];
mur = zeros(numMags,1);
M = zeros(numMags,3);
for i = 1:numMags
    fac{i} = minConvexHull(magnets{i}.vertices);
    mur(i) = magnets{i}.mur;
    M(i,:) = magnets{i}.magnetisation;
    tritemp = triangulateFaces(fac{i});
    tris = [tris;tritemp+max(tris(:))];
    pts = [pts;magnets{i}.vertices];
    whichmag = [whichmag;i*ones(size(tritemp,1),1)];
end
tris = tris(2:end,:);

% Subdivide mesh and calculate properties of each element
if meshDensity > 1
    [pts,tris] = subdivideMesh(pts,tris,meshDensity);
    whichmag = repelem(whichmag,meshDensity^2);
end
mur = mur(whichmag);

% Calculating field with unity charge density
centres = meshFaceCentroids(pts,tris);
fieldpts = [centres;obspt];
for i = 1:length(tris)
    Btemp = triangleField(pts(tris(i,:),:),1,fieldpts)/mu0;
    
    Bx(:,i) = Btemp(:,1);
    By(:,i) = Btemp(:,2);
    Bz(:,i) = Btemp(:,3);
end

% Calculate the A matrix and solve for sigma
Bxonmag = Bx(1:length(tris),:);
Byonmag = By(1:length(tris),:);
Bzonmag = Bz(1:length(tris),:);
norms = meshFaceNormals(pts,tris);
M = M(whichmag,:);
sigma0 = dot(M',norms')';
areas = meshFaceAreas(pts,tris);
nx = repmat(norms(:,1)',length(norms),1)';
ny = repmat(norms(:,2)',length(norms),1)';
nz = repmat(norms(:,3)',length(norms),1)';
beta = Bxonmag.*nx + Byonmag.*ny + Bzonmag.*nz;
r = diag((mur-1)./mur);
A = eye(size(beta)) - r*beta;
sigma = A^(-1)*sigma0;

% Calculate total force and torque on each magnet
Bextx = mu0*(Bxonmag.*(1-eye(size(Bxonmag))))*sigma;
Bexty = mu0*(Byonmag.*(1-eye(size(Byonmag))))*sigma;
Bextz = mu0*(Bzonmag.*(1-eye(size(Bzonmag))))*sigma;
f = diag(areas.*sigma)*[Bextx,Bexty,Bextz];
r = torquepts(whichmag)-centres;
t = diag(areas.*sigma)*cross(r,[Bextx,Bexty,Bextz]);
F = zeros(length(magnets),3);
T = zeros(length(magnets),3);
for i = 1:length(magnets)
    F(i,:) = sum(f(whichmag==i,:));
    T(i,:) = sum(t(whichmag==i,:));
end

% Test code:
myBxonmag = Bxonmag;
myByonmag = Byonmag;
myBzonmag = Bzonmag;
[a,b] = meshgrid(whichmag,whichmag');
mask = a~=b;
myBxonmag = myBxonmag.*mask;
myByonmag = myByonmag.*mask;
myBzonmag = myBzonmag.*mask;
myBextx = mu0*myBxonmag*sigma;
myBexty = mu0*myByonmag*sigma;
myBextz = mu0*myBzonmag*sigma;
myf = diag(areas.*sigma)*[myBextx,myBexty,myBextz];
myr = torquepts(whichmag)-centres;
myt = diag(areas.*sigma)*cross(myr,[myBextx,myBexty,myBextz]);
myF = zeros(length(magnets),3);
for i = 1:length(magnets)
    myF(i,:) = sum(myf(whichmag==i,:));
    myT(i,:) = sum(myt(whichmag==i,:));
end

% Calculate field at the desired points
if size(obspt,1) > 0
    B = mu0*[Bx(length(tris)+1:end,:)*sigma,By(length(tris)+1:end,:)*sigma,Bz(length(tris)+1:end,:)*sigma];
else
    B = [];
end

figure;
hold on;
for i = 1:length(tris)
    c = (sigma(i)-min(sigma))/(max(sigma)-min(sigma));
    a = drawMesh(pts,tris(i,:),[c,0.5-abs(0.5-c),1-c]);
%     a.EdgeColor = 'none';
%     a.FaceAlpha = 0.2;
end
% drawMesh(pts,tris,'FaceAlpha',0.1);
xlabel('x');
ylabel('y');
zlabel('z');
hold on;
quiver3(centres(:,1),centres(:,2),centres(:,3),myf(:,1),myf(:,2),myf(:,3));

t = toc;

end



