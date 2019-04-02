
% Function to estimate the force and torque between two polyhedral
% permanent magnets, both with constant uniform magnetisation. Let one
% magnet be magnet A and the other be magnet B. This function calculates
% the force and torque on magnet B.
%
% Inputs:
% verticesA: The vertices of one magnet in an (n x 3) matrix.
% verticesB: The vertices of the magnet we are calculating the force and
% torque on in an (n x 3) matrix.
% magA: The magnetisation vector of magnet A in Teslas.
% magB: The magnetisation vector of magnet B in Teslas.
% nFFT: The exponent of 2 which will be used as the number of gridpoints
% for the FFT. Larger integers will produce more accurate results but take
% longer to compute. Recommended range: 3 to 7.
%
% Outputs:
% F: A vector representing the x, y, and z forces on magnet B.
%
% James O'Connell 1st April 2019

function F = polyhedronForceFFT(verticesA,verticesB,magA,magB,nFFT)

tic;

F = [0,0,0]';

Fac = minConvexHull(verticesB);
[Ver,~] = surfToMesh(verticesB(:,1),verticesB(:,2),verticesB(:,3));
norms = meshFaceNormals(Ver,Fac);
MdotN = 1/(pi*4e-7)*dot(repmat(magB,size(norms,1),1)',norms')';
Fac = Fac(abs(MdotN)>eps);
norms = meshFaceNormals(Ver,Fac);

for i = 1:length(Fac)
    
    Fface = [0,0,0]';
    
    facepts = Ver(Fac{i},:);
    
    % Rotate the face to be parallel to the XY plane
    n = norms(i,:)/norm(norms(i,:));
    thetay = atan2(n(1),-n(3));
    thetax = atan2(n(2),-sqrt(n(1)^2+n(3)^2));
    Ry = [cos(thetay),0,sin(thetay);0,1,0;-sin(thetay),0,cos(thetay)];
    Rx = [1,0,0;0,cos(thetax),-sin(thetax);0,sin(thetax),cos(thetax)];
    Rn = Rx*Ry;
    R = Rn^-1;
    facepts = facepts*R;
    z = facepts(1,3);
    Apts = verticesA*R;
    Amag = magA*R;
    
    % Set up FFT grid and perform FFT
    x1 = min(facepts(:,1));
    x2 = max(facepts(:,1));
    y1 = min(facepts(:,2));
    y2 = max(facepts(:,2));
    n = 2^nFFT;
    % To avoid Gibb's phenomenon, take an area slightly larger
    deltad = 0.01;
    xstart = x1-deltad*(x2-x1);
    xend = x2+deltad*(x2-x1);
    ystart = y1-deltad*(y2-y1);
    yend = y2+deltad*(y2-y1);
    x = linspace(xstart,xend,n);
    y = linspace(ystart,yend,n);
    [X,Y] = meshgrid(x,y);
    B = polyhedronField(Apts,Amag,[X(:),Y(:),repmat(z,size(X(:)))]);
    size(X(:))
    Bx = reshape(B(:,1),size(X));
    ffx = fft2(Bx);
    By = reshape(B(:,2),size(X));
    ffy = fft2(By);
    Bz = reshape(B(:,3),size(X));
    ffz = fft2(Bz);
    
    % Calculate spatial frequencies
    freqsx = 1/(x(2)-x(1))*(0:(n-1))/n;
    freqsy = 1/(y(2)-y(1))*(0:(n-1))/n;
    freqsx = fftshift(freqsx);
    freqsy = fftshift(freqsy);
    freqsx = freqsx - freqsx(1);
    freqsy = freqsy - freqsy(1);
    [Freqsx,Freqsy] = meshgrid(freqsx,freqsy);
    % Coefficients of x and y in the exponent of e:
    A = 2*pi*1i*Freqsx;
    B = 2*pi*1i*Freqsy;
    
    % Decompose the facet into trapezia:
    trapVertices = trapDecomp(facepts(:,1:2));
    for j = 1:2:length(trapVertices)-3
        
        mytrap = trapVertices(j:j+3,:);
%         plot(mytrap(:,1),mytrap(:,2),'ro');
        
        Ftrapx = solveFFTintegral(ffx,A,B,xstart,ystart,mytrap,n);
        Ftrapy = solveFFTintegral(ffy,A,B,xstart,ystart,mytrap,n);
        Ftrapz = solveFFTintegral(ffz,A,B,xstart,ystart,mytrap,n);
        
        Fface = Fface + MdotN(i)*[Ftrapx;Ftrapy;Ftrapz];
        
        
    end
    Fface = Fface'*Rn;
    F = F + Fface';
    
end

toc;






end