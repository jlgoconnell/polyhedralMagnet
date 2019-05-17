
% Gradient descent method for optimising frusta
%
% James O'Connell 16th May 2019

clear;
clc;
close all;

V = 1e-4; % Frustum volume in m^3 (equal to 100cm^3)

M = [0,0,1]; % Magnetisation vector in Teslas

% Initial guess:
a(1) = 65; % Wall andgle in degrees
a(2) = 0.05; % Height in m

delta = [0.001,0.0001]; % Step size for calculating derivatives

gamma = 10*delta; % Step size for iteration

d = 0.01; % Separation distance in m

i = 1;

while i <= 100
    
    F = zeros(2);
    a(1) = deg2rad(a(1));
    
    b = (sqrt(V/a(2)*tan(a(1))^2-a(2)^2/3)+a(2))/tan(a(1));
    bu = b-2*a(2)/tan(a(1));

    if b > 0 && b == real(b) && bu > 0 && bu == real(bu) % Make sure the frustum can physically exist
        verticesA = [-b/2,-b/2,0;-b/2,b/2,0;b/2,-b/2,0;b/2,b/2,0; ...
            -bu/2,-bu/2,-a(2);-bu/2,bu/2,-a(2);bu/2,-bu/2,-a(2);bu/2,bu/2,-a(2)];
        verticesB = -verticesA + repmat([0,0,d],size(verticesA,1),1);

        [myF,~,~] = polyhedronForce(verticesB,verticesA,M,M,12,mean(verticesB));
        F(2,1) = myF(3);
    else
        F(2,1) = NaN;
    end
    F(2,2) = F(2,1);
    
    atemp = [a(1)+delta(1),a(2)];
    b = (sqrt(V/atemp(2)*tan(atemp(1))^2-atemp(2)^2/3)+atemp(2))/tan(atemp(1));
    bu = b-2*atemp(2)/tan(atemp(1));

    if b > 0 && b == real(b) && bu > 0 && bu == real(bu) % Make sure the frustum can physically exist
        verticesA = [-b/2,-b/2,0;-b/2,b/2,0;b/2,-b/2,0;b/2,b/2,0; ...
            -bu/2,-bu/2,-atemp(2);-bu/2,bu/2,-atemp(2);bu/2,-bu/2,-atemp(2);bu/2,bu/2,-atemp(2)];
        verticesB = -verticesA + repmat([0,0,d],size(verticesA,1),1);

        [myF,~,~] = polyhedronForce(verticesB,verticesA,M,M,12,mean(verticesB));
        F(1,1) = myF(3);
    else
        F(1,1) = NaN;
    end
    
    atemp = [a(1),a(2)+delta(2)];
    b = (sqrt(V/atemp(2)*tan(atemp(1))^2-atemp(2)^2/3)+atemp(2))/tan(atemp(1));
    bu = b-2*atemp(2)/tan(atemp(1));

    if b > 0 && b == real(b) && bu > 0 && bu == real(bu) % Make sure the frustum can physically exist
        verticesA = [-b/2,-b/2,0;-b/2,b/2,0;b/2,-b/2,0;b/2,b/2,0; ...
            -bu/2,-bu/2,-atemp(2);-bu/2,bu/2,-atemp(2);bu/2,-bu/2,-atemp(2);bu/2,bu/2,-atemp(2)];
        verticesB = -verticesA + repmat([0,0,d],size(verticesA,1),1);

        [myF,~,~] = polyhedronForce(verticesB,verticesA,M,M,12,mean(verticesB));
        F(1,2) = myF(3);
    else
        F(1,2) = NaN;
    end
    
    gradF = (F(1,:)-F(2,:))./delta;
    
    F
    
    a(1) = rad2deg(a(1));
    
    a = a + [0.0001,0.00001].*gradF;
    
    i = i + 1;
end

a(2) = a(2) * 1000;

a