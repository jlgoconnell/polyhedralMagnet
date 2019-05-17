

function [F] = calcFrustaForce(h,theta,M,V,d)

% X = load('frustumvars.mat');
% V = X.V;
% d = X.d;
% M = X.M;
% 
% h = x(1);
% theta = x(2);

theta = deg2rad(theta);

b = (sqrt(V/h*tan(theta)^2-h^2/3)+h)/tan(theta);
bu = b-2*h/tan(theta);

if b > 0 && bu > 0 && real(b) == b && real(bu) == bu

    verticesA = [-b/2,-b/2,0;-b/2,b/2,0;b/2,-b/2,0;b/2,b/2,0; ...
                    -bu/2,-bu/2,-h;-bu/2,bu/2,-h;bu/2,-bu/2,-h;bu/2,bu/2,-h];
    verticesB = -verticesA + repmat([0,0,d],size(verticesA,1),1);

    [myF,~,~] = polyhedronForce(verticesB,verticesA,M,M,16,mean(verticesB));
    F = -myF(3);

else
    F = NaN;
end


end