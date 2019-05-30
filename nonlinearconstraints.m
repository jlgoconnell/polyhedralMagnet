
function [c,ceq] = nonlinearconstraints(h,theta,M,V,d)

c = [(h-sqrt(V/h*tand(theta)-h^2/3))/tand(theta),h^2/3-V/h*tand(theta)];
ceq = [];

end