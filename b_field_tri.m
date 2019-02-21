
% This function calculates the B-field around a triangle with vertices at
% (0,0), (x,y1), and (x,y2). Note that the line opposite the origin-vertex
% must be parallel to the y-axis.

function B = b_field_tri(x,z_obs,z_tri,y1,y2)

c = z_obs-z_tri;
if abs(c) < eps
    c = sign(c)*1e-6;%eps;
end

Da1 = sqrt(x^2+y1^2);
Da2 = sqrt(x^2+y2^2);
Dac1 = sqrt(x^2+c^2+y1^2);
Dac2 = sqrt(x^2+c^2+y2^2);

Bx = 1/2*(log((Dac2+y2)*(Dac1-y1)/((Dac2-y2)*(Dac1+y1))) + y1/Da1*log((Dac1+Da1)/(Dac1-Da1)) - y2/Da2*log((Dac2+Da2)/(Dac2-Da2)));
By = x/2*(1/Da2*log((Dac2+Da2)/(Dac2-Da2)) - 1/Da1*log((Dac1+Da1)/(Dac1-Da1)));
Bz = c/abs(c)*atan2(x*y2-x*y1,y1*y2+x^2) + atan2(x*c*y1*Dac2-x*c*y2*Dac1,c^2*y1*y2+x^2*Dac1*Dac2);

B = [Bx,By,Bz];

end