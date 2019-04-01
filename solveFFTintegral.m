
% Function to calculate the FFT integrals solved by hand
% 
% James O'Connell 1st April 2019

function F = solveFFTintegral(ff,A,B,xstart,ystart,trappts,n)

x1 = trappts(1,1);
x2 = trappts(3,1);
y11 = trappts(1,2);
y12 = trappts(3,2);
y21 = trappts(2,2);
y22 = trappts(4,2);

m1 = (y12-y11)/(x2-x1);
m2 = (y22-y21)/(x2-x1);
c1 = y11-m1*x1;
c2 = y21-m2*x1;

% Solve integral equations using left and right parts of the expression:
Fl = zeros(size(A));
Fr = zeros(size(A));
% A ~= 0, B ~= 0:
[IJ1] = find(A+B*m1~=0);
[IJ2] = find(A+B*m2~=0);
Fl(IJ2) = ff(IJ2)./B(IJ2).* ...
    ((exp(A(IJ2)*(x2-xstart)+B(IJ2)*(y22-ystart))-exp(A(IJ2)*(x1-xstart)+...
    B(IJ2)*(y21-ystart)))./(A(IJ2)+B(IJ2)*m2));
Fr(IJ1) = -ff(IJ1)./B(IJ1).* ...
    ((exp(A(IJ1)*(x2-xstart)+B(IJ1)*(y12-ystart))-exp(A(IJ1)*(x1-xstart)+...
    B(IJ1)*(y11-ystart)))./(A(IJ1)+B(IJ1)*m1));
[IJ1] = find(A+B*m1==0);
[IJ2] = find(A+B*m2==0);
Fl(IJ2) = ff(IJ2)./(B(IJ2).* ...
    exp(A(IJ2)*xstart+B(IJ2)*ystart)).*(x2-x1).*exp(B(IJ2)*c2);
Fr(IJ1) = ff(IJ1)./(B(IJ1).* ...
    exp(A(IJ1)*xstart+B(IJ1)*ystart)).*(x2-x1).*exp(B(IJ1)*c1);
% A = 0, B ~= 0:
if m1 == 0
    Fr(2:end,1) = -ff(2:end,1)./B(2:end,1).*(x2-x1).*exp(B(2:end,1)*(c1-ystart));
else
    Fr(2:end,1) = -ff(2:end,1)./(B(2:end,1).^2.*exp(B(2:end,1)*ystart)).* ...
        ((exp(B(2:end,1)*y12)-exp(B(2:end,1)*y11))/m1);
end
if m2 == 0
    Fl(2:end,1) = ff(2:end,1)./B(2:end,1).*(x2-x1).*exp(B(2:end,1)*(c2-ystart));
else
    Fl(2:end,1) = ff(2:end,1)./(B(2:end,1).^2.*exp(B(2:end,1)*ystart)).* ...
        ((exp(B(2:end,1)*y22)-exp(B(2:end,1)*y21))/m2);
end
% A ~= 0, B = 0:
Fl(1,2:end) = ff(1,2:end)./A(1,2:end).*(exp(A(1,2:end)*(x2-xstart)).*(y22-y12-(m2-m1)./A(1,2:end))-exp(A(1,2:end)*(x1-xstart)).*(y21-y11-(m2-m1)./A(1,2:end)));
Fr(1,2:end) = 0;
% A = 0, B = 0:
Fl(1,1) = ff(1,1)*((m2-m1)/2*(x2^2-x1^2)+(c2-c1)*(x2-x1));
Fr(1,1) = 0;

FF = Fl+Fr;
F = real(sum(FF(:))/n^2);





end