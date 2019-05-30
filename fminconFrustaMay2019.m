
clear;
clc;
close all;

V = 8e-6;
M = [0,0,1];

x0 = [V^(1/3),75];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [0,0];
ub = [];

myd = 0.0005:0.0005:0.1;
F = zeros(length(myd),0);

for i = 1:length(myd)
    d = myd(i);
    
    fun = @(x) calcFrustaForce(x(1),x(2),M,V,d);

    [x,F(i)] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,@(x) nonlinearconstraints(x(1),x(2),M,V,d));
    x(1) = x(1)*1000;

    X(i,:) = x;
    x0 = [x(1)/1000,x(2)];
    
    i/length(myd)*100
end
F = -F';

yyaxis left
plot(myd*1000,X(:,1))
hold on
grid on
yyaxis right
plot(myd*1000,X(:,2))

pcscale = 0.4;
l = (V/pcscale)^(1/3);
h = (pcscale^2*V)^(1/3);
magnet_fixed = magnetdefine('type','cuboid','dim',[l l h],'magn',1,'magdir','z');
magnet_float = magnetdefine('type','cuboid','dim',[l l h],'magn',1,'magdir','z');
N = length(myd);
offset = repmat([0; 0; h],[1 N]);
displ_range = offset+[0;0;1]*myd;
Fcuboid = magnetforces(magnet_fixed,magnet_float,displ_range);
Fcuboid = -Fcuboid(3,:)';

pcincrease = (F-Fcuboid)./Fcuboid*100;




