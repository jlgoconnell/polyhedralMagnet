
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

myd = 0.001:0.0005:0.06;
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
F = -F;

yyaxis left
plot(myd*1000,X(:,1))
hold on
grid on
yyaxis right
plot(myd*1000,X(:,2))