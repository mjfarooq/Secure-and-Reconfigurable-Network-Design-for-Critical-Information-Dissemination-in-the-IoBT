function [c,ceq] = const(x,alpha,T1,T2,Tc)
c(1) =  (1/(alpha*(1 - Tc))) - x(1)^2 * x(2) *pi*x(3)^2  -  x(2)*pi*x(4)^2;
c(2) = (1/(alpha*(1 - T1))) - x(1)^2 * x(2) * pi*x(3)^2;
c(3) = (1/(alpha*(1 - T2))) - x(2) * pi * x(4)^2;

ceq = [];