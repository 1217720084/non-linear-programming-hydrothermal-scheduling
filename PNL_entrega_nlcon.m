function [c,ceq] = EX_PNL_nlcon(x)
    a0 = 7.352460e+02;
    a1 = 3.496580e-03;
    a2 = -1.974370e-07;
    a3 = 6.917050e-12;
    a4 = -9.773650e-17;
    b0 = 6.71633e+02;
    b1 = 1.01738e-03;
    b2 = -1.79972e-07;
    b3 = 2.51328e-11;
    b4 = 0.00000e+00;
    k = 0.008633;   
    T=length(x)/4; %numero de variaveis =3
    d=1312;
    %restrições não lineares igualdades
    ceq=zeros(1,T);
    for t=1:1:T
        %ceq(t)=g(t)+k*q(t)*(a0+a1*x(t)+a2*x(t)^2+a3*x(t)^3+a4*x(t)^4-b0 -b1*q(t)-b2*q(t)^2-b3*q(t)^3-b4*q(t)^4 )-d(t) = 0
        ceq(t)=x(t+2*T)+k*x(t+T)*(a0+a1*x(t)+a2*(x(t)^2)+a3*(x(t)^3)+a4*(x(t)^4)-b0-b1*x(t+T)-b2*(x(t+T)^2)-b3*(x(t+T)^3)-b4*(x(t+T)^4))-d;
    end 
    %restrições não lineares desigualdades
    c = []; %<=0
end