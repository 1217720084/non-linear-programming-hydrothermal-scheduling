clc
clear
%%
%definição de variáveis
x(1) = 0; % x(t=1)
x(2) = 0; % x(t=2)
x(3) = 0; % x(t=3)
x(4) = 0; % x(t=4)
x(5) = 0; % x(t=5)
x(6) = 0; % x(t=6)
x(7) = 0; % x(t=7)
x(8) = 0; % x(t=8)
x(9) = 0; % x(t=9)
x(10) = 0; % x(t=10)
x(11) = 0; % x(t=11)
x(12) = 0; % x(t=12)

x(13) = 0; % g(t=1)
x(14) = 0; % g(t=2)
x(15) = 0; % g(t=3)
x(16) = 0; % g(t=4)
x(17) = 0; % g(t=5)
x(18) = 0; % g(t=6)
x(19) = 0; % g(t=7)
x(20) = 0; % g(t=8)
x(21) = 0; % g(t=9)
x(22) = 0; % g(t=10)
x(23) = 0; % g(t=11)
x(24) = 0; % g(t=12)

x(25) = 0; % q(t=1)
x(26) = 0; % q(t=2)
x(27) = 0; % q(t=3)
x(28) = 0; % q(t=4)
x(29) = 0; % q(t=5)
x(30) = 0; % q(t=6)
x(31) = 0; % q(t=7)
x(32) = 0; % q(t=8)
x(33) = 0; % q(t=9)
x(34) = 0; % q(t=10)
x(35) = 0; % q(t=11)
x(36) = 0; % q(t=12)

y = [740 606 500 412 433 509 721 1234 1770 1658 1480 1014];
%d = ones(1,12)*1312;
%%
%definição de limites
xmin = 5733;
xmax = 22950;
qmin = 200;
qmax = 1692;

lb = zeros(36,1);
lb([1:12],1) = xmin; %lower bound x
%o problema não define g, mas é >0
lb([25:36],1) = qmin; %lower bound q

ub = inf(36,1);
ub([1:12],1) = xmax; %upper bound x
%o problema não define g, mas é ilimitado
ub([25:36],1) = qmax; %upper bound q

%%
%restrições lineares = 0:
%x(t+1)-x(t)+beta*q(t)=beta*y(t)
beta = 2.628;
x_inicial = 22950;
x_final = 22950;

Aeq = zeros(14,36);
%x(1) = 22950 %t=1
Aeq(1,1) = 1;
%x(2) - x(1) + beta*x(25)= beta*y(1); %t=2
Aeq(2,2) = 1;
Aeq(2,1) = -1;
Aeq(2,25) = beta;
%x(3) - x(2) + beta*x(26)= beta*y(2); %t=3
Aeq(3,3) = 1;
Aeq(3,2) = -1;
Aeq(3,26) = beta;
%x(4) - x(3) + beta*x(27)= beta*y(3); %t=4
Aeq(4,4) = 1;
Aeq(4,3) = -1;
Aeq(4,27) = beta;
%x(5) - x(4) + beta*x(28)= beta*y(4); %t=5
Aeq(5,5) = 1;
Aeq(5,4) = -1;
Aeq(5,28) = beta;
%x(6) - x(5) + beta*x(29)= beta*y(5); %t=6
Aeq(6,6) = 1;
Aeq(6,5) = -1;
Aeq(6,29) = beta;
%x(7) - x(6) + beta*x(30)= beta*y(6); %t=7
Aeq(7,7) = 1;
Aeq(7,6) = -1;
Aeq(7,30) = beta;
%x(8) - x(7) + beta*x(31)= beta*y(7); %t=8
Aeq(8,8) = 1;
Aeq(8,7) = -1;
Aeq(8,31) = beta;
%x(9) - x(8) + beta*x(32)= beta*y(8); %t=9
Aeq(9,9) = 1;
Aeq(9,8) = -1;
Aeq(9,32) = beta;
%x(10) - x(9) + beta*x(33)= beta*y(9); %t=10
Aeq(10,10) = 1;
Aeq(10,9) = -1;
Aeq(10,33) = beta;
%x(11) - x(10) + beta*x(34)= beta*y(10); %t=11
Aeq(11,11) = 1;
Aeq(11,10) = -1;
Aeq(11,34) = beta;
%x(12) - x(11) + beta*x(35) = beta*y(11); %t=12
Aeq(12,12) = 1;
Aeq(12,11) = -1;
Aeq(12,35) = beta;
%x(12) = 22950 %t=12
Aeq(13,12) = 1;

beq = zeros(14,1);
beq(1,1) = x_inicial;
beq(2:12,1) = beta*y(1:11);
beq(13,1) = x_final;
%restrições lineares <= 0
A = [];
b = [];
%%
%restrições não lineares
nonlincon = @EX_PNL_nlcon;
%min soma de x([14:25])
J = @(x) 0.02*((x(13)^2)+(x(14)^2)+(x(15)^2)+(x(16)^2)+(x(17)^2)+(x(18)^2)+(x(19)^2)+(x(20)^2)+(x(21)^2)+(x(22)^2)+(x(23)^2)+(x(24)^2));
options = optimset('MaxFunEvals', 999999, 'MaxIter', 99999); %necessário aumentar o numero de iterações para convergir
[X,fval,ef,output,lambda] = fmincon(J,x,A,b,Aeq,beq,lb,ub,nonlincon,options);
fval
t = 1:1:12;
A = [transpose(t) transpose(X(1:12)) transpose(X(13:24)) transpose(X(25:36))];
T = array2table(A,'VariableNames',{'t', 'x','g','q'})
output
%%
custo = @(g) 0.02*g^2;
Reservatorio=X([1:12]);
Reservatorio(1)=22950;
Vazao=X(25:36);
geracao=X(13:24);
Custo_final=zeros(1,12);
for i=12:-1:1
    if(i==12)
        Custo_final(i)=custo(geracao(i));
    else
        Custo_final(i)=custo(geracao(i))+Custo_final(i+1);
    end
end

for i=1:1:size(Reservatorio,2)
    plot(i,Reservatorio(i),'g.','MarkerSize',20);
    hold on;
    text(i,X(i)+110,strcat('c*=',num2str(fix(Custo_final(i)))));
    text(i,X(i),strcat('q*=',num2str(fix(Vazao(i)))));
    if(i<12)
        plot([i i+1], [Reservatorio(i) Reservatorio(i+1)], 'g');
    end;
end
axis([0 13 5733-1000 22950+1000]);
xlabel('Estágio t');
ylabel('Volume Reservatório [hm³]');