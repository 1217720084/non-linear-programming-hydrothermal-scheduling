%%
%===============DEFINIÇÃO DE DADOS==========================
%===========================================================
clc
clear
%polinomio cota x volume
a0 = 7.352460e+02;
a1 = 3.496580e-03;
a2 = -1.974370e-07;
a3 = 6.917050e-12;
a4 = -9.773650e-17;
%polinomio de jusante
b0 = 6.71633e+02;
b1 = 1.01738e-03;
b2 = -1.79972e-07;
b3 = 2.51328e-11;
b4 = 0.00000e+00;
k = 0.008633;  
beta=2.6284; %conversão m3/s -> hm3 (365/12)*24*60*60/1E6
qmin=200; %mínima histórica deck NEWAVE é 102m3/s, mas a principio não ha limite
qmax=1692; %vazão dos dois cj de geradores 6*211+2*213 m3/s
Xmin=5737.5; %5737.5 é múltiplo de 1% do reservatório, próximo ao 5733 dado deck NEWAVE
Xmax=22950;
Y_DATA = csvread('VAZOES_FURNAS_CSV.txt'); %HISTÓRICO DECK 11/19 
Y_DATA=Y_DATA(:,2:13);
Y_MLT=mean(Y_DATA(1:2019-1931,:));

Y_AUX=[];
for index=1:1:length(Y_DATA)
    Y_AUX=[Y_AUX Y_DATA(index,:)]; 
end

Y_HISTORICO=Y_AUX(1:(2019-1931)*12); %01/31 até 12/18

% y=Y_AUX(1:(2019-1931)*12); %01/31 até 12/18
% y=Y_AUX(1:(1940-1931)*12); %01/31 até 12/18

% y=[2740 2606 2500 2412 2433 2509 2721 1234 1770 1658 1480 1014];


% y=[740 606 500 412 433 509 721 1234 1770 1658 1480 1014];
% y=[y y y y];
y=[Y_MLT(6:end) Y_MLT Y_MLT Y_MLT];

%%
%===============DEFINIÇÃO DE VARIÁVEIS======================
%===========================================================
T=length(y); 
X_otimo=zeros(1,T);
Q_otimo=zeros(1,T);
G_otimo=zeros(1,T);
Vertimento_otimo=zeros(1,T);
x=[X_otimo Q_otimo G_otimo Vertimento_otimo];
%%
%===============DEFINIÇÃO DE LIMITES========================
%===========================================================
lb = [Xmin*ones(1,length(X_otimo)) qmin*ones(1,length(Q_otimo)) zeros(1,length(G_otimo)) zeros(1,length(Vertimento_otimo))];
ub = [Xmax*ones(1,length(X_otimo)) qmax*ones(1,length(Q_otimo)) 1312*ones(1,length(G_otimo)) 9999*ones(1,length(Vertimento_otimo))];
%%
%===========DEFINIÇÃO DE RESTRIÇÕES LINEARES================
%===========================================================
%restrições lineares igualdades
Aeq=zeros(T+1,4*T); %T restrições de transição de estado e 1 de condição inicial
beq=zeros(T+1,1);
%x(t+1)-x(t)+beta*q(t)+beta*v(t)=beta*y(t) transição de estado
for t=1:1:T
    Aeq(t,t+1)=1; %x(t+1)
    Aeq(t,t)=-1; %x(t)
    Aeq(t,t+T)=beta; %q(t)
    Aeq(t,t+3*T)=beta; %vertimento(t)
    beq(t,1)=beta*y(t); %y(t)
end
%X(t=1)=Xmax estado inicial
Aeq(T+1,1)=1;
beq(T+1,1)=Xmax/2;
% X(t=T)=Xmax estado final
% Aeq(T+1,T)=1;
% beq(T+1)=Xmax;
%restrições lineares desigualdades
A = [];
b = [];
%%
%=======DEFINIÇÃO DE RESTRIÇÕES NÃO LINEARES================
%===========================================================
nonlincon = @PNL_entrega_nlcon;
%%
%=============DEFINIÇÃO DA FUNÇÃO OBJETIVO==================
%===========================================================
J=@(x) 0.02*sum(x(2*T+1:3*T).^2);
% J = @(x) 0.02*((x(25)^2)+(x(26)^2)+(x(27)^2)+(x(28)^2)+(x(29)^2)+(x(30)^2)+(x(31)^2)+(x(32)^2)+(x(33)^2)+(x(34)^2)+(x(35)^2)+(x(36)^2));
%%
%=================DADOS ENTRADA FMINCON=====================
%===========================================================
tic
tstart=tic;
options = optimset('MaxFunEvals', 9999999, 'MaxIter', 999, 'Algorithm','interior-point'); %necessário aumentar o numero de iterações para convergir
[X,fval,ef,output,lambda] = fmincon(J,x,A,b,Aeq,beq,lb,ub,nonlincon,options);
fval;
output;
tempo_de_exec_PNL=toc(tstart)
%%
%=================DESMEMBRAR VARIAVEIS======================
%===========================================================
X_otimo=X(1:T);
Q_otimo=X(T+1:2*T);
G_otimo=X(2*T+1:3*T);
Vertimento_otimo=X(3*T+1:4*T);
%%
%====================CÁLCULO DE CUSTO=======================
%===========================================================

custo = @(g) 0.02*g^2;
for p=length(X_otimo):-1:1
    if(p==length(X_otimo))
        Custo_otimo(p)=custo(G_otimo(p));
    else
        Custo_otimo(p)=Custo_otimo(p+1)+custo(G_otimo(p));
    end
end
%%
%====LEITURA DE DADOS HYDROLAB PARA PLOT COMPARATIVO========
%===========================================================
% X_DATA = csvread('VOL_PDD_HYDROLAB.txt');
% X_otimo=[X_otimo;X_DATA.'];
% 
% Q_DATA = csvread('Q_PDD_HYDROLAB.txt');
% Q_otimo=[Q_otimo;Q_DATA.'];
% csvwrite('RG_VOL_PNL.txt',X_otimo(1,:).');
% csvwrite('RG_TURB_PNL.txt',Q_otimo(1,:).');
% csvwrite('RG_VERT_PNL.txt',Vertimento_otimo(1,:).');



%%
%===============SAÍDA DE DADOS PARA GRÁFICO=================
%===========================================================

x_axis=1:T;
% x_axis=datetime(1930,12,1)+calmonths(1:T);

subplot(3,1,1);
[ax_volume nivel_reservatorio afluencia ] = plotyy(x_axis,X_otimo,x_axis,y);
set(nivel_reservatorio,'linestyle','-','marker','>');
set(ax_volume(1),'ylim',[5000 24000]);
ax_volume(1).YTick = [5000:(24000-5000)/5:24000];
title(strcat('Volume do reservatório e afluência por estágio-',num2str(T),' meses'))
ylabel(ax_volume(1), 'Volume do Reservatório [hm³]');
xlabel('Estágio t');% 
ylabel(ax_volume(2), 'Afluência [m³/s]');
set(afluencia,'linestyle','--');
legend('PDD','PDD Hydrolab','Afluência');
text(0+0.2,X_otimo(1,12)+250,strcat('Cota mínima'));
text(0+0.2,X_otimo(1,1)+250,strcat('Cota máxima'));
hold on
plot(x_axis,[ones(1,T)*Xmin;ones(1,T)*Xmax],'Color',[0 0 0],'LineStyle',':')
hold off

subplot(3,1,2);
[ax_turb turbinagem afluencia ] = plotyy(x_axis,Q_otimo,x_axis,y);
set(turbinagem,'linestyle','-','marker','>');
title(strcat('Turbinagem e afluência por estágio-',num2str(T),' meses'));
ylabel(ax_turb(1), 'Turbinagem [m³/s]');
xlabel('Estágio t');% 
ylabel(ax_turb(2), 'Afluência [m³/s]');
set(afluencia,'linestyle','--');
legend('PDD','PDD Hydrolab','Afluência');
hold on
plot(x_axis,[ones(1,T)*qmin;ones(1,T)*qmax],'Color',[0 0 0],'LineStyle',':')
hold off

subplot(3,1,3);
[ax_vert vert afluencia ] = plotyy(x_axis,Vertimento_otimo,x_axis,y);
set(vert,'linestyle','-','marker','>');
legend('PDD','PDD Hydrolab','Afluência');
title(strcat('Vertimento e afluência por estágio-',num2str(T),' meses'));
ylabel(ax_vert(1), 'Vertimento [m³/s]');
xlabel('Estágio t');% 
ylabel(ax_vert(2), 'Afluência [m³/s]');
set(afluencia,'linestyle','--');


linkaxes([ax_volume, ax_turb], 'x')


figure
plot(x_axis,Custo_otimo);
legend('PDD Custo Total','PDD Hydrolab Custo Total');
title(strcat('Custo total por estágio-',num2str(T),' meses'));
ylabel('Custo [R$]');
xlabel('Estágio t');
csvwrite(strcat(num2str(T),'_VOL_PNL.txt'),X_otimo.');
csvwrite(strcat(num2str(T),'_TURB_PNL.txt'),Q_otimo.');
csvwrite(strcat(num2str(T),'_VERT_PNL.txt'),Vertimento_otimo.');
csvwrite(strcat(num2str(T),'_CUSTO_PNL.txt'),Custo_otimo.');
csvwrite(strcat(num2str(T),'_PRODT_PNL.txt'),G_otimo.');
csvwrite(strcat(num2str(T),'_TEMPO_EXEC_PNL.txt'),tempo_de_exec_PNL);


