clear
close all
clc

% INPUT E PARAMETRI

fluido = 'R290';     % fluido di lavoro
p_1 = 1200000;
p_2 = 3360000;
p_3 = p_2;
p_4 = 3192000;
p_5 = 1263158;
p_6 = p_5;
%delta_p_r = 0.03;    % perdita di pressione relativa
M_r = 0.13216;      % portata effettivamente erogata dalla pompa [kg/s]
Q_ev_zero = 50000;   % potenza termica fornita all'evaporatore iniziale [W]
Q_cond_zero = 46700; % potenza termica estratta al condensatore iniziale [W]
eta_pp = 0.75;       % rendimento isoentropico della pompa
eta_rig = 0.7;       % efficienza del rigeneratore ((T_3 - T_2)/(T_5 - T_2))
eta_pol = 1;          % rendimento politropico
A_leak = 4.6e-6;     % area della sezione di perdita [m^2]
A_su = 27.43e-6;     % area della sezione di aspirazione [m^2]
rv = 2.5;           % rapporto volumetrico della macchina
V_s_exp = 100e-6;  % cilindrata dell'espansore [m^3]
T_loss = 0.47;       % coppia meccanica persa [Nm]
scelta_rigeneratore = 0; % 0=rigeneratore assente; 1=rigeneratore presente

% CALCOLO RISULTATI DELL'IMPIANTO

% calcolo pressioni e set point
% p_6 = p_5*(1 - (delta_p_r*scelta_rigeneratore)); % pressione al punto 6 [Pa]
% p_1 = p_6*(1 - delta_p_r); % pressione al punto 1 [Pa]
% p_3 = p_4*(1 + delta_p_r); % pressione al punto 3 [Pa]
% p_2 = p_3*(1 + (delta_p_r*scelta_rigeneratore)); % pressione al punto 2 [Pa]

T_1_cli = py.CoolProp.CoolProp.PropsSI('T','P',p_1,'Q',0,fluido); % temperatura al punto 1 sulla cli [K]
T_1_sp = T_1_cli -5; % temperatura al punto 1 di set point [K]


T_4_cls = py.CoolProp.CoolProp.PropsSI('T','P',p_4,'Q',1,fluido); % temperatura al punto 4 sulla cls [K]
T_4_sp = T_4_cls +8; % temperatura al punto 4 di set point [K]

% ricerca delle potenze ottimali
fun = @(X) fun_cont_pol(M_r,p_1,p_2,p_3,p_4,p_5,p_6,X,T_4_sp,T_1_sp,eta_rig,eta_pp,eta_pol,A_leak,A_su,T_loss,V_s_exp,rv,fluido,scelta_rigeneratore);

tic

options = optimoptions('fmincon','Display','none');
X = fmincon(fun,[Q_ev_zero Q_cond_zero],[],[],[],[],[0 0],[51000 48000],[],options); % ricerca potenze termiche agli scambiatori di calore ottimali [W]
Q_ev_ott = X(1); % potenza termica fornita all'evaporatore ottimale [W]
Q_cond_ott = X(2); % potenza termica estratta al condensatore ottimale [W]

% calcolo degli output
[W_sh,W_pp,N,T_1,T_2,T_3,T_4,T_5,T_6,M_leak] = fun_orc_pol(M_r,p_1,p_2,p_3,p_4,p_5,p_6,X,eta_rig,eta_pp,eta_pol,A_leak,A_su,T_loss,V_s_exp,rv,fluido,scelta_rigeneratore);

tempo = toc/60; % tempo di calcolo [min]

W_ut = W_sh - W_pp; % Potenza utile [W]
rendimento = W_ut/Q_ev_ott; % Rendimento

% calcolo dei titoli
x_1 = py.CoolProp.CoolProp.PropsSI('Q','P',p_1,'T',T_1,fluido); % titolo al punto 1
x_2 = py.CoolProp.CoolProp.PropsSI('Q','P',p_2,'T',T_2,fluido); % titolo al punto 2
x_3 = py.CoolProp.CoolProp.PropsSI('Q','P',p_3,'T',T_3,fluido); % titolo al punto 3
x_4 = py.CoolProp.CoolProp.PropsSI('Q','P',p_4,'T',T_4,fluido); % titolo al punto 4
x_5 = py.CoolProp.CoolProp.PropsSI('Q','P',p_5,'T',T_5,fluido); % titolo al punto 5
x_6 = py.CoolProp.CoolProp.PropsSI('Q','P',p_6,'T',T_6,fluido); % titolo al punto 6

% calcolo delle entropie
s_1 = py.CoolProp.CoolProp.PropsSI('S','P',p_1,'T',T_1,fluido); % entropia specifica al punto 1 [J/kgK]
s_2 = py.CoolProp.CoolProp.PropsSI('S','P',p_2,'T',T_2,fluido); % entropia specifica al punto 2 [J/kgK]
s_3 = py.CoolProp.CoolProp.PropsSI('S','P',p_3,'T',T_3,fluido); % entropia specifica al punto 3 [J/kgK]
s_4 = py.CoolProp.CoolProp.PropsSI('S','P',p_4,'T',T_4,fluido); % entropia specifica al punto 4 [J/kgK]
s_5 = py.CoolProp.CoolProp.PropsSI('S','P',p_5,'T',T_5,fluido); % entropia specifica al punto 5 [J/kgK]
s_6 = py.CoolProp.CoolProp.PropsSI('S','P',p_6,'T',T_6,fluido); % entropia specifica al punto 6 [J/kgK]

% calcolo delle entalpie
h_1 = py.CoolProp.CoolProp.PropsSI('H','P',p_1,'T',T_1,fluido); % entropia specifica al punto 1 [J/kgK]
h_2 = py.CoolProp.CoolProp.PropsSI('H','P',p_2,'T',T_2,fluido); % entropia specifica al punto 2 [J/kgK]
h_3 = py.CoolProp.CoolProp.PropsSI('H','P',p_3,'T',T_3,fluido); % entropia specifica al punto 3 [J/kgK]
h_4 = py.CoolProp.CoolProp.PropsSI('H','P',p_4,'T',T_4,fluido); % entropia specifica al punto 4 [J/kgK]
h_5 = py.CoolProp.CoolProp.PropsSI('H','P',p_5,'T',T_5,fluido); % entropia specifica al punto 5 [J/kgK]
h_6 = py.CoolProp.CoolProp.PropsSI('H','P',p_6,'T',T_6,fluido); % entropia specifica al punto 6 [J/kgK]

% COSTRUZIONE DIAGRAMMA T-s
% punti sulle curve limite
x_cli = 0; % titolo sulla cli
x_cls = 1; % titolo sulla cls

T_crit = py.CoolProp.CoolProp.PropsSI('Tcrit',fluido); % temperatura critica del fluido [K]
T_range = linspace(288.15,T_crit); % range di temperatura fino a T_crit [K]

for i = 1:length(T_range)
    s_cli(i) = py.CoolProp.CoolProp.PropsSI('S','T',T_range(i),'Q',x_cli,fluido); % entropia sulla cli [J/kgK]
    s_cls(i) = py.CoolProp.CoolProp.PropsSI('S','T',T_range(i),'Q',x_cls,fluido); % entropia sulla cls [J/kgK]
end

% plot curve

figure; hold on
plot(s_cli,T_range,'c')
plot(s_cls,T_range,'c')
xlabel('s [J/kgK]')
ylabel('T [K]')
title('Curva T-s')
grid on

% plot punti
plot(s_1,T_1,'b.','MarkerSize',20) % punto 1
plot(s_2,T_2,'g.','MarkerSize',20) % punto 2
plot(s_3,T_3,'r.','MarkerSize',20) % punto 3
plot(s_4,T_4,'m.','MarkerSize',20) % punto 4
plot(s_5,T_5,'y.','MarkerSize',20) % punto 5
plot(s_6,T_6,'k.','MarkerSize',20) % punto 6

legend('cli','cls','punto 1','punto 2','punto 3','punto 4','punto 5','punto 6') % legenda