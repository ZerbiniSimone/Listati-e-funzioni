clear; clc; close all;

% Dati iniziali
fluido = 'R290';     % fluido di lavoro
p_4 = 3192000;       % Pa
p_5 = 1263158;       % Pa
h_4 = 6.761019860295737e+05; % J/kg
M_r = 0.13216;       % kg/s
A_su = 27.43e-6;     % m^2
A_leak = 4.6e-6;     % m^2
rv = 2.5;            % rapporto volumetrico
V_s_exp = 100e-6;    % m^3
T_loss = 0.47;       % Nm
eta_pol = 1;

% Calcoli preliminari espansore
h_su = h_4;
s_su = py.CoolProp.CoolProp.PropsSI('S','P',p_4,'H',h_4,fluido);
fun_1 = @(p_su_1)abs(M_r - A_su*py.CoolProp.CoolProp.PropsSI('D','S',s_su,'P',p_su_1,fluido)*sqrt(2*(h_su - py.CoolProp.CoolProp.PropsSI('H','S',s_su,'P',p_su_1,fluido))));
options = optimoptions('fmincon', 'Display','none');
p_su_1 = fmincon(fun_1,p_4,[],[],[],[],0,p_4,[],options);
v_su_1 = 1/(py.CoolProp.CoolProp.PropsSI('D','S',s_su,'P',p_su_1,fluido));
s_su_1 = py.CoolProp.CoolProp.PropsSI('S','S',s_su,'P',p_su_1,fluido);
h_su_1 = h_su;

p_su_2 = p_su_1;
h_su_2 = h_su_1;
v_su_2 = v_su_1;
s_su_2 = s_su_1;

% Calcolo portata persa (M_leak) e interna (M_in)
cp_su_2_zero = py.CoolProp.CoolProp.PropsSI('C','H',h_su_2,'P',p_su_2,fluido); 
R_su_2_un = py.CoolProp.CoolProp.PropsSI('GAS_CONSTANT','H',h_su_2,'P',p_su_2,fluido); 
massa_molare_zero = py.CoolProp.CoolProp.PropsSI('M','H',h_su_2,'P',p_su_2,fluido);
R_su_2_init = R_su_2_un/massa_molare_zero;
gamma_su_2_init = cp_su_2_zero/(cp_su_2_zero - R_su_2_init);
p_ex_2 = p_5;
p_crit_leak = p_su_2*((2/(gamma_su_2_init + 1))^(gamma_su_2_init/(gamma_su_2_init - 1)));
p_thr_leak = max(p_ex_2,p_crit_leak);
v_thr_leak = 1/(py.CoolProp.CoolProp.PropsSI('D','P',p_thr_leak,'S',s_su_2,fluido));
h_thr_leak = py.CoolProp.CoolProp.PropsSI('H','P',p_thr_leak,'S',s_su_2,fluido);
M_leak = A_leak/v_thr_leak*sqrt(2*(h_su_2 - h_thr_leak));
M_in = M_r - M_leak;
N = (60*v_su_2*M_in)/V_s_exp;
W_loss = (2*pi*N)/60*T_loss;

% Inizio Analisi Multi-Stadio 
n_stadi = 1:20;
h_ris = zeros(size(n_stadi));
W_ris = zeros(size(n_stadi));
T_ris = zeros(size(n_stadi));

% Ciclo esterno: varia il numero di stadi totali 'n'
for idx = 1:length(n_stadi)
    n = n_stadi(idx);
    
    % Inizializzazione delle condizioni di partenza per l'espansione a 'n' stadi
    p_stadio = p_su_2;
    h_stadio = h_su_2;
    v_stadio = v_su_2;
    
    % Calcolo del rapporto volumetrico per ogni stadio
    rv_stadio = rv^(1/n);
    
    % Ciclo interno: esegue l'espansione per ogni stadio da 1 a 'n'
    for i = 1:n
        % Aggiornamento delle proprietà del fluido all'inizio dello stadio corrente
        cp_curr = py.CoolProp.CoolProp.PropsSI('C','H',h_stadio,'P',p_stadio,fluido);
        massa_molare = py.CoolProp.CoolProp.PropsSI('M','H',h_stadio,'P',p_stadio,fluido);
        R_curr = R_su_2_un / massa_molare;
        gamma_curr = cp_curr / (cp_curr - R_curr);
        
        % Calcolo esponente della politropica per lo stadio corrente
        m = gamma_curr / (gamma_curr - eta_pol * (gamma_curr - 1));
        
        % Calcolo dello stato a fine stadio
        p_next = p_stadio * (1/rv_stadio)^m;
        v_next = v_stadio * rv_stadio;
        d_next = 1 / v_next;
        h_next = py.CoolProp.CoolProp.PropsSI('H','P',p_next,'D',d_next,fluido);
        
        % Aggiornamento dello stato per lo stadio successivo
        p_stadio = p_next;
        v_stadio = v_next;
        h_stadio = h_next;
    end
    
    % Alla fine del ciclo interno, 'h_sc' è l'entalpia finale 'h_ad'
    p_ad = p_stadio;
    v_ad = v_stadio;
    h_ad = h_stadio;
    
    % Calcolo delle prestazioni complessive per 'n' stadi
    W_in = M_in*((h_su_2 - h_ad) + v_ad*(p_ad - p_ex_2));
    W_sh = W_in - W_loss;
    
    h_ex_2 = h_su_2 - W_in/M_in;
    h_ex_1 = (h_su_2*M_leak + h_ex_2*M_in)/M_r;
    
    h_5 = h_ex_1;
    T_5 = py.CoolProp.CoolProp.PropsSI('T','P',p_5,'H',h_5,fluido);
    
    % Risultati per il grafico
    h_ris(idx) = h_5;
    W_ris(idx) = W_sh;
    T_ris(idx) = T_5;
end

% Creazione Grafici
figure('Name', 'Analisi Espansione Multi-Stadio');

% Grafico Entalpia
subplot(3,1,1);
plot(n_stadi, h_ris, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4);
title('Andamento di h_5 in funzione del numero di stadi');
xlabel('n');
ylabel('h_5 [J/kg]');
grid on;

% Grafico Lavoro
subplot(3,1,2);
plot(n_stadi, W_ris, 'g-s', 'LineWidth', 1.5, 'MarkerSize', 4);
title('Andamento di W_{sh} all''albero in funzione del numero di stadi');
xlabel('n');
ylabel('W_{sh} [W]');
grid on;

% Grafico Temperatura, invariata perché nella campana
subplot(3,1,3);
plot(n_stadi, T_ris, 'b-d', 'LineWidth', 1.5, 'MarkerSize', 4);
title('Andamento di T_5 in funzione del numero di stadi');
xlabel('n');
ylabel('T_5 [K]');
grid on;

sgtitle('Analisi delle Prestazioni al Variare degli Stadi di Espansione', 'FontWeight', 'bold');