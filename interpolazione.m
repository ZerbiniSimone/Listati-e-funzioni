clear
close all
clc

Q_ev = 45000; % Calore fornito richiesto [W]

% Dati
Q_ev_vettore = [49800 43575 37350 31125 24900]; % Calore fornito all'evaporatore [W]
M_r_vettore = [0.013216 0.011564 0.009912 0.00826 0.006608]; % Portata massica del fluido [kg/s]

% Calcolo portata per interpolazione
M_r = interp1(Q_ev_vettore, M_r_vettore, Q_ev, 'linear'); % Portata massica relativa al calore fornito immesso [kg/s]

% Dati
W_sh_vettore = [3775.2 3094.5 2413.9 1733.3 1052.7]; % Potenza meccanica [W]
Q_cond_vettore = [46707 41080 35452 29824 24197]; % Calore scaricato al condensatore
eta_vettore = [0.0601 0.0553 0.0489 0.04 0.0265]; % Efficienza dell'impianto

% Calcolo potenza meccanica, calore estratto ed efficienza per interpolazione
W_sh = interp1(M_r_vettore, W_sh_vettore, M_r, 'linear'); % Potenza meccanica relativa alla portata calcolata [W]
Q_cond = interp1(M_r_vettore, Q_cond_vettore, M_r, 'linear'); % Calore estratto relativo alla portata calcolata [W]
eta = interp1(M_r_vettore, eta_vettore, M_r, 'linear'); % Efficienza relativa alla portata calcolata

