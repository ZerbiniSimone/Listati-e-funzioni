clear
close all
clc

portata = [0.13216 0.11564 0.09912 0.08260 0.06608]; % matrice portata parziale [kg/s]

T_5 = [318.5 319.4 320.5 322.1 324.5]; % matrice temperatura uscita espansore [K]

W_sh = [3775.2 3094.5 2413.9 1733.3 1052.7]; % matrice potenza meccanica [W]

N_exp = [740.7 607.1 473.6 340.1 206.5]; % matrice velocità di rotazione [rpm]

Q_ev = [49800 43575 37350 31125 24900]; % matrice potenza termica fornita all'evaporatore [W]

Q_cond = [46707 41080 35452 29824 24197]; % matrice potenza termica estratta al condensatore [W]

W_ut = [2991.8 2409.1 1826.4 1243.7 661]; % matrice potenza utile [W]

eta = [0.0601 0.0553 0.0489 0.04 0.0265]; % matrice rendimento

figure;
set(gcf, 'Name', 'Espansore', 'NumberTitle', 'off');

subplot(1,3,1);
plot(portata, T_5, '-.*c');
title('Temperatura di scarico espansore');
xlabel('M_r [kg/s]');
ylabel('T_5 [K]');
grid on

subplot(1,3,2);
plot(portata, W_sh, '-.*k');
title('Potenza meccanica erogata');
xlabel('M_r [kg/s]');
ylabel('W_s_h [W]');
grid on

subplot(1,3,3)
plot(portata, N_exp, '-.*m');
title('Velocità di rotazione');
xlabel('M_r [kg/s]');
ylabel('N_e_x_p [rpm]');
grid on

figure;
set(gcf, 'Name', 'Potenze termiche', 'NumberTitle', 'off');
subplot(1,2,1);

plot(portata, Q_ev, '-.*g');
title('Potenza termica fornita');
xlabel('M_r [kg/s]');
ylabel('Q_e_v [W]');
grid on

subplot(1,2,2);
plot(portata, Q_cond, '-.*b');
title('Potenza termica estratta');
xlabel('M_r [kg/s]');
ylabel('Q_c_o_n_d [W]');
grid on

figure;
set(gcf, 'Name', 'Prestazioni', 'NumberTitle', 'off');

subplot(1,2,1);
plot(portata, W_ut, '-.*r');
title('Potenza utile');
xlabel('M_r [kg/s]');
ylabel('W_u_t [W]');
grid on

subplot(1,2,2);
plot(portata, eta, '-.*g');
title('Rendimento');
xlabel('M_r [kg/s]');
ylabel('η_i_m_p');
grid on


