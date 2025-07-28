clear; clc; close all;

% Dati iniziali
fluido = 'R290';
p_4 = 3192000;
p_5 = 1263158;
h_4 = 6.761019860295737e+05;
M_r = 0.13216;
A_leak = 4.6e-6;
A_su = 27.43e-6;
rv = 2.5;
V_s_exp = 100e-6;
T_loss = 0.47;

for eta_pol = 1.00:-0.05:0.20
    % Calcolo prestazioni espansore per ogni eta_pol
    h_su = h_4;
    s_su = py.CoolProp.CoolProp.PropsSI('S','P',p_4,'H',h_4,fluido);

    p_zero = p_4;
    fun_1 = @(p_su_1)abs(M_r - A_su*py.CoolProp.CoolProp.PropsSI('D','S',s_su,'P',p_su_1,fluido)*sqrt(2*(h_su - py.CoolProp.CoolProp.PropsSI('H','S',s_su,'P',p_su_1,fluido))));
    options = optimoptions('fmincon', 'Display','none');
    p_su_1 = fmincon(fun_1,p_zero,[],[],[],[],0,p_4,[],options);
    h_su_1 = h_su;
    s_su_1 = py.CoolProp.CoolProp.PropsSI('S','S',s_su,'P',p_su_1,fluido);
    v_su_1 = 1/(py.CoolProp.CoolProp.PropsSI('D','S',s_su,'P',p_su_1,fluido));

    h_su_2 = h_su_1;
    s_su_2 = s_su_1;
    v_su_2 = v_su_1;
    p_su_2 = p_su_1;

    cp_su_2 = py.CoolProp.CoolProp.PropsSI('C','H',h_su_2,'P',p_su_2,fluido);
    R_su_2_un = py.CoolProp.CoolProp.PropsSI('GAS_CONSTANT','H',h_su_2,'P',p_su_2,fluido);
    massa_molare = py.CoolProp.CoolProp.PropsSI('M','H',h_su_2,'P',p_su_2,fluido);
    R_su_2 = R_su_2_un/massa_molare;
    gamma_su_2 = cp_su_2/(cp_su_2 - R_su_2);

    p_ex_2 = p_5;
    p_crit_leak = p_su_2*((2/(gamma_su_2 + 1))^(gamma_su_2/(gamma_su_2 - 1)));
    p_thr_leak = max(p_ex_2,p_crit_leak);
    v_thr_leak = 1/(py.CoolProp.CoolProp.PropsSI('D','P',p_thr_leak,'S',s_su_2,fluido));
    h_thr_leak = py.CoolProp.CoolProp.PropsSI('H','P',p_thr_leak,'S',s_su_2,fluido);

    M_leak = A_leak/v_thr_leak*sqrt(2*(h_su_2 - h_thr_leak));
    M_in = M_r - M_leak;
    N = (60*v_su_2*M_in)/V_s_exp;

    m = gamma_su_2 / (gamma_su_2 - eta_pol * (gamma_su_2 - 1));
    p_ad = p_su_2*(1/rv)^m;
    v_ad = v_su_2*rv;
    d_ad = 1 / v_ad;
    h_ad = py.CoolProp.CoolProp.PropsSI('H','P',p_ad,'D',d_ad,fluido);

    W_in = M_in*((h_su_2 - h_ad) + v_ad*(p_ad - p_ex_2));
    W_loss = (2*pi*N)/60*T_loss;
    W_sh = W_in - W_loss;

    h_ex_2 = h_su_2 - W_in/M_in;
    h_ex_1 = (h_su_2*M_leak + h_ex_2*M_in)/M_r;

    h_5 = h_ex_1;
    T_5 = py.CoolProp.CoolProp.PropsSI('T','P',p_5,'H',h_5,fluido);

    % Visualizzazione i risultati per ogni eta_pol
    fprintf(['eta_pol = %.2f | W_sh = %.2f W | T_5 = %.2f K | h_5 = %.2f J/kg\n'], eta_pol, W_sh, T_5, h_5);
end
