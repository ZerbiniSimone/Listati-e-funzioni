function[W_sh, W_pp, N, T_1, T_2, T_3, T_4, T_5, T_6, M_leak] = fun_orc_pol(M_r, p_1, p_2, p_3, p_4, p_5, p_6, X, eta_rig, eta_pp, eta_pol, A_leak, A_su, T_loss, V_s_exp, rv, fluido, scelta_rigeneratore)
% FUNZIONE CHE RACCHIUDE IL MODELLO DELL'IMPIANTO ORC

% definizione potenze termiche
Q_ev = X(1);   % calore fornito all'evaporatore [W]
Q_cond = X(2); % calore estratto al condensatore [W]
% fine definizione potenze termiche

% definizione dati iniziali per ciclo iterativo
h_1 = py.CoolProp.CoolProp.PropsSI('H','P',p_1,'Q',0,fluido); % entalpia specifica punto 1 sulla cli [J/kg]
T_1_cli = py.CoolProp.CoolProp.PropsSI('T','P',p_1,'Q',0,fluido); % temperatura punto 1 sulla cli [K]
T_4_cls = py.CoolProp.CoolProp.PropsSI('T','P',p_4,'Q',1,fluido); % temperatura punto 4 sulla cls [K]
T_5_1 = (T_4_cls + T_1_cli)/2; % temperatura punto 5 di tentativo [K]
delta_T_5 = 100;
delta_h_1 = 100;
% fine definizione dati iniziali

if scelta_rigeneratore == 0 % Rigeneratore assente
    while delta_T_5>=1 || delta_h_1>=10
       % calcolo prestazioni pompa
          v_1 = 1/(py.CoolProp.CoolProp.PropsSI('D','P',p_1,'H',h_1,fluido)); % volume specifico punto 1 [m^3/kg]
          W_pp = M_r*(v_1*(p_2-p_1))/eta_pp; % potenza assorbita dalla pompa [W]
          h_2 = h_1 + W_pp/M_r; % entalpia specifica punto 2 [J/kg]
       % fine calcolo prestazioni pompa

          h_3 = h_2;

       % calcolo prestazioni evaporatore
          h_4 = h_3 + Q_ev/M_r; % entalpia specifica punto 4 [J/kg]
       % fine calcolo prestazioni evaporatore

       % calcolo prestazioni espansore
        % Caduta di pressione all'alimentazione (su-su,1)
          h_su = h_4; % entalpia specifica all'alimentazione [J/kg]
          s_su = py.CoolProp.CoolProp.PropsSI('S','P',p_4,'H',h_4,fluido); % entropia all'alimentazione [J/kg]

          p_zero = p_4; % valore iniziale
          fun_1 = @(p_su_1)abs(M_r - A_su*py.CoolProp.CoolProp.PropsSI('D','S',s_su,'P',p_su_1,fluido)*sqrt(2*(h_su - py.CoolProp.CoolProp.PropsSI('H','S',s_su,'P',p_su_1,fluido))));
          options = optimoptions('fmincon', 'Display','none'); % rende silenziosa l'operazione
          p_su_1 = fmincon(fun_1,p_zero,[],[],[],[],0,p_4,[],options); % pressione sezione su,1 [Pa]
          h_su_1 = py.CoolProp.CoolProp.PropsSI('H','S',s_su,'P',p_su_1,fluido); % entalpia specifica sezione su,1 [J/kg]
          s_su_1 = py.CoolProp.CoolProp.PropsSI('S','S',s_su,'P',p_su_1,fluido); % entropia specifica sezione su,1 [J/kgK]
          v_su_1 = 1/(py.CoolProp.CoolProp.PropsSI('D','S',s_su,'P',p_su_1,fluido)); % volume specifico sezione su,1 [m^3/kg]

        % Perdite interne
          h_su_2 = h_su_1; % entalpia specifica sezione su,2 [J/kg]
          s_su_2 = s_su_1; % entropia specifica sezione su,2 [J/kgK]
          v_su_2 = v_su_1; % volume specifico sezione su,2 [m^3/kg]
          p_su_2 = p_su_1; % pressione sezione su,2 [Pa]

          cp_su_2 = py.CoolProp.CoolProp.PropsSI('C','H',h_su_2,'P',p_su_2,fluido); % calore specifico sezione su,2 [J/kgK]
          R_su_2_un = py.CoolProp.CoolProp.PropsSI('GAS_CONSTANT','H',h_su_2,'P',p_su_2,fluido); % costante universale dei gas [J/molK]
          massa_molare = py.CoolProp.CoolProp.PropsSI('M','H',h_su_2,'P',p_su_2,fluido); % massa molare [kg/mol]
          R_su_2 = R_su_2_un/massa_molare; % costante specifica del gas [J/kgK]
          gamma_su_2 = cp_su_2/(cp_su_2 - R_su_2); % esponente isoentropico

          p_ex_2 = p_5; % pressione sezione ex,2 [Pa]
          p_crit_leak = p_su_2*((2/(gamma_su_2 + 1))^(gamma_su_2/(gamma_su_2 - 1))); % pressione critica [Pa]
          p_thr_leak = max(p_ex_2,p_crit_leak); % pressione sezione di gola [Pa]
          v_thr_leak = 1/(py.CoolProp.CoolProp.PropsSI('D','P',p_thr_leak,'S',s_su_2,fluido)); % volume specifico sezione di gola [m^3/kg]
          h_thr_leak = py.CoolProp.CoolProp.PropsSI('H','P',p_thr_leak,'S',s_su_2,fluido); % entalpia specifica sezione di gola [J/kg]

          M_leak = A_leak/v_thr_leak*sqrt(2*(h_su_2 - h_thr_leak)); % portata persa [kg/s]

        % Potenza erogata (su,2-ex,2)
          M_in = M_r - M_leak; % portata interna [kg/s]
          N = (60*v_su_2*M_in)/V_s_exp; % velocità di rotazione espansore [giri/min]

          m = gamma_su_2 / (gamma_su_2 - eta_pol * (gamma_su_2 - 1)); % esponente della politropica
          p_ad = p_su_2*(1/rv)^m; % pressione sezione ad [Pa]
          v_ad = v_su_2*rv; % volume specifico sezione ad [m^3/kg]
          d_ad = 1 / v_ad; % densità sezione ad [kg/m^3]
          h_ad = py.CoolProp.CoolProp.PropsSI('H','P',p_ad,'D',d_ad,fluido); % entalpia specifica sezione ad [J/kg]

          W_in = M_in*((h_su_2 - h_ad) + v_ad*(p_ad - p_ex_2)); % potenza interna [W]
          W_loss = (2*pi*N)/60*T_loss; % potenza meccanica persa [W]

          W_sh = W_in - W_loss; % potenza meccanica all'albero [W]

        % Mescolamento adiabatico (ex,2-ex,1)
          h_ex_2 = h_su_2 - W_in/M_in; % entalpia specifica sezione ex,2 [J/kg]
          h_ex_1 = (h_su_2*M_leak + h_ex_2*M_in)/M_r; % entalpia specifica sezione ex,1 [J/kg]

        % Scarico (ex,1-ex)
          h_5 = h_ex_1; % entalpia specifica punto 5 [J/kg]
          T_5 = py.CoolProp.CoolProp.PropsSI('T','P',p_5,'H',h_5,fluido); % temperatura punto 5 [K]
       % fine calcolo prestazioni espansore

          h_6 = h_5; % entalpia specifica punto 6 [J/kg]
   
       % calcolo prestazioni condensatore
          h_1_ric = h_6 - Q_cond/M_r; % ricalcolo dell'entalpia punto 1 [J/kg]
       % fine calcolo prestazioni condensatore

       % calcolo del delta e aggiornamento variabili
          delta_T_5 = abs(T_5 - T_5_1);
          delta_h_1 = abs(h_1_ric - h_1);

          h_1 = h_1_ric;
          T_5_1 = T_5;   
    end
    T_1 = py.CoolProp.CoolProp.PropsSI('T','P',p_1,'H',h_1_ric,fluido); % temperatura punto 1 [K]
    T_2 = py.CoolProp.CoolProp.PropsSI('T','P',p_2,'H',h_2,fluido); % temperatura punto 2 [K]
    T_3 = T_2;
    T_4 = py.CoolProp.CoolProp.PropsSI('T','P',p_4,'H',h_4,fluido); % temperatura punto 4 [K]
    T_6 = T_5;

elseif scelta_rigeneratore ==1 % Rigeneratore presente
    while delta_T_5>=1 || delta_h_1>=10
       % calcolo prestazioni pompa
          v_1 = 1/(py.CoolProp.CoolProp.PropsSI('D','P',p_1,'H',h_1,fluido)); % volume specifico punto 1 [m^3/kg]
          W_pp = M_r*(v_1*(p_2-p_1))/eta_pp; % potenza assorbita dalla pompa [W]
          h_2 = h_1 + W_pp/M_r; % entalpia specifica punto 2 [J/kg]
          T_2 = py.CoolProp.CoolProp.PropsSI('T','P',p_2,'H',h_2,fluido); % temperatura punto 2 [K]
       % fine calcolo prestazioni pompa

       % calcolo prestazioni rigeneratore
          T_3 = T_2 + eta_rig*(T_5_1-T_2); % temperatura punto 3 [K]
          h_3 = py.CoolProp.CoolProp.PropsSI('H','P',p_3,'T',T_3,fluido); % entalpia specifica punto 3 [J/kg]  
       % fine calcolo prestazioni rigeneratore

       % calcolo prestazioni evaporatore
          h_4 = h_3 + Q_ev/M_r; % entalpia specifica punto 4 [J/kg]
       % fine calcolo prestazioni evaporatore

       % calcolo prestazioni espansore
        % Caduta di pressione all'alimentazione (su-su,1)
          h_su = h_4; % entalpia specifica all'alimentazione [J/kg]
          s_su = py.CoolProp.CoolProp.PropsSI('S','P',p_4,'H',h_4,fluido); % entropia all'alimentazione [J/kg]

          p_zero = p_4; % valore iniziale
          fun_1 = @(p_su_1)abs(M_r - A_su*py.CoolProp.CoolProp.PropsSI('D','S',s_su,'P',p_su_1,fluido)*sqrt(2*(h_su - py.CoolProp.CoolProp.PropsSI('H','S',s_su,'P',p_su_1,fluido))));
          options = optimoptions('fmincon', 'Display','none'); % rende silenziosa l'operazione
          p_su_1 = fmincon(fun_1,p_zero,[],[],[],[],0,p_4,[],options); % pressione sezione su,1 [Pa]
          h_su_1 = py.CoolProp.CoolProp.PropsSI('H','S',s_su,'P',p_su_1,fluido); % entalpia specifica sezione su,1 [J/kg]
          s_su_1 = py.CoolProp.CoolProp.PropsSI('S','S',s_su,'P',p_su_1,fluido); % entropia specifica sezione su,1 [J/kgK]
          v_su_1 = 1/(py.CoolProp.CoolProp.PropsSI('D','S',s_su,'P',p_su_1,fluido)); % volume specifico sezione su,1 [m^3/kg]

        % Perdite interne
          h_su_2 = h_su_1; % entalpia specifica sezione su,2 [J/kg]
          s_su_2 = s_su_1; % entropia specifica sezione su,2 [J/kgK]
          v_su_2 = v_su_1; % volume specifico sezione su,2 [m^3/kg]
          p_su_2 = p_su_1; % pressione sezione su,2 [Pa]

          cp_su_2 = py.CoolProp.CoolProp.PropsSI('C','H',h_su_2,'P',p_su_2,fluido); % calore specifico sezione su,2 [J/kgK]
          R_su_2_un = py.CoolProp.CoolProp.PropsSI('GAS_CONSTANT','H',h_su_2,'P',p_su_2,fluido); % costante universale dei gas [J/molK]
          massa_molare = py.CoolProp.CoolProp.PropsSI('M','H',h_su_2,'P',p_su_2,fluido); % massa molare [kg/mol]
          R_su_2 = R_su_2_un/massa_molare; % costante specifica del gas [J/kgK]
          gamma_su_2 = cp_su_2/(cp_su_2 - R_su_2); % esponente isoentropico

          p_ex_2 = p_5; % pressione sezione ex,2 [Pa]
          p_crit_leak = p_su_2*((2/(gamma_su_2 + 1))^(gamma_su_2/(gamma_su_2 - 1))); % pressione critica [Pa]
          p_thr_leak = max(p_ex_2,p_crit_leak); % pressione sezione di gola [Pa]
          v_thr_leak = 1/(py.CoolProp.CoolProp.PropsSI('D','P',p_thr_leak,'S',s_su_2,fluido)); % volume specifico sezione di gola [m^3/kg]
          h_thr_leak = py.CoolProp.CoolProp.PropsSI('H','P',p_thr_leak,'S',s_su_2,fluido); % entalpia specifica sezione di gola [J/kg]

          M_leak = A_leak/v_thr_leak*sqrt(2*(h_su_2 - h_thr_leak)); % portata persa [kg/s]

        % Potenza erogata (su,2-ex,2)
          M_in = M_r - M_leak; % portata interna [kg/s]
          N = (60*v_su_2*M_in)/V_s_exp; % velocità di rotazione espansore [giri/min]

          m = gamma_su_2 / (gamma_su_2 - eta_pol * (gamma_su_2 - 1)); % esponente della politropica
          p_ad = p_su_2*(1/rv)^m; % pressione sezione ad [Pa]
          s_ad = s_su_2; % entropia specifica nella sezione ad [J/kgK]
          v_ad = v_su_2*rv; % volume specifico sezione ad [m^3/kg]
          d_ad = 1 / v_ad; % densità sezione ad [kg/m^3]
          h_ad = py.CoolProp.CoolProp.PropsSI('H','P',p_ad,'D',d_ad,fluido); % entalpia specifica sezione ad [J/kg]

          W_in = M_in*((h_su_2 - h_ad) + v_ad*(p_ad - p_ex_2)); % potenza interna [W]
          W_loss = (2*pi*N)/60*T_loss; % potenza meccanica persa [W]

          W_sh = W_in - W_loss; % potenza meccanica all'albero [W]

        % Mescolamento adiabatico (ex,2-ex,1)
          h_ex_2 = h_su_2 - W_in/M_in; % entalpia specifica sezione ex,2 [J/kg]
          h_ex_1 = (h_su_2*M_leak + h_ex_2*M_in)/M_r; % entalpia specifica sezione ex,1 [J/kg]

        % Scarico (ex,1-ex)
          h_5 = h_ex_1; % entalpia specifica punto 5 [J/kg]
          T_5 = py.CoolProp.CoolProp.PropsSI('T','P',p_5,'H',h_5,fluido); % temperatura punto 5 [K]
       % fine calcolo prestazioni espansore
    
       % calcolo prestazioni rigeneratore
          h_6 = h_5 - (h_3 - h_2); % entalpia specifica punto 6 [J/kg]
       % fine calcolo prestazioni rigeneratore

       % calcolo prestazioni condensatore
          h_1_ric = h_6 - Q_cond/M_r; % ricalcolo dell'entalpia punto 1 [J/kg]
       % fine calcolo prestazioni condensatore

       % calcolo del delta e aggiornamento variabili
          delta_T_5 = abs(T_5 - T_5_1);
          delta_h_1 = abs(h_1_ric - h_1);

          h_1 = h_1_ric;
          T_5_1 = T_5;
    end

    T_1 = py.CoolProp.CoolProp.PropsSI('T','P',p_1,'H',h_1_ric,fluido); % temperatura punto 1 [K]
    T_4 = py.CoolProp.CoolProp.PropsSI('T','P',p_4,'H',h_4,fluido); % temperatura punto 4 [K]
    T_6 = py.CoolProp.CoolProp.PropsSI('T','P',p_6,'H',h_6,fluido); % temperatura punto 6 [K]

end
end