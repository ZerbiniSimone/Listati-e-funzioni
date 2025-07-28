clear; clc; close all;

fluido = 'R290';     % fluido di lavoro
p_4 = 3192000;
p_5 = 1263158;
h_4 = 6.598526557827813e+05;
%delta_p_r = 0.03;    % perdita di pressione relativa
M_r = 0.13216;      % portata effettivamente erogata dalla pompa [kg/s]
A_leak = 4.6e-6;     % area della sezione di perdita [m^2]
A_su = 27.43e-6;     % area della sezione di aspirazione [m^2]
rv = 2.5;           % rapporto volumetrico della macchina
V_s_exp = 100e-6;  % cilindrata dell'espansore [m^3]
T_loss = 0.47;       % coppia meccanica persa [Nm]
eta_pol = 1;

% calcolo prestazioni espansore
        % Caduta di pressione all'alimentazione (su-su,1)
          h_su = h_4; % entalpia specifica all'alimentazione [J/kg]
          s_su = py.CoolProp.CoolProp.PropsSI('S','P',p_4,'H',h_4,fluido); % entropia all'alimentazione [J/kg]

          p_zero = p_4; % valore iniziale
          fun_1 = @(p_su_1)abs(M_r - A_su*py.CoolProp.CoolProp.PropsSI('D','S',s_su,'P',p_su_1,fluido)*sqrt(2*(h_su - py.CoolProp.CoolProp.PropsSI('H','S',s_su,'P',p_su_1,fluido))));
          options = optimoptions('fmincon', 'Display','none'); % rende silenziosa l'operazione
          p_su_1 = fmincon(fun_1,p_zero,[],[],[],[],0,p_4,[],options); % pressione sezione su,1 [Pa]
          h_su_1 = h_su; % entalpia specifica sezione su,1 [J/kg]
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
     

        
          