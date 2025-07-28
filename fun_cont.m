function controllo = fun_cont(M_r, p_1, p_2, p_3, p_4, p_5, p_6, X, T_4_sp, T_1_sp, eta_rig, eta_pp, A_leak, A_su, T_loss, V_s_exp, rv, fluido, scelta_rigeneratore)

% FUNZIONE CHE CALCOLA LO SCARTO TRA LE TEMPERATURE CALCOLATE E QUELLE DI SET POINT

[~,~,~,T_1,~,~,T_4,~,~,~] = fun_orc(M_r,p_1,p_2,p_3,p_4,p_5,p_6,X,eta_rig,eta_pp,A_leak,A_su,T_loss,V_s_exp,rv,fluido,scelta_rigeneratore);

controllo_T_4 = abs(T_4 - T_4_sp); % scarto della T_4
controllo_T_1 = abs(T_1 - T_1_sp); % scarto della T_1

controllo = sqrt(controllo_T_1^2 + controllo_T_4^2)/2; % scarto complessivo

end