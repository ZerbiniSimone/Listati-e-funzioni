p_5 = 1263158;
fluido = 'R290';

h_test_1 = 598426.79;
T_test_1 = py.CoolProp.CoolProp.PropsSI('T','P',p_5,'H',h_test_1,fluido)
X_test_1= py.CoolProp.CoolProp.PropsSI('Q', 'P', p_5, 'H', h_test_1, fluido)

h_test_2 = 608131.95;
T_test_2 = py.CoolProp.CoolProp.PropsSI('T','P',p_5,'H',h_test_2,fluido)
X_test_2 = py.CoolProp.CoolProp.PropsSI('Q', 'P', p_5, 'H', h_test_2, fluido)