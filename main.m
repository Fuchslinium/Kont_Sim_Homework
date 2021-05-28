%% Becherwurf

% Parameter
par.m_Ball = 2.7e-3;   % Masse des Balles in kg
par.d_Ball = 0.04;     % Durchmesser des Balles in m
par.d_Becher = 0.095;  % Durchmesser des Bechers in m
par.s_Becher = 5e-3;   % Abstand der Becher in m
par.g = 9.81;          % Fallbeschleunigung in m/s²
par.l = 2.1;           % Länge des Tisches in m
par.cw = 0.47;         % Luftwiderstandsbeiwert des Balles
par.rho = 1.2;         % Dichte der Luft in kg/m³

% Anfangsbedingungen
x0 = 0;                % Nullposition x in  m
y0 = 0.05;             % Nullposition y in m
v_0 = 2;               % Anfangsgeschwindigkeit des Balles in m/s
alpha_0 = 20*pi/180;   % Anfangswinkel des Balles in rad

AB = [x0; y0; v_0*cos(alpha_0); v_0*sin(alpha_0)];  % Anfangsbedingungen für x, y, v_x, v_y 

% ODE Solver