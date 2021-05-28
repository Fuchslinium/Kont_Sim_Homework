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

f0 = [x0; y0; v_0*cos(alpha_0); v_0*sin(alpha_0)];  % Anfangsbedingungen für x, y, v_x, v_y 

% ODE Solver

% FUNKTIONEN---------------------------------------------------------------
function f = equations(t,y,par)

% Parameter entpacken
m_Ball = par.m_Ball;
g = par.g;
cw = par.cw;
d_Ball = par.d_Ball;
rho = par.rho;

A = d_Ball^2*pi/4;    % Querschnittsfläche des Balles

f(1) = y(3);
f(2) = y(4);
f(3) = ((1/8)*rho*d_Ball^2*pi*cw*sqrt(y(3)+y(4))*y(3))/m_Ball;
f(4) = (((1/8)*rho*d_Ball^2*pi*cw*sqrt(y(3)+y(4))*y(3))/m_Ball)-g;

end