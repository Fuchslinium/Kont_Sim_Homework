%% Becherwurf Lorenz Version
clearvars
clc
% Parameter
par.m_Ball = 2.7e-3;   % Masse des Balles in kg
par.d_Ball = 0.04;     % Durchmesser des Balles in m
par.d_Becher = 0.095;  % Durchmesser des Bechers in m
par.s_Becher = 5e-3;   % Abstand der Becher in m
par.g = 9.81;          % Fallbeschleunigung in m/s²
par.l = 2.1;           % Länge des Tisches in m
par.cw = 0.47;         % Luftwiderstandsbeiwert des Balles
par.rho = 1.2;         % Dichte der Luft in kg/m³
t_step = 1/80;

% Anfangsbedingungen
x0 = 0;                 % Nullposition x in  m
y0 = 1.8;               % Nullposition y in m
z0 = 0;
v_0 = 9;               % Anfangsgeschwindigkeit des Balles in m/s
alpha_0 = 20*pi/180;    % Anfangswinkel des Balles in rad
alpha_y0 = 5*pi/180;    % Anfangswinkel des Balles um die y Achse in rad ausgehend von der x Achse
vy_0 = v_0*sin(alpha_0); % Anfangsgeschwindigkeit in y Richtung
vx_0 = v_0*cos(alpha_0)*cos(alpha_y0); % Anfangsgeschwindigkeit in y Richtung
vz_0 = v_0*cos(alpha_0)*sin(alpha_y0); % Anfangsgeschwindigkeit in z Richtung

f0 = [x0; y0; z0; vx_0; vy_0; vz_0];  % Anfangsbedingungen für x, y, z, v_x, v_y, v_z

% ODE Solver

tend = 5;             % Endzeit in s
tspan = [0 tend];        

options = odeset('RelTol', 1e-5, 'AbsTol', 1e-8);           % Select options for the solver
sol = ode45(@(t,y) equations(t,y,par), tspan, f0, options);     % Call the solver

%% Quick Visualization/Animation-------------------------------------------
t = tspan(1):t_step:tspan(2);
results = deval(sol, t);


figure(1)
plot(results(1,:),results(2,:), '-o')
ylim([0 5*y0]);xlim([0 2.8])
xlabel('x-Position in m')
ylabel('y-Position in m')


figure(2)
plot(results(1,:),results(3,:), '-o')
ylim([0 1]);xlim([0 2.8])
xlabel('x-Position in m')
ylabel('z-Position in m')

figure(3)
plot(sol.x,sol.y(2,:), '-o')
xlabel('solx')
ylabel('soly')


%% FUNKTIONEN---------------------------------------------------------------


% Explizite Differentialgleichung
function f = equations(~,y,par)

% Unpack parameters
m_Ball = par.m_Ball;
g = par.g;
cw = par.cw;
d_Ball = par.d_Ball;
rho = par.rho;

f(1) = y(4);
f(2) = y(5);
f(3) = y(6);
f(4) = ((1/8)*rho*(d_Ball^2)*pi*cw*sqrt((y(4)^2)+(y(5)^2)+(y(6)^2))*y(4))/m_Ball;
f(5) = ((1/8)*rho*(d_Ball^2)*pi*cw*sqrt((y(4)^2)+(y(5)^2)+(y(6)^2))*y(5))/m_Ball-g;
f(6) = ((1/8)*rho*(d_Ball^2)*pi*cw*sqrt((y(4)^2)+(y(5)^2)+(y(6)^2))*y(6))/m_Ball;

f = f';
end
