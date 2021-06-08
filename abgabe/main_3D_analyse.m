%% Becherwurf 3D Analyse
% Diese Programm dient zur Analyse der Flugbahn eines Balles unter gegebenen
% Anfangsbedingungen im dreidimensionalen Raum. Neben dem Einfluss der
% Anfangsbedingungen, welche nachfolgend verändert werden können, kann auch
% der Einfluss eines rotierenden Balls analysiert werden. Die
% Drehfrequenz wird auf die kurze Zeitdauer des Fluges als konstant
% angenommen und kann folgend bei den Parametereinstellungen verändert
% werden.


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
par.spinn_y = 3;       % Spinn um die y Achse in U/s, vereinfacht als konstant

% Anfangsbedingungen
x0 = 0;                 % Nullposition x in  m
y0 = 1.8;               % Nullposition y in m
z0 = 0;
v_0 = 6.5;               % Anfangsgeschwindigkeit des Balles in m/s
alpha_0 = 70*pi/180;    % Anfangswinkel des Balles in rad
alpha_y0 = 22*pi/180;    % Anfangswinkel des Balles um die y Achse in rad ausgehend von der x Achse
vy_0 = v_0*sin(alpha_0); % Anfangsgeschwindigkeit in y Richtung
vx_0 = v_0*cos(alpha_0)*cos(alpha_y0); % Anfangsgeschwindigkeit in y Richtung
vz_0 = v_0*cos(alpha_0)*sin(alpha_y0); % Anfangsgeschwindigkeit in z Richtung

f0 = [x0; y0; z0; vx_0; vy_0; vz_0];  % Anfangsbedingungen für x, y, z, v_x, v_y, v_z

% ODE Solver

tend = 5;             % Endzeit in s
tspan = [0 tend];        

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-8, 'Events', @(t,y) TischKontaktEvent(t,y,par)); % Select options for the solver
sol = ode45(@(t,y) equations(t,y,par), tspan, f0, options);     % Call the solver

%% Quick Visualization/Animation-------------------------------------------

% Diagramm xsol ysol-------------------------------------------------------
figure(13)
plot(sol.x,sol.y(2,:), '-o')
xlabel('solx');ylabel('soly');
title('Analyse Loesung');

%Diagramm 3D Ansicht-------------------------------------------------------
figure (14); clf(14);
plot3(sol.y(1,:),sol.y(3,:),sol.y(2,:), '-o'); grid on; hold on

% Becher zeichnen auf richtige Position
d=par.d_Becher;
l=par.l;
s=par.s_Becher;
alpha_Becher = asin(1/2);
X1 = l+d/2; Y1 = -0.9; Z1 = 0;
X2 = l+d/2+(d+s)*cos(alpha_Becher); Z2 = (d+s)*cos(alpha_Becher);
X3 = l+d/2+(d+s)*cos(alpha_Becher); Z3 = -(d+s)*cos(alpha_Becher);
X4 = l+d/2+2*(d+s)*cos(alpha_Becher); Z4 = 2*(d+s)*cos(alpha_Becher);
X5 = l+d/2+2*(d+s)*cos(alpha_Becher); Z5 = 0;
X6 = l+d/2+2*(d+s)*cos(alpha_Becher); Z6 = -2*(d+s)*cos(alpha_Becher);
h = 0.15;
r = 0.02;

[X,Y,Z] = cylinder(r);Z = -Z*h;
surf(X+X1,Y+Z1,Z);surf(X+X2,Y+Z2,Z);surf(X+X3,Y+Z3,Z);
surf(X+X4,Y+Z4,Z);surf(X+X5,Y+Z5,Z);surf(X+X6,Y+Z6,Z);
hold off

ylim([-1 1]);xlim([0 3]);zlim([-0.2 10]); 
xticks([0:0.5:3]);yticks([-1:0.5:1])
xlabel('Position x'); ylabel('Position z'); zlabel('Position y');
title('Flugbahn 3D');

% Diagramm XY--------------------------------------------------------------
figure(15);clf(15);
plot(sol.y(1,:),sol.y(2,:), '-o');hold on;
r1 = rectangle('Position', [l, -0.15, d, 0.15],'FaceColor','r');
r2 = rectangle('Position', [l+(d+s)*cos(alpha_Becher), -0.15, d, 0.15],'FaceColor','r');
r3 = rectangle('Position', [l+2*(d+s)*cos(alpha_Becher), -0.15, d, 0.15],'FaceColor','r');
x = [0 2.5];
y = [-0.15 -0.15];
line(x,y,'Color','k','LineWidth',1);
xlim([0 2.6]);ylim([-0.2;10]);
hold off
xlabel('Position x');ylabel('Position y');
title('Flugbahn XY Ebene');

% Diagramm XZ--------------------------------------------------------------
figure(16);clf(16);
plot(sol.y(1,:),sol.y(3,:), '-o');hold on;
circle(X1,Z1,d/2);circle(X2,Z2,d/2);circle(X3,Z3,d/2);
circle(X4,Z4,d/2);circle(X5,Z5,d/2);circle(X6,Z6,d/2);
rectangle('Position', [0, -0.5, 2.5, 1],'EdgeColor','k');
xlim([0 2.6]);ylim([-1 1]);
hold off
xlabel('Position x');ylabel('Position z');
title('Flugbahn XZ Ebene');


%% Animation 3D------------------------------------------------------------
figure(17); clf(17);
title('Animation Wurf');
for k=1:1:length(sol.x)
    subplot(2,2,1);
    plot(sol.y(1,k),sol.y(2,k),'k.','MarkerSize',20);hold on;
    plot(sol.y(1,:),sol.y(2,:), '--k');
    r1 = rectangle('Position', [l, -0.15, d, 0.15],'FaceColor','r');
    r2 = rectangle('Position', [l+(d+s)*cos(alpha_Becher), -0.15, d, 0.15],'FaceColor','r');
    r3 = rectangle('Position', [l+2*(d+s)*cos(alpha_Becher), -0.15, d, 0.15],'FaceColor','r');
    x = [0 2.5];
    y = [-0.15 -0.15];
    line(x,y,'Color','k','LineWidth',1);
    hold off
    xlim([0 2.6]);ylim([-0.2;10]);
    xlabel('Position x'); ylabel('Position y');
    title('XY')

    subplot(2,2,2);
    plot(sol.y(1,k),sol.y(3,k),'k.','MarkerSize',20);hold on;
    plot(sol.y(1,:),sol.y(3,:), '--k')
    circle(X1,Z1,d/2);circle(X2,Z2,d/2);circle(X3,Z3,d/2);
    circle(X4,Z4,d/2);circle(X5,Z5,d/2);circle(X6,Z6,d/2);
    rectangle('Position', [0, -0.5, 2.5, 1],'EdgeColor','k');
    xlim([0 2.6]);ylim([-1 1]);
    hold off
    xlabel('Position x'); ylabel('Position z');
    title('XZ')

    subplot(2,2,[3,4]);
    plot3(sol.y(1,k),sol.y(3,k),sol.y(2,k),'k.', 'MarkerSize',20);grid on; hold on
    plot3(sol.y(1,:),sol.y(3,:),sol.y(2,:), '--k');
    [X,Y,Z] = cylinder(r);Z = -Z*h;
    surf(X+X1,Y+Z1,Z,'FaceColor','r');surf(X+X2,Y+Z2,Z,'FaceColor','r');surf(X+X3,Y+Z3,Z,'FaceColor','r');
    surf(X+X4,Y+Z4,Z,'FaceColor','r');surf(X+X5,Y+Z5,Z,'FaceColor','r');surf(X+X6,Y+Z6,Z,'FaceColor','r');
    hold off
    ylim([-1 1]);xlim([0 3]);zlim([-0.2 10]); 
    xticks([0:0.5:3]);yticks([-1:0.5:1])
    xlabel('Position x'); ylabel('Position z'); zlabel('Position y');
    title('XYZ')
    pause(1/10);
end

%% FUNKTIONEN---------------------------------------------------------------


% Explizite Differentialgleichung
function f = equations(~,y,par)

% Unpack parameters
m_Ball = par.m_Ball;
g = par.g;
cw = par.cw;
d_Ball = par.d_Ball;
rho = par.rho;
omega_y = par.spinn_y*2*pi;

f(1) = y(4);
f(2) = y(5);
f(3) = y(6);
f(4) = -((1/8)*rho*(d_Ball^2)*pi*cw*sqrt((y(4)^2)+(y(5)^2)+(y(6)^2))*y(4))/m_Ball...
    +(rho*(d_Ball^3)*((pi^2)/8)*omega_y*y(6))/m_Ball;
f(5) = -((1/8)*rho*(d_Ball^2)*pi*cw*sqrt((y(4)^2)+(y(5)^2)+(y(6)^2))*y(5))/m_Ball-g;
f(6) = -((1/8)*rho*(d_Ball^2)*pi*cw*sqrt((y(4)^2)+(y(5)^2)+(y(6)^2))*y(6))/m_Ball...
    -(rho*(d_Ball^3)*((pi^2)/8)*omega_y*y(4))/m_Ball;

% Die Formel für den Magnuseffekt des drehenden Balls wurde aus
%https://www.grc.nasa.gov/www/k-12/VirtualAero/BottleRocket/airplane/beach.html
%entnommen

f = f';
end

% Event Tisch/Becherkontakt = ODE Abbruch
function [value,isterminal,direction] = TischKontaktEvent(~,y,par)

value = y(2)-par.d_Ball;
isterminal = 1;
direction = 0; 

end

% Funktion Kreis plotten
function circle(x,y,r)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)
ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
plot(x+xp,y+yp,'Color','r');
end
