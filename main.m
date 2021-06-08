%% Becherwurf 
clear all
clc
clf

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
z0 = 0;                 % Nullposition z in m
v_0 = 3;               % Anfangsgeschwindigkeit des Balles in m/s
alpha_0 = 15*pi/180;    % Anfangswinkel des Balles in rad
alpha_y0 = 5*pi/180;    % Anfangswinkel des Balles um die y Achse in rad ausgehend von der x Achse
vy_0 = v_0*sin(alpha_0); % Anfangsgeschwindigkeit in y Richtung
vx_0 = v_0*cos(alpha_0)*cos(alpha_y0); % Anfangsgeschwindigkeit in y Richtung
vz_0 = v_0*cos(alpha_0)*sin(alpha_y0); % Anfangsgeschwindigkeit in z Richtung

f0 = [x0; y0; z0; vx_0; vy_0; vz_0];  % Anfangsbedingungen für x, y, z, v_x, v_y, v_z

% ODE Solver

tend = 1;             % Endzeit in s
tspan = [0 tend];     % Zeitvektor in s   

options = odeset('RelTol', 1e-5, 'AbsTol', 1e-8);           % Select options for the solver
sol = ode45(@(t,y) equations(t,y,par), tspan, f0, options);     % Call the solver

% Plot der Flugkurve
figure(1)
plot(sol.y(1,:),sol.y(2,:), '-o')
hold on
ylim([0 2*y0])
xlabel('x-Position in m')
ylabel('y-Position in m')
hold on

% Plot Simulink
results_sim = sim('Hausuebung_Simulink');
x_sim = results_sim.x_sim.Data;
y_sim = results_sim.y_sim.Data;
% figure(9);clf(9);
plot(x_sim(:), y_sim(:), '-o')
legend('Position des Balles, Matlab Solver', 'Lösung mit Simulink')
title("Verlgeich Matlab und Simulink")
hold off

%% Parameteranalyse
% varying alpha_0

figure(2)
alpha_vary = alpha_0-10*pi/180:5*pi/180:alpha_0+10*pi/180;
for index = 1:length(alpha_vary)
    
    v_0x = v_0*cos(alpha_vary(index));
    v_0y = v_0*sin(alpha_vary(index));
   
    f1 = [x0; y0; z0; v_0x; v_0y; vz_0];  % Anfangsbedingungen für x, y, v_x, v_y 

    sol = ode45(@(t,y) equations(t,y,par), tspan, f1, options);
    plot(sol.y(1,:),sol.y(2,:));
    hold on
end

ylim([0 2*y0])
xlabel('x-Position in m')
ylabel('y-Position in m')
s_l = par.l + par.s_Becher + (3/2)*par.d_Becher;    % Mittlerer Becher als Ziel
l = line([s_l s_l],[0 5*y0],'Color','red','LineStyle','--');
legend('\alpha_0-10°','\alpha_0-5°','\alpha_0','\alpha_0+5°','\alpha_0+10°', 'Ziel')
title("varying \alpha_0")
hold off

%varying v0
figure(3)

v_vary = v_0-1 : 0.5 : v_0+1;
for index = 1:length(v_vary)
    
    v_0x = v_vary(index)*cos(alpha_0);
    v_0y = v_vary(index)*sin(alpha_0);
   
    f1=[x0; y0; z0; v_0x; v_0y; vz_0];  % Anfangsbedingungen für x, y, v_x, v_y 

    sol = ode45(@(t,y) equations(t,y,par), tspan, f1, options);
    plot(sol.y(1,:),sol.y(2,:));
    hold on
end

ylim([0 2*y0])
xlabel('x-Position in m')
ylabel('y-Position in m')
s_l = par.l + par.s_Becher + (3/2)*par.d_Becher;    % Mittlerer Becher als Ziel
l = line([s_l s_l],[0 5*y0],'Color','red','LineStyle','--');
legend('v0-1','v0-0.5','v0', 'v0+0.5','v0+1','Ziel')
title("varying v_0")
hold off

%% Visualisierung

% Anfangslösung berechnen:

sol = ode45(@(t,y) equations(t,y,par), tspan, f0, options);     % Call the solver

soly_new= sol.y(2,:); %extracting solution fom sol.y matrix
soly_new=soly_new(soly_new>0); %cutting solution vector to stop animation when ball passes cup level

solx_new= sol.y(1,:);

solx_new(:,length(soly_new)+1:end)=[]; %matching x vector size to y vector size

distance_t= interp1(soly_new,solx_new,0,"linear", "extrap"); % travelled distance when passing the top level of the cups

%adding extrapolated value at cup level so animation stops at cup level
soly_new(end+1)=0;
solx_new(end+1) = distance_t;

figure(4)

%plot axis limits and scaling
xlim([0,max(sol.y(1,:))]);
ylim([-0.5, max(sol.y(2,:))+0.5]);
daspect([1,1,1]);

%added parameters
h_becheranim=0.3; %cup-height for animation
danim = par.d_Ball;
ranim = danim/2;

%ball
h = rectangle('Position',[solx_new(1)-ranim soly_new(1)-ranim danim danim],'Curvature',[1,1]);
h.LineWidth = 3;

%cups 
cup1_anim = rectangle("Position",[par.l 0-h_becheranim par.d_Becher h_becheranim] );
cup2_anim = rectangle("Position",[par.l+par.s_Becher+par.d_Becher 0-h_becheranim par.d_Becher h_becheranim] );
cup3_anim = rectangle("Position",[par.l+2*par.s_Becher+2*par.d_Becher 0-h_becheranim par.d_Becher h_becheranim] );
cup1_anim.EdgeColor = "w";
cup2_anim.EdgeColor = "w";
cup3_anim.EdgeColor = "w";
cup1_anim.FaceColor = "r";
cup2_anim.FaceColor = "r";
cup3_anim.FaceColor = "r";


%Tabletop
tabletop_anim = rectangle("Position",[0 -0.01-h_becheranim par.l+3*par.s_Becher+3*par.d_Becher 0.01]);
tabletop_anim.EdgeColor = [0.2 0 0];
tabletop_anim.FaceColor = [0.2 0 0];
tabletop_anim.LineWidth = 3;

%trajectory
a_line = animatedline;


% animation of Ball and trajectory
for m =  2:length(solx_new)
    h.Position = [solx_new(m)-ranim;soly_new(m)-ranim;danim;danim];
   
    addpoints(a_line,solx_new(m),soly_new(m));
    drawnow
   
    pause(1/length(soly_new));
      
end

% event information
pause (0.1)

if distance_t >= par.l+3*par.s_Becher+3*par.d_Becher %too short event
   t_line= line([par.l+3*par.s_Becher+3*par.d_Becher,distance_t],[0,0]);
   txt = '\rightarrow too far!';
   message = text(2.2,0.1,txt);
   message.FontSize=14;
   
elseif distance_t <= par.l %too far event
    t_line= line([par.l,distance_t],[0,0]);
    txt = '\leftarrow too short!';
    message = text(1.8,0.1,txt);
    message.FontSize=14;
   
else
    txt = 'good throw!'; %within cup range event
    message = text(2,0.1,txt);
    message.FontSize=14;
    
end


% zoomed in Result
pause(0.1)
xlim([1.5 2.7])
ylim([-0.5 0.5])

%% Quick Visualization-----------------------------------------------------
t = tspan(1):t_step:tspan(2);
results = deval(sol, t);

figure(5)
plot3(results(1,:),results(3,:),results(2,:),'LineWidth', 4)
xlabel('x-Position in m')
ylabel('z-Position in m')
ylim([0 0.5])
zlabel('y-Position in m')
zlim([0 2*y0])
title("3D Plot")
grid on

% figure(6)
% plot(results(1,:),results(2,:), '-o')
% ylim([0 2*y0]);xlim([0 2.8])
% xlabel('x-Position in m')
% ylabel('y-Position in m')
% 
% 
% figure(7)
% plot(results(1,:),results(3,:), '-o')
% ylim([0 1]);xlim([0 2.8])
% xlabel('x-Position in m')
% ylabel('z-Position in m')
% 
% figure(8)
% plot(sol.y(1,:),sol.y(2,:), '-o')
% xlabel('solx')
% ylabel('soly')



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
f(4) = -((1/8)*rho*(d_Ball^2)*pi*cw*sqrt((y(4)^2)+(y(5)^2)+(y(6)^2))*y(4))/m_Ball;
f(5) = -((1/8)*rho*(d_Ball^2)*pi*cw*sqrt((y(4)^2)+(y(5)^2)+(y(6)^2))*y(5))/m_Ball-g;
f(6) = -((1/8)*rho*(d_Ball^2)*pi*cw*sqrt((y(4)^2)+(y(5)^2)+(y(6)^2))*y(6))/m_Ball;

f = f';
end
