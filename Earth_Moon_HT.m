%% Earth - Moon Hohmann Transfer
% 
%% Setup
clc;close all;clear
% Radius of target orbits (WRT Earth)
ro = 1E7; % Radius of initial orbit [m]
rf = 3.7E8; % Radius of target orbit [m]

% Constants
G = 6.6743E-11; % Gravitational constant [N*m^2/kg^2]
Me = 5.972E24; % Mass of Earth [kg]
Mm = 7.346E22; % Mass of Moon [kg]
Re = 6.3781E6; % Radius of Earth [m]
Rm = 1.7371E6; % Radius of Moon [m]
Rem = 3.844E8; % Distance between Earth and Moon [m]
Tau = 2*pi/sqrt(G*(Me+Mm))*Rem^1.5; % Earth-Moon orbital period [s]
omega = 2*pi/Tau; % Angular velocity of Earth-Moon system [rad/s]
xe = Me/(Me+Mm)*Rem; % Earth distance from B0 [m]
xm = Mm/(Me+Mm)*Rem; % Moon distance from B0 [m]

% Delta V Burns
a = (ro+rf)/2;
Tau_i = 2*pi*ro/sqrt(G*Me/ro);
Tau_ht = 2*pi*a/sqrt(G*Me/a);
al = pi-Tau_ht/Tau_i*pi*(ro/rf)^1.5;
dV1 = sqrt(G*Me*(2/ro-1/a)) - sqrt(G*Mm/ro);

%% EoM
R1mag = @(x) sqrt((xe+x(1))^2+x(2)^2+x(3)^2);
R2mag = @(x) sqrt((x(1)-xm)^2+x(2)^2+x(3)^2);
Fe_co = @(x) G*Me/R1mag(x)^3;
Fm_co = @(x) G*Mm/R2mag(x)^3;

% State vector: [x y z x' y' z']
eqn = @(t,x) [x(4);x(5);x(6);
    2*omega*x(5)+omega^2*x(1)-Fe_co(x)*(x(1)+xe)-Fm_co(x)*(x(1)-xm);
    omega^2*x(2)-2*omega*x(4)-(Fe_co(x)+Fm_co(x))*x(2);
    -(Fe_co(x)+Fm_co(x))*x(3)];

%% Initial orbit
x01 = [ro;0;0;0;sqrt(G*Me/ro);0]; % Initial conditions, circular orbit
tspan1 = [0 2*pi*ro/sqrt(G*Me/ro)]; % Full orbit
opt1 = odeset('Events',@(t,x)eventalpha(t,x)); % Stop at alpha

sol1 = ode45(eqn,tspan1,x01,opt1);

%% Transfer orbit
% Transfer orbit burn
x0tr = sol1.y(:,end);
uv = x0tr(4:6)/norm(x0tr(4:6));
x0tr(4:6) = x0tr(4:6)+dV1*uv; % Add dV
tspantr = [0 2*pi*sqrt(a^3/(G*Me))]+sol1.x(end);

soltr = ode45(eqn,tspantr,x0tr,opttr);

%% Final orbit
% Circularization burn
x02 = soltr.y(:,end);


sol2 = ode45(eqn,tspan2,x02);

%% Plotting & Animation

