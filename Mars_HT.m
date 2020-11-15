%% Mars Hohmann Transfer 
% Mars-centered Hohmann Transfer between two circular orbits
%% Setup
clc;close all;clear
% Radius of orbits
R1 = 5E6; % Altitude 1 [m]
R2 = 7E6; % Altitude 2 [m]
h = 10; % Time step [s]

% Constants
Mm = 6.39E23; % Mass of Mars [kg]
Rm = 3.3895E6; % Radius of Mars [m]
G = 6.6743E-11; % Gravitational constant [N*m^2/kg^2]

% EoM (r,theta,rd,thetad)
eqn = @(t,x) [x(3);x(4);x(4)^2*x(1)-G*Mm/x(1)^2;-2*x(3)*x(4)/x(1)];

% Velocity calculations
a = (R1+R2)/2;
dV1 = sqrt(G*Mm*(2/R1-1/a)) - sqrt(G*Mm/R1);
dV2 = sqrt(G*Mm/R2) - sqrt(G*Mm*(2/R2-1/a));

% Initial conditions 
x01 = [R1;0;-pi/2;sqrt(G*Mm/R1)/R1];
tspan1 = [0 2*pi*R1/sqrt(G*Mm/R1)];

%% Initial orbit
sol1RK = RK4(eqn,tspan1,x01,h);
sol1ABM = ABM4(eqn,tspan1,x01,h);

%% Transfer orbit
% Add dV
x0trRK = sol1RK.y(:,end);
x0trRK(4) = x0trRK(4)+dV1/R1;
x0trABM = sol1ABM.y(:,end);
x0trABM(4) = x0trABM(4)+dV1/R1;

tspantr = [0 pi*sqrt(a^3/(G*Mm))]+tspan1(2);

soltrRK = RK4(eqn,tspantr,x0trRK,h);
soltrABM = ABM4(eqn,tspantr,x0trABM,h);

%% Final orbit
% Add more dV
x02RK = soltrRK.y(:,end);
x02RK(4) = x02RK(4)+dV2/R2;
x02ABM = soltrABM.y(:,end);
x02ABM(4) = x02RK(4)+dV2/R2;

tspan2 = [0 pi*R2/sqrt(G*Mm/R2)]+tspantr(2);

sol2RK = RK4(eqn,tspan2,x02RK,h);
sol2ABM = ABM4(eqn,tspan2,x02ABM,h);

%% Plotting
% Combining 
t_RK = [sol1RK.x soltrRK.x sol2RK.x];
t_ABM = [sol1ABM.x soltrABM.x sol2RK.x];
R_RK = [sol1RK.y(1,:) soltrRK.y(1,:) sol2RK.y(1,:)];
Th_RK = [sol1RK.y(2,:) soltrRK.y(2,:) sol2RK.y(2,:)];
R_ABM = [sol1ABM.y(1,:) soltrABM.y(1,:) sol2ABM.y(1,:)];
Th_ABM = [sol1ABM.y(2,:) soltrABM.y(2,:) sol2ABM.y(2,:)];

x_RK = R_RK.*cos(Th_RK);
y_RK = R_RK.*sin(Th_RK);
x_ABM = R_ABM.*cos(Th_ABM);
y_ABM = R_ABM.*sin(Th_ABM);

% Mars circle
th = 0:0.01:2*pi;
xm = Rm*cos(th);
ym = Rm*sin(th);

% Trajectory plot
figure
plot(x_RK,y_RK,x_ABM,y_ABM,'--',xm,ym,'k--')
legend('RK','ABM','Mars')

