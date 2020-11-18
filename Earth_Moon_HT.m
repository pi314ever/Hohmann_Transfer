%% Earth - Moon Hohmann Transfer
% 
%% Setup
clc;close all;clear
% Radius of target orbits (WRT Earth)
ro = 1E7; % Radius of initial orbit [m]
rf = 3.5E8; % Radius of target orbit [m]

% Constants
G = 6.6743E-11; % Gravitational constant [N*m^2/kg^2]
Me = 5.972E24; % Mass of Earth [kg]
Mm = 7.346E22; % Mass of Moon [kg]
Re = 6.3781E6; % Radius of Earth [m]
Rm = 1.7371E6; % Radius of Moon [m]
Rem = 3.844E8; % Distance between Earth and Moon [m]
Tau = 2*pi/sqrt(G*(Me+Mm))*Rem^1.5; % Earth-Moon orbital period [s]
omega = 2*pi/Tau; % Angular velocity of Earth-Moon system [rad/s]
xe = Mm/(Me+Mm)*Rem; % Earth distance from B0 [m]
xm = Me/(Me+Mm)*Rem; % Moon distance from B0 [m]

% Delta V Burns
a = (ro+rf)/2;
al = pi-pi*(a/Rem)^1.5;
dV1 = sqrt(G*Me*(2/ro-1/a)) - sqrt(G*Me/ro);
% dV1 = dV1/3.5;

%% EoM
R1mag = @(x) sqrt((xe+x(1,:)).^2+x(2,:).^2+x(3,:).^2);
R2mag = @(x) sqrt((x(1,:)-xm).^2+x(2,:).^2+x(3,:).^2);
Fe_co = @(x) G*Me/R1mag(x).^3;
Fm_co = @(x) G*Mm/R2mag(x).^3;

% State vector: [x y z x' y' z']
eqn = @(t,x) [x(4);x(5);x(6);
    2*omega*x(5)+omega^2*x(1)-Fe_co(x)*(x(1)+xe)-Fm_co(x)*(x(1)-xm);
    omega^2*x(2)-2*omega*x(4)-(Fe_co(x)+Fm_co(x))*x(2);
    -(Fe_co(x)+Fm_co(x))*x(3)];

%% Initial orbit
x01 = [-ro-xe;0;0;0;-sqrt(G*Me/ro);0]; % Initial conditions, circular orbit
tspan1 = [0 2*pi*ro/sqrt(G*Me/ro)]*8; % Full orbit
opt1 = odeset('Events',@(t,x)eventalpha(t,x,al,xe,xm,Re,Rm)); % Stop at alpha

sol1 = ode45(eqn,tspan1,x01,opt1);

%% Transfer orbit
% Transfer orbit burn
x0tr = sol1.y(:,end);
uv = x0tr(4:6)/norm(x0tr(4:6));
x0tr(4:6) = x0tr(4:6)+dV1*uv; % Add dV
tspantr = [0 2*pi*sqrt(a^3/(G*Me))]*4+sol1.x(end);
opttr = odeset('Events',@(t,x)eventmoon(t,x,xe,xm,Re,Rm)); % Stop at closest approach

soltr = ode45(eqn,tspantr,x0tr,opttr);

%% Final orbit
% Circularization burn
x02 = soltr.y(:,end);
uv2 = x02(4:6)/norm(x02(4:6));
x02(4:6) = sqrt(G*Mm/R2mag(x02))*uv2;
tspan2 = [0 pi/6*rf/sqrt(G*Me/rf)/100]*100+soltr.x(end);
opt2 = odeset('Events',@(t,x)eventhit(t,x,xe,xm,Re,Rm,Rem)); % Stop if hit Moon or Earth

sol2 = ode45(eqn,tspan2,x02,opt2);

%% What if no dV2
x02_nodV = soltr.y(:,end);
tspan2_nodV = tspan2;
opt2_nodV = opt2;

sol2_nodV = ode45(eqn,tspan2,x02,opt2);

%% Deval for constant time step
h = 10; % Time step [s]
t1 = sol1.x(1):h:sol1.x(end-1);
ttr = soltr.x(1):h:soltr.x(end-1);
t2 = sol2.x(1):h:sol2.x(end);
t2_nodV = sol2_nodV.x(1):h:sol2_nodV.x(end);
dsol1 = deval(sol1,t1);
dsoltr = deval(soltr,ttr);
dsol2 = deval(sol2,t2);
dsol2_nodV = deval(sol2_nodV.t2_nodV);

%% Extract coordinates
x1 = dsol1(1,:);
xtr = dsoltr(1,:);
x2 = dsol2(1,:);
y1 = dsol1(2,:);
ytr = dsoltr(2,:);
y2 = dsol2(2,:);
z1 = dsol1(3,:);
ztr = dsoltr(3,:);
z2 = dsol2(3,:);

x = [x1 xtr x2];
y = [y1 ytr y2];
z = [z1 ztr z2];
t = [t1 ttr t2];

%% Barycenter Trajectory plot
[xs,ys,zs] = sphere;
% Earth coords
Xe = xs*Re-xe;
Ye = ys*Re;
Ze = zs*Re;
% Moon coords
Xm = xs*Rm+xm;
Ym = ys*Rm;
Zm = zs*Rm;
figure
plot3(x,y,z,'linewidth',2)
hold on
plot3(-xe,0,0,'kx','markersize',5)
plot3(xm,0,0,'kx','markersize',5)
plot3(x(1),y(1),z(1),'x','markersize',5)
plot3(xtr(1),ytr(1),ztr(1),'x','markersize',5)
plot3(x2(1),y2(1),z2(1),'x','markersize',5)
surf(Xe,Ye,Ze);surf(Xm,Ym,Zm)
legend('Trajectory','Earth','Start','dV1','dV2')
xlabel('x')
ylabel('y')
zlabel('z')
axis equal

%% 2D Plot
th = 0:0.01:2*pi;
xet = Re*cos(th)-xe;
yet = Re*sin(th);
xmt = Rm*cos(th)+xm;
ymt = Rm*sin(th);
figure
plot(xet,yet,'k--',xmt,ymt,'k--')
hold on
plot(x(1),y(1),'x')
plot(xtr(1),ytr(1),'x')
plot(x2(1),y2(1),'x')
legend('E','M','Start','dV1','dV2')
axis equal
comet(x,y)

%% Newtonian plot
ang = t*omega;
figure
plot(x.*cos(ang)-y.*sin(ang),x.*sin(ang)+y.*cos(ang))
axis equal


%% Animation


