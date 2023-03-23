clear all; close all; clc;

m2 = 7.347673 * 10^(22);    % Mass of Moon (kg)
m1 = 5.97219 * 10^(24);     % Mass of Earth (kg)
mu = m2 / (m2 + m1);        % Earth-Moon Mass Ratio

r12 = 3.825 * 10^(5);       % Earth-Moon Distance (km)

C = 2.98799707128764;       % L4 Jacobi Constant


%% Calculations

% Short / Long Period about L4

s1 = 0.29820789544347032;
s2 = 0.95450094347526815;

dyn_short = dynICs_short(384.388174,0);
tau_short = 0:0.001:2*pi/s2;
x0 = [384.388174;
      0;
      -dyn_short(1);
      dyn_short(2)];
out_short = short(tau_short,x0);

dyn_long = dynICs_long(384.388174,0);
tau_long = 0:0.001:2*pi/s1;
x0 = [384.388174;
      0;
      -dyn_long(1);
      dyn_long(2)];
out_long = long(tau_long,x0);

% L4 ZVC
x1 = -mu;
x2 = 1 - mu;

[X,Y] = meshgrid(0.475:0.00001:0.5,0.855:0.00001:0.875);
Z = X.^2 + Y.^2 + 2.*(1-mu)./sqrt((X-x1).^2+Y.^2)...
    + 2.*mu./sqrt((X-x2).^2+Y.^2);



% Transfer
SOI_earth = 924000;     % km
SOI_moon = 66100;       % km
center_earth = [0 0];
center_moon = [r12 0];
moon_orb = 382500;
moon_rad = 1740;
earth_rad = 6378;

r0 = [1/2;sqrt(3)/2;0]*r12;
intercept = r0 + [384.388174;0;0];
r_int = norm(intercept);

G = 6.674 * 10^(-20); % km^3/kg-s^2
m1 = 5.97219 * 10^(24); % kg
mu = G * m1;

rp = 6678; % LEO
ra = r_int;

a = (rp+ra)/2;
e = 1 - rp/a;
p = a*(1 - e^2);
P = 2*pi*sqrt(a^3/mu);
P_half = P/2;
P_moon = 27.32*24*60^2;
n_frac_rad = 2*pi*P_half/P_moon;  % rad
n_frac_deg = n_frac_rad*180/pi;   % degrees

t = [];
x_sc = [];
y_sc = [];
x_moon = [];
y_moon = [];
f_sc = 0:0.001:pi;
f_moon = pi-n_frac_rad-pi/3:0.0003622:pi-pi/3;
for i = 1:length(f_sc)
    r_sc = p/(1+e*cos(f_sc(i)));
    x_sc(i) = r_sc*cos(f_sc(i));
    y_sc(i) = r_sc*sin(f_sc(i));
end
for i = 1:length(f_moon)
    x_moon(i) = r12*cos(f_moon(i));
    y_moon(i) = r12*sin(f_moon(i));
end

x_sc_full = [];
y_sc_full = [];
x_moon_full = [];
y_moon_full = [];
f = 0:0.0001:2*pi;
for i = 1:length(f)
    r_sc = p/(1+e*cos(f(i)));
    x_sc_full(i) = r_sc*cos(f(i));
    y_sc_full(i) = r_sc*sin(f(i));

    x_moon_full(i) = r12*cos(f(i));
    y_moon_full(i) = r12*sin(f(i));
end


%% Simulation

tspan = [0 105];
% x0 = [0.492506;0.85;0;0;0;0];
% x0 = [-5;5;0;-2;2;0];
x0 = [0;0.00008;0;0;0;0];
[~,dr] = ode89(@p9,tspan,x0);

% figure(1)
% plot(dr(:,1),dr(:,2))
% xlabel('x-Position (non-dim)')
% ylabel('y-Position (non-dim)')
% 
% figure(2)
% hold on
% plot(out_long(1,:),out_long(2,:))
% plot(out_short(1,:),out_short(2,:))
% legend('Long Period (s_1)','Short Period (s_2)')
% xlim([-900 900])
% ylim([-700 700])
% xlabel('x-Position (km)')
% ylabel('y-Position (km)')
% 
% figure(3)
% hold on
% plot(dr(:,1),dr(:,2),'r--')
% plot(out_short(1,:)/r12,out_short(2,:)/r12,'m')
% plot(out_long(1,:)/r12,out_long(2,:)/r12,'b')
% xlim([-900 900]/r12)
% ylim([-700 700]/r12)

figure(4)
contourf(X,Y,2*C-Z,[C C])
xlim([0.475 0.5])
ylim([0.855 0.875])
xlabel('x-Position')
ylabel('y-Position')
title('Zero Velocity Surface for L4')

% figure(5)
% hold on
% viscircles(center_earth,SOI_earth,'Color','k','LineStyle',':','LineWidth',0.5);
% viscircles(center_moon,SOI_moon,'Color','k','LineStyle',':','LineWidth',0.5);
% viscircles(center_earth,moon_orb,'Color','k','LineStyle','--','LineWidth',0.75);
% viscircles(center_earth,earth_rad,'Color','k');
% viscircles(center_moon,moon_rad,'Color','k');
% xlim([-1.2e6 1.2e6])
% ylim([-1.2e6 1.2e6])

figure(6)
hold on
plot(y_sc,-x_sc,'Color','k','LineWidth',2)
plot(y_moon,-x_moon,'Color','k','LineWidth',2)
plot(y_sc_full,-x_sc_full,'k--')
plot(y_moon_full,-x_moon_full,'k--')
viscircles([y_moon(end) -x_moon(end)],SOI_moon,'Color','k','LineStyle',':','LineWidth',0.5);
xlim([-4.5e5 4.5e5])
ylim([-4.5e5 4.5e5])
xlabel('x-Position (km)')
ylabel('y-Position (km)')



%% Functions

function dyn = dynICs_short(zeta0,eta0)

m2 = 7.347673 * 10^(22);    % Mass of Moon (kg)
m1 = 5.97219 * 10^(24);     % Mass of Earth (kg)
mu = m2 / (m2 + m1);        % Earth-Moon Mass Ratio

% s1 = 0.29820789544347032;
s2 = 0.95450094347526815;

Uxx = 3/4;
Uxy = 3*sqrt(3)/2 * (mu-1/2);

Gamma = (s2^2 + Uxx) / (4*s2^2 + Uxy^2);

zeta0dot = 1/2 * (Uxy*zeta0 + eta0/Gamma);
eta0dot = -1/2 * ((s2^2 + Uxx)*zeta0 + Uxy*eta0);

dyn = [zeta0dot;eta0dot];

end

function dyn = dynICs_long(zeta0,eta0)

m2 = 7.347673 * 10^(22);    % Mass of Moon (kg)
m1 = 5.97219 * 10^(24);     % Mass of Earth (kg)
mu = m2 / (m2 + m1);        % Earth-Moon Mass Ratio

s1 = 0.29820789544347032;
% s2 = 0.95450094347526815;

Uxx = 3/4;
Uxy = 3*sqrt(3)/2 * (mu-1/2);

Gamma = (s1^2 + Uxx) / (4*s1^2 + Uxy^2);

zeta0dot = 1/2 * (Uxy*zeta0 + eta0/Gamma);
eta0dot = -1/2 * ((s1^2 + Uxx)*zeta0 + Uxy*eta0);

dyn = [zeta0dot;eta0dot];

end

function state = short(tau,in)

% s1 = 0.29820789544347032;
s2 = 0.95450094347526815;

zeta0 = in(1);
eta0 = in(2);
zeta0dot = in(3);
eta0dot = in(4);

zeta = [];
eta = [];
for i = 1:length(tau)
    zeta(i) = zeta0*cos(s2*tau(i)) + zeta0dot/s2*sin(s2*tau(i));
    eta(i) = eta0*cos(s2*tau(i)) + eta0dot/s2*sin(s2*tau(i));
end

state = [zeta;eta];

end

function state = long(tau,in)

s1 = 0.29820789544347032;
% s2 = 0.95450094347526815;

zeta0 = in(1);
eta0 = in(2);
zeta0dot = in(3);
eta0dot = in(4);

zeta = [];
eta = [];
for i = 1:length(tau)
    zeta(i) = zeta0*cos(s1*tau(i)) + zeta0dot/s1*sin(s1*tau(i));
    eta(i) = eta0*cos(s1*tau(i)) + eta0dot/s1*sin(s1*tau(i));
end

state = [zeta;eta];

end

function dr_dt = p9(~,in)

m2 = 7.347673 * 10^(22);    % Mass of Moon (kg)
m1 = 5.97219 * 10^(24);     % Mass of Earth (kg)
mu = m2 / (m2 + m1);        % Earth-Moon Mass Ratio

dr = [in(1);in(2);in(3)];
dr_dot = [in(4);in(5);in(6)];

Omega = [0;0;1];
Otilde = [0 -Omega(3) Omega(2);Omega(3) 0 -Omega(1);-Omega(2) Omega(1) 0];

F = (1-mu)/4 * [-1 3*sqrt(3) 0;3*sqrt(3) 5 0;0 0 -4] + ...
    mu/4 * [-1 -3*sqrt(3) 0;-3*sqrt(3) 5 0;0 0 -4];

M = [zeros(3) eye(3);F-Otilde^2 -2*Otilde];
dr_dt = M * [dr;dr_dot];

end
