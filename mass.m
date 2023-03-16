clear all
close all
clc

%% Defaults
FS = 18;
LW = 1.5;
set(0, 'defaultAxesFontSize', FS)
set(0, 'defaultLineLineWidth', LW)
set(0, 'defaultAxesTickLabelInterpreter', 'Latex')
set(0, 'defaultTextInterpreter', 'Latex')
set(0, 'defaultLegendInterpreter', 'Latex')

%% Data
D = 1;
d = 0.02;
rho_h2o = 1e3;
rho_atm = 1.225;
patm = 1e5;
p0 = patm;
g = 9.85;
h = 0.25;
L = 10;
T0 = 300;
A = 0.25*pi*D^2;
m = 100;
gamma = 1.4;
delta = d/D;

%% Initial conditions and parameters
xi0 = p0/patm;
eta0 = 1;
rho0 = rho_atm;
Pi0 = rho0/rho_atm;

Ca = rho_h2o.*A*L./m;
Cp = patm*A/(m*g);

Lambda = h/L;

tc = L/delta*sqrt(rho_atm/patm/gamma)*((gamma + 1)/2)^((gamma + 1)/...
    (2*gamma - 2));

alpha0 = L/(g*tc^2);
beta = alpha0;

a = sqrt(2/(gamma - 1))*((gamma + 1)/2)^((gamma + 1)/(2*gamma - 2));
b = ((gamma + 1)/2)^(gamma/(gamma - 1));

K = 1/gamma*(1 - Cp/Ca + 1/Ca);

%% Calculations (variable velocity)
tf = 0.7*tc;
tspan = linspace(0, tf/tc, 1e3);
y0 = [1, 0, 1];

[t, y] = ode23s(@(t, y) ode_mass(t, y, Ca, Cp, gamma, Lambda, beta,...
    alpha0), tspan, y0);

eta = y(:,1);
etaDot = y(:,2);
xi = y(:,3);

zm = eta*L;
vm = etaDot*L/tc;
pg = xi*patm;

%% Theoretical approximate solution (alpha, beta << 1)
A = 1/(Cp/Ca - 1/gamma*(1 + 1/Ca - Lambda));
% eta_ = 1 - (Cp*(xi_ - 1) - 1)/Ca;
eta_ = 1 - t;
xi_ = 1 - ((eta - 1)*Ca + 1)/Cp;


%% Energy
% Energy obtained through a turbine expansion
Wt = gradient((eta - Lambda).*xi.^(1/gamma))./gradient(t)*...
    gamma/(gamma - 1).*(1 - xi.^((gamma - 1)/gamma));


%% Dimensional plots
% figure,
% subplot(311), plot(t*tc, zm), hold on
% ylim([0 L])
% % plot(t(end)*tc, etaf*L, 'ro'), hold on
% % plot(t(end)*tc, etaft*L, 'gs')
% % legend('Exact', 'Num. std.', 'Theo. std.')
% ylabel('$z$ (m)')
% % 
% subplot(312), plot(t*tc, vm)
% ylabel('$\dot{z}$ (m/s)')
% % 
% % subplot(313), plot(t*tc, am)
% % ylabel('$\ddot{z}$ (m/s)')
% % xlabel('$t$ (s)')
% % saveas(gcf, 'z', 'fig');
% % saveas(gcf, 'z','epsc')
% % saveas(gcf, 'z','png')
% % 
% % figure,
% subplot(313), plot(t*tc, pg),hold on
% % plot(t(end)*tc, xif*patm, 'ro'),hold on
% % plot(t(end)*tc, xift*patm, 'gs')
% % legend('Exact', 'Num. std.', 'Theo. std.')
% ylabel('$p$ (Pa)')
% % 
% % subplot(212), plot(t*tc, pDot)
% % ylabel('$\dot{p}$ (Pa/s)')
% xlabel('$t$ (s)')
% % saveas(gcf, ['figs/dim_lambda_' num2str(lambda)], 'fig');
% % saveas(gcf, 'p','epsc')
% % saveas(gcf, 'p','png')

%% Dimensionless plots
figure,
subplot(311), plot(t, eta, t, eta_, 'r--'), hold on
% ylim([0 1])
% plot(t(end), etaf, 'ro'), hold on
% plot(t(end), etaft, 'gs')
legend('RK45', 'Theo. std.', 'location', 'southwest')
ylabel('$\eta$')

subplot(312), plot(t, etaDot)
ylabel('$\dot{\eta}$')

% subplot(313), plot(t, etaDdot)
% ylabel('$\ddot{\eta}$')
% xlabel('$\tau$')
% saveas(gcf, 'eta', 'fig');
% saveas(gcf, 'eta','epsc')
% saveas(gcf, 'eta','png')

% figure,
subplot(313), plot(t, xi, t, xi*0 + b, 'k--'), hold on
% plot(t(end), xif, 'ro'), hold on
% plot(t(end), xift, 'gs')
% legend('Exact', 'Num. std.', 'Theo. std.')
ylabel('$\xi$')

% subplot(212), plot(t, xiDot)
% ylabel('$\dot{\xi}$')
xlabel('$\tau$')
% 
% saveas(gcf, ['figs/lambda_' num2str(lambda)], 'fig');
% saveas(gcf, 'xi','epsc')
% saveas(gcf, 'xi','png')

%%
figure,
subplot(311), plot(t, eta, t, eta_, 'r--'), hold on
% ylim([0 1])
% plot(t(end), etaf, 'ro'), hold on
% plot(t(end), etaft, 'gs')
legend('RK45', 'Theo. std.', 'location', 'southwest')
ylabel('$\eta$')

subplot(312), plot(t, etaDot)
ylabel('$\dot{\eta}$')

% subplot(313), plot(t, etaDdot)
% ylabel('$\ddot{\eta}$')
% xlabel('$\tau$')
% saveas(gcf, 'eta', 'fig');
% saveas(gcf, 'eta','epsc')
% saveas(gcf, 'eta','png')

% figure,
subplot(313), plot(t, xi, t, xi_, 'r--', t, xi*0 + b, 'k--'), hold on
% plot(t(end), xif, 'ro'), hold on
% plot(t(end), xift, 'gs')
% legend('Exact', 'Num. std.', 'Theo. std.')
ylabel('$\xi$')

% subplot(212), plot(t, xiDot)
% ylabel('$\dot{\xi}$')
xlabel('$\tau$')
% 
% saveas(gcf, ['figs/lambda_' num2str(lambda)], 'fig');
% saveas(gcf, 'xi','epsc')
% saveas(gcf, 'xi','png')


figure,
subplot(221), plot(xi, eta)
ylabel('$\eta$')
axis square

subplot(222), plot(etaDot, eta)
xlabel('$\dot{\eta}$')
axis square

subplot(223), plot(xi, etaDot)
xlabel('$\xi$')
ylabel('$\dot{\eta}$')
axis square

% saveas(gcf, ['figs/eta_vs_xi_lambda_' num2str(lambda)], 'fig');

%%
figure,
subplot(211), plot(t, eta, t, eta_, 'r--'), hold on
% ylim([0 1])
% plot(t(end), etaf, 'ro'), hold on
% plot(t(end), etaft, 'gs')
legend('RK45', 'Theo. std.', 'location', 'southwest')
ylabel('$\eta$')
xlabel('$\tau$')

figure,
plot(xi, eta, xi, eta_, 'r--'), hold on
% ylim([0 1])
% plot(t(end), etaf, 'ro'), hold on
% plot(t(end), etaft, 'gs')
legend('RK45', 'Theo. std.', 'location', 'southwest')
ylabel('$\eta$')
xlabel('$\xi$')
axis square

%% Energy plots
%Potencia y energía desarrollada por la masa cayendo
figure,
subplot(211), plot(t*tc, 3.6*gradient(Wt)./gradient(t))
ylabel('$P_t$ (kW)')

subplot(212), plot(t*tc, 3.6*Wt)
ylabel('$W_t$ (kWh)')
% saveas(gcf, ['figs/dim_Wt_lambda_' num2str(lambda)], 'fig');

%

figure,
subplot(211), plot(t*tc, gradient(Wt)./gradient(t))
ylabel('$\dot{\mathcal{W}}_t$')

subplot(212), plot(t*tc, Wt)
ylabel('$\mathcal{W}_t$')
% saveas(gcf, ['figs/Wt_lambda_' num2str(lambda)], 'fig');


%% Temperature
% T = (p0./pg).^((1 - gamma)/gamma)*T0;
% theta = (T - T0)./(max(T) - T0);
% 
% figure,
% subplot(311), plot(t*tc, T)
% ylabel('$T$ (K)')
% subplot(312), plot(t*tc, rhog)
% ylabel('$\rho_g$ (kg/m$^3$)')
% subplot(313), plot(t*tc, G*G0),
% ylabel('$G$ (kg/s)')
% xlabel('$t$ (s)')
% 
% figure,
% subplot(311), plot(t, theta)
% ylabel('$\theta$')
% subplot(312), plot(t, rhog/rho0)
% ylabel('$\Pi$')
% subplot(313), plot(t, G),
% ylabel('$\mathcal{G}$')
% xlabel('$\tau$')
% 
% %%
% figure, plot(t*tc, A*(zm - h).*rhog, 'r--'), hold on
% plot(t*tc, cumtrapz(t*tc, G*G0)), hold on
% 
% 
