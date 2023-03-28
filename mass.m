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
L = 12;
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
tf = 0.75*tc;
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

%% Subsonic theoretical approximate solution (alpha, beta << 1)
% Coefficients
A = Cp/Ca*(1 + 1/gamma);
B = 1/gamma*(1 + Cp/Ca + 1/Ca);
u0 = sqrt(b^((gamma - 1)/gamma) - 1);

% u_1 = -0.5*a*(gamma - 1)/gamma*t/(B - A);
u_1 = 2.^(2/3).*(A+(-1).*B).*((-1)+gamma).*gamma.*((-3).*a.*(2.*A+B.*(( ...
  -2)+gamma)).^2.*((-1)+gamma).^2.*gamma.^2.*t+((2.*A+B.*((-2)+ ...
  gamma)).^3.*((-1)+gamma).^3.*gamma.^4.*(16.*(A+(-1).*B).^3.* ...
  gamma.^2+9.*a.^2.*(2.*A+B.*((-2)+gamma)).*((-1)+gamma).*t.^2)).^( ...
  1/2)).^(-1/3)+(-1).*2.^(-2/3).*(2.*A+B.*((-2)+gamma)).^(-1).* ...
  gamma.^(-1).*((-3).*a.*(2.*A+B.*((-2)+gamma)).^2.*((-1)+gamma) ...
  .^2.*gamma.^2.*t+((2.*A+B.*((-2)+gamma)).^3.*((-1)+gamma).^3.* ...
  gamma.^4.*(16.*(A+(-1).*B).^3.*gamma.^2+9.*a.^2.*(2.*A+B.*((-2)+ ...
  gamma)).*((-1)+gamma).*t.^2)).^(1/2)).^(1/3);

% xi_1 = (u_1.^2 + 1).^(gamma/(gamma - 1));
xi_1 = 1 + gamma/(gamma-1)*u_1.^2;
eta_1 = 1 - (Cp*(xi_1 - 1) - 1)/Ca;

% Critical time (xi = b)
t1 = t(xi_1 < b);
tb = t1(end);
tbb = -A*2*gamma/(gamma + 1)*(1 - b^((1 + gamma)/(2*gamma))) +...
    B*2*gamma/(1 - gamma)*(1 - b^((1 - gamma)/(2*gamma)));
etab = 1 - (Cp*(b - 1) - 1)/Ca;

%% Numerically, solving the nonlinear equation
ff = @(t) 2.^(2/3).*(A+(-1).*B).*((-1)+gamma).*gamma.*((-3).*a.*(2.*A+B.*(( ...
  -2)+gamma)).^2.*((-1)+gamma).^2.*gamma.^2.*t+((2.*A+B.*((-2)+ ...
  gamma)).^3.*((-1)+gamma).^3.*gamma.^4.*(16.*(A+(-1).*B).^3.* ...
  gamma.^2+9.*a.^2.*(2.*A+B.*((-2)+gamma)).*((-1)+gamma).*t.^2)).^( ...
  1/2)).^(-1/3)+(-1).*2.^(-2/3).*(2.*A+B.*((-2)+gamma)).^(-1).* ...
  gamma.^(-1).*((-3).*a.*(2.*A+B.*((-2)+gamma)).^2.*((-1)+gamma) ...
  .^2.*gamma.^2.*t+((2.*A+B.*((-2)+gamma)).^3.*((-1)+gamma).^3.* ...
  gamma.^4.*(16.*(A+(-1).*B).^3.*gamma.^2+9.*a.^2.*(2.*A+B.*((-2)+ ...
  gamma)).*((-1)+gamma).*t.^2)).^(1/2)).^(1/3);

gg = @(t) (1 + gamma/(gamma-1)*ff(t).^2 - b);
t0 = 0.1; % starting point
tb_ = fzero(gg, t0);

%% Sonic theoretical approximate solution (alpha, beta << 1)
% Coefficients
C_ = A*2*gamma/(1 + gamma) - 2*gamma/(1 - gamma)*B;
D_ = A*b^((1 - gamma)/(2*gamma)) - B*b^((1 - 3*gamma)/(2*gamma));

% First order solution
% xi_2 = b + (t - tb)/D_;

% Second order solution
phi = A*b^((1-gamma)/(2*gamma)) - B*b^((1 - 3*gamma)/(2*gamma));
psi = A*(1-gamma)/(4*gamma)*b^((1-3*gamma)/(2*gamma)) -...
    B*(1 - 3*gamma)/(4*gamma)*b^((1 - 5*gamma)/(2*gamma));
xi_2 = b + (-phi + sqrt(phi^2 + 4*psi*(t - tb)))/(2*psi);


eta_2 = 1 - (Cp*(xi_2 - 1) - 1)/Ca;

tf = 1 - (Ca/Cp + 1 + 1/Cp - b)*D_;



%% Critical time
% tbb = (b - 1)*D_;
% xibb = 1 + tbb/D_;

tbb = psi*(2*b - 1 - b^2) + phi*(b - 1);
xibb = b;
etabb = 1 - (Cp*(xibb - 1) - 1)/Ca;


%% Final values
eta_f = 0;
xi_f = Ca/Cp + 1 + 1/Cp;

%% Complete theoretical solution (alpha, beta << 1)
xi_ = [xi_1(xi_1 < b); xi_2(xi_2 >= b)];
eta_ = [eta_1(xi_1 < b); eta_2(xi_2 >= b)];

t_ = [t(xi_1 < b); t(xi_2 >= b)];  

figure, plot(t_, xi_, tb, b, 'o')
figure, plot(t_, eta_, tb, etab, 'o')

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
subplot(311), plot(t, eta, t_, eta_, 'r--', tb, etab, 'o',...
    tbb, etab, 's'), hold on
% ylim([0 1])
% plot(t(end), etaf, 'ro'), hold on
% plot(t(end), etaft, 'gs')
legend('RK45', 'Theo.', 'location', 'southwest')
ylabel('$\eta$')

subplot(312), plot(t, etaDot, t, t*0 - 1, 'r--')
ylabel('$\dot{\eta}$')

% subplot(313), plot(t, etaDdot)
% ylabel('$\ddot{\eta}$')
% xlabel('$\tau$')
% saveas(gcf, 'eta', 'fig');
% saveas(gcf, 'eta','epsc')
% saveas(gcf, 'eta','png')

% figure,
subplot(313), plot(t, xi, t, xi*0 + b, 'k--', tb, b, 'o',...
    tbb, b, 's'), hold on
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
subplot(311), plot(t, eta, t_, eta_, 'r--', tb, etab, 'o',...
    tbb, etabb, 's'), hold on
% ylim([0 1])
% plot(t(end), etaf, 'ro'), hold on
% plot(t(end), etaft, 'gs')
legend('RK45', 'Theo.', 'location', 'southwest')
ylabel('$\eta$')

subplot(312), plot(t, etaDot, t, t*0 - 1, 'r--')
ylabel('$\dot{\eta}$')

% subplot(313), plot(t, etaDdot)
% ylabel('$\ddot{\eta}$')
% xlabel('$\tau$')
% saveas(gcf, 'eta', 'fig');
% saveas(gcf, 'eta','epsc')
% saveas(gcf, 'eta','png')

% figure,
subplot(313), plot(t, xi, t_, xi_, 'r--',...
    t, xi*0 + b, 'k--', tb, b, 'o', tbb, xibb, 's'), hold on
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

%% Phase portraits
figure,
subplot(221), plot(xi, eta, xi_, eta_, 'r--')
ylabel('$\eta$')
axis square

subplot(222), plot(etaDot, eta, gradient(eta_)./gradient(t_), eta_, 'r--')
xlabel('$\dot{\eta}$')
axis square

subplot(223), plot(xi, etaDot, xi_, gradient(eta_)./gradient(t_), 'r--')
xlabel('$\xi$')
ylabel('$\dot{\eta}$')
axis square

% saveas(gcf, ['figs/eta_vs_xi_lambda_' num2str(lambda)], 'fig');

%%
figure,
subplot(211), plot(t, eta, t_, eta_, 'r--', tb, etab, 'o',...
    tbb, etabb, 's'), hold on
% ylim([0 1])
% plot(t(end), etaf, 'ro'), hold on
% plot(t(end), etaft, 'gs')
legend('RK45', 'Theo.', 'location', 'southwest')
ylabel('$\eta$')
xlabel('$\tau$')

figure,
plot(xi, eta, xi_, eta_, 'r--', tb, b, 'o',...
    tbb, xibb, 's'), hold on
% ylim([0 1])
% plot(t(end), etaf, 'ro'), hold on
% plot(t(end), etaft, 'gs')
legend('RK45', 'Theo.', 'location', 'southwest')
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
