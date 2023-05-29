clear all
close all
clc

%% Defaults
FS = 18;
LW = 1.5;
set(0, 'defaultAxesFontSize', FS)
set(0, 'defaultLineLineWidth', LW)
set(0, 'defaultAxesLineWidth', 1.25)
set(0, 'defaultAxesTickLabelInterpreter', 'Latex')
set(0, 'defaultTextInterpreter', 'Latex')
set(0, 'defaultLegendInterpreter', 'Latex')

%% Data
D = 1;
d = 0.05;
rho_h2o = 1e3;
rho_atm = 1.225;
patm = 1e5;
p0 = patm;
g = 9.85;
h = 0.25;
L = 20;
T0 = 300;
A = 0.25*pi*D^2;
m = 1;
gamma = 1.4;

%% Initial conditions and parameters
xi0 = p0/patm;
eta0 = 1;
rho0 = rho_atm;
Pi0 = rho0/rho_atm;

AR = L/D;

Ca = rho_h2o.*A*L./m;
Co = rho_atm*A*L/m;
Cp = patm*A/(m*g);
Lambda = h/L;
delta = d/D;

tc = L/delta*sqrt(rho_atm/patm/gamma)*((gamma + 1)/2)^((gamma + 1)/...
    (2*gamma - 2));

alpha = L/(g*tc^2);
dc = 1000;
beta = dc*L/(m*g*tc);

a = sqrt(2/(gamma - 1))*((gamma + 1)/2)^((gamma + 1)/(2*gamma - 2));
b = ((gamma + 1)/2)^(gamma/(gamma - 1));

K = 1/gamma*(1 - Cp/Ca + 1/Ca);

%% Calculations (variable velocity)
%         tf = 10; %0.9*tc;
tspan = linspace(0, 1.5, 2*2048);
y0 = [1, 0, 1];
Opt = odeset('Events', @myEvent);
[t, y] = ode45(@(t, y) ode_mass(t, y, Ca, Cp, gamma, Lambda, beta,...
    alpha), tspan, y0, Opt);

tf = t(end);

eta = real(y(:,1));
etaDot = real(y(:,2));
xi = real(y(:,3));

zm = eta*L;
vm = etaDot*L/tc;
pg = xi*patm;

%% Characterization of the oscillation amplitude
tdd = 5*d;
[yupper, ylower] = envelope(eta(t <= tdd), 10, 'peak');
etaMean = 0.5*(yupper + ylower);
etaMax = max(yupper - etaMean);

%% Subsonic theoretical approximate solution (alpha, beta << 1)
% Coefficients
A = Cp/Ca*(1 + 1/gamma);
B = 1/gamma*(1 + Cp/Ca + 1/Ca);
u0 = sqrt(b^((gamma - 1)/gamma) - 1);

u_1 = 2.^(2/3).*(A+(-1).*B).*((-1)+gamma).*gamma.*((-3).*a.*(2.*A+B.*(( ...
    -2)+gamma)).^2.*((-1)+gamma).^2.*gamma.^2.*t+((2.*A+B.*((-2)+ ...
    gamma)).^3.*((-1)+gamma).^3.*gamma.^4.*(16.*(A+(-1).*B).^3.* ...
    gamma.^2+9.*a.^2.*(2.*A+B.*((-2)+gamma)).*((-1)+gamma).*t.^2)).^( ...
    1/2)).^(-1/3)+(-1).*2.^(-2/3).*(2.*A+B.*((-2)+gamma)).^(-1).* ...
    gamma.^(-1).*((-3).*a.*(2.*A+B.*((-2)+gamma)).^2.*((-1)+gamma) ...
    .^2.*gamma.^2.*t+((2.*A+B.*((-2)+gamma)).^3.*((-1)+gamma).^3.* ...
    gamma.^4.*(16.*(A+(-1).*B).^3.*gamma.^2+9.*a.^2.*(2.*A+B.*((-2)+ ...
    gamma)).*((-1)+gamma).*t.^2)).^(1/2)).^(1/3);

if sum(imag(u_1)) ~= 0
    u_1 = (1/2).*2.^(-2/3).*gamma.^(-1).*((-3).*a.*(2.*A+B.*((-2)+gamma)) ...
        .^2.*((-1)+gamma).^2.*gamma.^2.*t+((2.*A+B.*((-2)+gamma)).^3.*(( ...
        -1)+gamma).^3.*gamma.^4.*(16.*(A+(-1).*B).^3.*gamma.^2+9.*a.^2.*( ...
        2.*A+B.*((-2)+gamma)).*((-1)+gamma).*t.^2)).^(1/2)).^(-1/3).*((-4) ...
        .*(-2).^(1/3).*(A+(-1).*B).*((-1)+gamma).*gamma.^2+(1+(sqrt(-1)*( ...
        -1)).*3.^(1/2)).*(2.*A+B.*((-2)+gamma)).^(-1).*((-3).*a.*(2.*A+B.* ...
        ((-2)+gamma)).^2.*((-1)+gamma).^2.*gamma.^2.*t+((2.*A+B.*((-2)+ ...
        gamma)).^3.*((-1)+gamma).^3.*gamma.^4.*(16.*(A+(-1).*B).^3.* ...
        gamma.^2+9.*a.^2.*(2.*A+B.*((-2)+gamma)).*((-1)+gamma).*t.^2)).^( ...
        1/2)).^(2/3));
end

u_1 = real(u_1);

% In terms of eta
% XX = Ca/Cp;
% A1 = 1+(-1).*Cp.^(-1).*(1+Cp.^(-1).*(1+(-1).*gamma).*gamma.^(-1)).* ...
%     gamma.^(-1)+(-1).*gamma.^(-1).*XX; %u
% B1 = (1/3).*(gamma.^(-1)+(1+Cp.^(-1).*(1+(-1).*gamma).*gamma.^(-1)).* ...
%     gamma.^(-1)+((-1)+gamma).*gamma.^(-2).*XX); %u^3
% 
% K1 = a*XX/2*sqrt(gamma-1)/sqrt(gamma)*t;
% uu_1 = 6.^(-2/3).*B1.^(-1).*(9.*A1.*B1.^2.*(Cp.^(-1)).^(1/2)+9.*B1.^3.*( ...
%     Cp.^(-1)).^(3/2)+9.*B1.^2.*K1+3.^(1/2).*(B1.^3.*(4.*A1.^3+27.*B1.* ...
%     Cp.^(-3).*(B1+A1.*Cp+(Cp.^(-1)).^(-3/2).*K1).^2)).^(1/2)).^(-1/3).* ...
%     ((-2).*3.^(1/3).*A1.*B1+2.^(1/3).*(9.*A1.*B1.^2.*(Cp.^(-1)).^(1/2) ...
%     +9.*B1.^3.*(Cp.^(-1)).^(3/2)+9.*B1.^2.*K1+3.^(1/2).*(B1.^3.*(4.* ...
%     A1.^3+27.*B1.*Cp.^(-3).*(B1+A1.*Cp+(Cp.^(-1)).^(-3/2).*K1).^2)).^( ...
%     1/2)).^(2/3));

% Alternative perturbation solution with eta = 1 - phi, phi << 1
% phii_1 = (uu_1.^2 - 1/Cp)/XX;
%
xi_1 = 1 + gamma/(gamma - 1)*u_1.^2;
eta_1 = computeEta(xi_1, Ca, Cp); %1 - (Cp*(xi_1 - 1) - 1)/Ca;

% Critical time (xi = b)
t1 = t(xi_1 < b);
t2 = t(xi_1 >= b);
if isempty(t1)
    tb = 0;
else
    tb = t1(end);
end
% tbb = -A*2*gamma/(gamma + 1)*(1 - b^((1 + gamma)/(2*gamma))) +...
%     B*2*gamma/(1 - gamma)*(1 - b^((1 - gamma)/(2*gamma)));
% etab = computeEta(b, Ca, Cp);

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
eta_2 = computeEta(xi_1, Ca, Cp); %1 - (Cp*(xi_2 - 1) - 1)/Ca;

tf_ = 1 - (Ca/Cp + 1 + 1/Cp - b)*D_;

%% Critical time
tbb = psi*(2*b - 1 - b^2) + phi*(b - 1);
xibb = b;
etabb = computeEta(xibb, Ca, Cp);

%% Final values
% eta_f = 0;
% xi_f = Ca/Cp + 1 + 1/Cp;

%% Complete theoretical solution (alpha, beta << 1)
xi_ = [xi_1(xi_1 <= b); xi_2(xi_2 > b)];
eta_ = [eta_1(xi_1 <= b); eta_2(xi_2 > b)];

t_ = [t(xi_1 <= b); t(xi_2 > b)];

% figure, plot(t_, xi_, tb, b, 'o')
% figure, plot(t_, eta_, tb, etab, 'o')

%% Resultant force and system developed energy

F = -beta*etaDot - 1 + Cp*(xi - 1) - Ca*(1 - eta);
WtF = cumtrapz(t, F.*etaDot);
% WtF = lowpass(WtF, 1e-5, t(2)-t(1), 'ImpulseResponse', 'iir');
WtFF = WtF(end);

PtF = F.*etaDot;
% PtF = lowpass(PtF, 1e-5, t(2)-t(1), 'ImpulseResponse', 'iir');
PtFF = max(PtF);

% F_ = -1 + Cp*(xi_ - 1) - Ca*(1 - eta_);
% WtF_ = cumtrapz(t, F_.*gradient(eta_)./gradient(t));

%% Energy
% Energy obtained through a turbine expansion
Wt = gradient((eta - Lambda).*xi.^(1/gamma))./gradient(t)*...
    gamma/(gamma - 1).*(1 - xi.^((gamma - 1)/gamma));

Wt0 = lowpass(Wt, 1e-5, t(2)-t(1), 'ImpulseResponse', 'iir');

Wt_ = gradient((eta_).*xi_.^(1/gamma))./gradient(t_)*...
    gamma/(gamma - 1).*(1 - xi_.^((gamma - 1)/gamma));

%% Numerically, solving the nonlinear equation
% ff = @(t) 2.^(2/3).*(A+(-1).*B).*((-1)+gamma).*gamma.*((-3).*a.*(2.*A+B.*(( ...
%     -2)+gamma)).^2.*((-1)+gamma).^2.*gamma.^2.*t+((2.*A+B.*((-2)+ ...
%     gamma)).^3.*((-1)+gamma).^3.*gamma.^4.*(16.*(A+(-1).*B).^3.* ...
%     gamma.^2+9.*a.^2.*(2.*A+B.*((-2)+gamma)).*((-1)+gamma).*t.^2)).^( ...
%     1/2)).^(-1/3)+(-1).*2.^(-2/3).*(2.*A+B.*((-2)+gamma)).^(-1).* ...
%     gamma.^(-1).*((-3).*a.*(2.*A+B.*((-2)+gamma)).^2.*((-1)+gamma) ...
%     .^2.*gamma.^2.*t+((2.*A+B.*((-2)+gamma)).^3.*((-1)+gamma).^3.* ...
%     gamma.^4.*(16.*(A+(-1).*B).^3.*gamma.^2+9.*a.^2.*(2.*A+B.*((-2)+ ...
%     gamma)).*((-1)+gamma).*t.^2)).^(1/2)).^(1/3);
% 
% gg = @(t) (1 + gamma/(gamma-1)*ff(t).^2 - b);
% t0 = 0.1; % starting point
% % tb_ = fzero(gg, t0);
lambda = 1e-4;
myfun = @(eta) zeroGas(Cp, Ca, lambda, gamma, eta);
eta_ss = fzero(myfun, 0.5);
xi_ss = ((1 - lambda)./(eta_ss - lambda)).^gamma;


%% Dimensionless plots
figure,
subplot(311), plot(t, eta, t_, eta_, 'r--', ...
    tbb, etabb, 'ks'), hold on
ylim([0 1.25]), %xlim([0 1.0])
% plot(t(end), etaf, 'ro'), hold on
% plot(t(end), etaft, 'gs')
legend('RK45', 'Theo.', 'Theo.', 'location', 'best')
ylabel('$\eta$')

subplot(312), plot(t, etaDot, t_, gradient(eta_)./gradient(t_), 'r--')
%xlim([0 1.0])
ylabel('$\dot{\eta}$')

% subplot(313), plot(t, etaDdot)
% ylabel('$\ddot{\eta}$')
% xlabel('$\tau$')
% saveas(gcf, 'eta', 'fig');
% saveas(gcf, 'eta','epsc')
% saveas(gcf, 'eta','png')

% figure,
subplot(313), plot(t, xi, t_, xi_, 'r--', t, xi*0 + b, 'k--', ...
    tbb, b, 'ks'), hold on
%xlim([0 1.0])
% plot(t(end), xif, 'ro'), hold on
% plot(t(end), xift, 'gs')
% legend('Exact', 'Num. std.', 'Theo. std.')
ylabel('$\xi$')

% subplot(212), plot(t, xiDot)
% ylabel('$\dot{\xi}$')
xlabel('$\tau$')
%
print(gcf, 'eta_etaDot_xi_solution','-dpdf','-fillpage')

% saveas(gcf, ['figs/lambda_' num2str(lambda)], 'fig');
% saveas(gcf, 'xi','epsc')
% saveas(gcf, 'xi','png')

%% Energy plots
%Potencia y energía desarrollada por la masa cayendo
WtDot = gradient(Wt)./gradient(t);
WtDot_ = gradient(Wt_)./gradient(t_);

figure,
% subplot(211), plot(t(1:end-1)*tc, 3.6*WtDot(1:end-1))
% ylabel('$P_t$ (kW)')

% subplot(212),
plot(t*tc, 3.6*Wt)
ylabel('$W_t$ (kWh)')
xlabel('$t$ (s)')

% saveas(gcf, ['figs/dim_Wt_lambda_' num2str(lambda)], 'fig');

%
figure,
% subplot(211), plot(t(1:end-1), WtDot(1:end-1), t_, WtDot_, 'r--')
% ylabel('$\dot{\mathcal{W}}_t$')

% subplot(212),
plot(t, Wt, t_, Wt_, 'r--')
ylabel('$\mathcal{W}_t$')
xlabel('$\tau$')

% saveas(gcf, ['figs/Wt_lambda_' num2str(lambda)], 'fig');

figure,
subplot(211)
plot(t, WtF)
ylabel('$\mathcal{W}_m$')
% xlabel('$\tau/\tau_f$')
% axis square

subplot(212)
plot(t, PtF)
xlabel('$\tau/\tau_f$')
ylabel('$\mathcal{P}_m$')
% axis square


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

%%
figure,
subplot(211), plot(t, eta, t_, eta_, 'r--',...
    tbb, etabb, 'ks'), hold on
% ylim([0 1])
% plot(t(end), etaf, 'ro'), hold on
% plot(t(end), etaft, 'gs')
legend('RK45', 'Theo.', 'Theo.', 'location', 'southwest')
ylabel('$\eta$')

% subplot(312), plot(t, etaDot, t_, gradient(eta_)./gradient(t_), 'r--')
% ylabel('$\dot{\eta}$')

% subplot(313), plot(t, etaDdot)
% ylabel('$\ddot{\eta}$')
% xlabel('$\tau$')
% saveas(gcf, 'eta', 'fig');
% saveas(gcf, 'eta','epsc')
% saveas(gcf, 'eta','png')

% figure,
subplot(212), plot(t, xi, t_, xi_, 'r--',...
    t, xi*0 + b, 'k--', tbb, xibb, 'ks'), hold on
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

% save('subsonic.mat', 't', 'eta', 'etaDot', 'xi', 'b');


% figure,
% subplot(211), plot(t, eta, t, [ylower, yupper], 'k--',...
%     t, etaMean, 'r--'), hold on
% % ylim([0 1])
% % plot(t(end), etaf, 'ro'), hold on
% % plot(t(end), etaft, 'gs')
% legend('RK45', 'Env. up', 'Env. down', 'Mean', 'location', 'best')
% ylabel('$\eta$')
% 
% subplot(212), plot(t, etaMax),
% % plot(t(end), xif, 'ro'), hold on
% % plot(t(end), xift, 'gs')
% % legend('Exact', 'Num. std.', 'Theo. std.')
% ylabel('$\max{\eta}$')
% 
% % subplot(212), plot(t, xiDot)
% % ylabel('$\dot{\xi}$')
% xlabel('$\tau$')
% 

%% Phase portraits
figure('Renderer', 'painters', 'Position', [480 375 1024 600]),

subplot(131), plot(xi, eta, xi_, eta_, 'r--',...
    eta*0 + b, linspace(min(min(eta), min(eta_)),...
    max(max(eta), max(eta_)), length(eta)), 'k--')
xlabel('$\xi$')
ylabel('$\eta$')
% legend('RK45', 'Theo.', 'Theo.', 'location', 'southwest')
axis square

subplot(132), plot(etaDot, eta, gradient(eta_)./gradient(t_), eta_, 'r--')
xlabel('$\dot{\eta}$')
ylabel('$\eta$')
axis square

subplot(133), plot(xi, etaDot, xi_, gradient(eta_)./gradient(t_), 'r--',...
    etaDot*0 + b, linspace(min(min(etaDot), min(gradient(eta_)./gradient(t_))),...
    max(max(etaDot), max(gradient(eta_)./gradient(t_))), length(etaDot)), 'k--')
xlabel('$\xi$')
ylabel('$\dot{\eta}$')
axis square

print(gcf, 'phasePortraits','-dpdf','-fillpage')
print(gcf, 'phasePortraits','-dpdf','-fillpage')

% saveas(gcf, ['figs/eta_vs_xi_lambda_' num2str(lambda)], 'fig');

%%
% figure,
% subplot(211), plot(t, eta, t_, eta_, 'r--', tb, etab, 'ko',...
%     tbb, etabb, 'ks'), hold on
% % ylim([0 1])
% % plot(t(end), etaf, 'ro'), hold on
% % plot(t(end), etaft, 'gs')
% legend('RK45', 'Theo.', 'location', 'southwest')
% ylabel('$\eta$')
% xlabel('$\tau$')
%
% subplot(212), plot(xi, eta, xi_, eta_, 'r--', tb, b, 'ko',...
%     tbb, xibb, 'ks'), hold on
% % ylim([0 1])
% % plot(t(end), etaf, 'ro'), hold on
% % plot(t(end), etaft, 'gs')
% legend('RK45', 'Theo.', 'location', 'southwest')
% ylabel('$\eta$')
% xlabel('$\xi$')


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
%%
% fc = 30; Fs = 1/(t(2) - t(1));
% Wn = 2*(fc/Fs);
% 
% [b_,a_] = butter(8, Wn); %El filtrosss orden 8. LA fc = 5 Hz
% 
% xiMidline = filtfilt(b_, a_, xi); % La señal filtrá
% 
% % %         etaMax(i,j)
% figure, plot(t, xi, t, xiMidline, 'r')
