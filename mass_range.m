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
d = 0.002:0.001:0.05; %0.02;
rho_h2o = 1e3;
rho_atm = 1.225;
patm = 1e5;
p0 = patm;
g = 9.85;
h = 0.25;
L = 10:0.25:20; %12;
T0 = 300;
A = 0.25*pi*D^2;
m = 1;
gamma = 1.4;

for j = 1:length(d)
    for i = 1:length(L)
        %% Initial conditions and parameters
        xi0 = p0/patm;
        eta0 = 1;
        rho0 = rho_atm;
        Pi0 = rho0/rho_atm;
        
        AR = L(i)/D;
        
        Ca = rho_h2o.*A*L(i)./m;
        Co = rho_atm*A*L(i)/m;
        Cp = patm*A/(m*g);
        Lambda = h/L(i);
        delta = d(j)/D;
        
        tc(i,j) = L(i)/delta*sqrt(rho_atm/patm/gamma)*((gamma + 1)/2)^((gamma + 1)/...
            (2*gamma - 2));
        
        alpha = L(i)/(g*tc(i,j)^2);
        dc = 1;
        beta = dc*L(i)/(m*g*tc(i,j));
        
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
        
        tf(i,j) = t(end);
        
        eta = real(y(:,1));
        etaDot = real(y(:,2));
        xi = real(y(:,3));
        xiMean = 1/tf(i,j)*(cumtrapz(t, xi));
        I = find(xiMean >= b);
        if isempty(I)
            tb0(i,j) = NaN;
            etab0(i,j) = NaN; 
        else
            aux1 = t(I);
            aux2 = eta(I);
            tb0(i,j) = aux1(1);
            etab0(i,j) = aux2(1);           
        end
        
        zm = eta*L(i);
        vm = etaDot*L(i)/tc(i,j);
        pg = xi*patm;
        
        %% Characterization of the oscillation amplitude
        tdd = 5*d(j);
        [yupper, ylower] = envelope(eta(t <= tdd), 10, 'peak');
        etaMean = 0.5*(yupper + ylower);
        etaMax(i,j) = max(yupper - etaMean);
        
%         [yupper, ylower] = envelope(xi, 10, 'peak');
%         xiMean = 0.5*(yupper + ylower);
%         xiMax(i,j) = max(yupper - xiMean);        
        
%         etaMax(i,j)
%         plot(t, xi,...
%            t, [ylower, yupper], 'k--',...
%            t, xiMean, 'r')

        %% Characterization of the pressure amplitude
%         tdd = 5*d(j);

        fc = 64; Fs = 1/(t(2) - t(1));
        Wn = 2*(fc/Fs);

        [b_,a_] = butter(8, Wn); %El filtrosss orden 8. LA fc = 5 Hz

        xiMidline = filtfilt(b_, a_, xi); % La señal filtrá        plot(t, xi, t, xiMidline, 'r')

%         plot(t, xi, t, xiMidline, 'r')

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
        
%         if imag(u_1) ~= 0
%             u_1 = (-2).^(2/3).*(A+(-1).*B).*((-1)+gamma).*gamma.*((-3).*a.*(2.*A+B.* ...
%                 ((-2)+gamma)).^2.*((-1)+gamma).^2.*gamma.^2.*t+((2.*A+B.*((-2)+ ...
%                 gamma)).^3.*((-1)+gamma).^3.*gamma.^4.*(16.*(A+(-1).*B).^3.* ...
%                 gamma.^2+9.*a.^2.*(2.*A+B.*((-2)+gamma)).*((-1)+gamma).*t.^2)).^( ...
%                 1/2)).^(-1/3)+(1/2).*2.^(-2/3).*(1+sqrt(-1).*3.^(1/2)).*(2.*A+B.*( ...
%                 (-2)+gamma)).^(-1).*gamma.^(-1).*((-3).*a.*(2.*A+B.*((-2)+gamma)) ...
%                 .^2.*((-1)+gamma).^2.*gamma.^2.*t+((2.*A+B.*((-2)+gamma)).^3.*(( ...
%                 -1)+gamma).^3.*gamma.^4.*(16.*(A+(-1).*B).^3.*gamma.^2+9.*a.^2.*( ...
%                 2.*A+B.*((-2)+gamma)).*((-1)+gamma).*t.^2)).^(1/2)).^(1/3);
%         end
        
        
        % In terms of eta
%         XX = Ca/Cp;
%         A1 = 1+(-1).*Cp.^(-1).*(1+Cp.^(-1).*(1+(-1).*gamma).*gamma.^(-1)).* ...
%             gamma.^(-1)+(-1).*gamma.^(-1).*XX; %u
%         B1 = (1/3).*(gamma.^(-1)+(1+Cp.^(-1).*(1+(-1).*gamma).*gamma.^(-1)).* ...
%             gamma.^(-1)+((-1)+gamma).*gamma.^(-2).*XX); %u^3
%         
%         K1 = a*XX/2*sqrt(gamma-1)/sqrt(gamma)*t;
%         uu_1 = 6.^(-2/3).*B1.^(-1).*(9.*A1.*B1.^2.*(Cp.^(-1)).^(1/2)+9.*B1.^3.*( ...
%             Cp.^(-1)).^(3/2)+9.*B1.^2.*K1+3.^(1/2).*(B1.^3.*(4.*A1.^3+27.*B1.* ...
%             Cp.^(-3).*(B1+A1.*Cp+(Cp.^(-1)).^(-3/2).*K1).^2)).^(1/2)).^(-1/3).* ...
%             ((-2).*3.^(1/3).*A1.*B1+2.^(1/3).*(9.*A1.*B1.^2.*(Cp.^(-1)).^(1/2) ...
%             +9.*B1.^3.*(Cp.^(-1)).^(3/2)+9.*B1.^2.*K1+3.^(1/2).*(B1.^3.*(4.* ...
%             A1.^3+27.*B1.*Cp.^(-3).*(B1+A1.*Cp+(Cp.^(-1)).^(-3/2).*K1).^2)).^( ...
%             1/2)).^(2/3));
        
            
        % Alternative perturbation solution with eta = 1 - phi, phi << 1
        % phii_1 = (uu_1.^2 - 1/Cp)/XX;
        %
        xi_1 = 1 + gamma/(gamma-1)*u_1.^2;
        eta_1 = computeEta(xi_1, Ca, Cp); %1 - (Cp*(xi_1 - 1) - 1)/Ca;
        
        % Critical time (xi = b)
        t1 = t(xi_1 < b);
        t2 = t(xi_1 >= b);
        if isempty(t1)
            tb = 0;
        else
            tb = t1(end);
        end
%         tbb = -A*2*gamma/(gamma + 1)*(1 - b^((1 + gamma)/(2*gamma))) +...
%             B*2*gamma/(1 - gamma)*(1 - b^((1 - gamma)/(2*gamma)));
%         etab = 1 - (Cp*(b - 1) - 1)/Ca;
        
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
        etabb = 1 - (Cp*(xibb - 1) - 1)/Ca;              
               
        %% Final values
%         eta_f = 0;
%         xi_f = Ca/Cp + 1 + 1/Cp;
        
        %% Complete theoretical solution (alpha, beta << 1)
        xi_ = [xi_1(xi_1 < b); xi_2(xi_2 >= b)];
        eta_ = [eta_1(xi_1 < b); eta_2(xi_2 >= b)];
        
        t_ = [t(xi_1 < b); t(xi_2 >= b)];
        
        % figure, plot(t_, xi_, tb, b, 'o')
        % figure, plot(t_, eta_, tb, etab, 'o')
        
        %% Resultant force and system developed energy
        
        F = -beta*etaDot - 1 + Cp*(xi - 1) - Ca*(1 - eta);
        WtF = cumtrapz(t, F.*etaDot);
        WtF = lowpass(WtF, 1e-5, t(2)-t(1), 'ImpulseResponse', 'iir');
        WtFF(i,j) = WtF(end);
        
        %% Energy
        % Energy obtained through a turbine expansion
        Wt = gradient((eta - Lambda).*xi.^(1/gamma))./gradient(t)*...
            gamma/(gamma - 1).*(1 - xi.^((gamma - 1)/gamma));
        Wt_ = gradient((eta_ - Lambda).*xi_.^(1/gamma))./gradient(t_)*...
            gamma/(gamma - 1).*(1 - xi_.^((gamma - 1)/gamma));
        
        Wt0 = lowpass(Wt, 2e-6, t(2)-t(1), 'ImpulseResponse', 'iir');
        Wtt(i,j) = real(Wt0(end));
        Wtt_(i,j) = real(Wt_(end-3));
        
        %% Numerically, solving the nonlinear equation
        %         ff = @(t) 2.^(2/3).*(A+(-1).*B).*((-1)+gamma).*gamma.*((-3).*a.*(2.*A+B.*(( ...
        %             -2)+gamma)).^2.*((-1)+gamma).^2.*gamma.^2.*t+((2.*A+B.*((-2)+ ...
        %             gamma)).^3.*((-1)+gamma).^3.*gamma.^4.*(16.*(A+(-1).*B).^3.* ...
        %             gamma.^2+9.*a.^2.*(2.*A+B.*((-2)+gamma)).*((-1)+gamma).*t.^2)).^( ...
        %             1/2)).^(-1/3)+(-1).*2.^(-2/3).*(2.*A+B.*((-2)+gamma)).^(-1).* ...
        %             gamma.^(-1).*((-3).*a.*(2.*A+B.*((-2)+gamma)).^2.*((-1)+gamma) ...
        %             .^2.*gamma.^2.*t+((2.*A+B.*((-2)+gamma)).^3.*((-1)+gamma).^3.* ...
        %             gamma.^4.*(16.*(A+(-1).*B).^3.*gamma.^2+9.*a.^2.*(2.*A+B.*((-2)+ ...
        %             gamma)).*((-1)+gamma).*t.^2)).^(1/2)).^(1/3);
        %
        %         gg = @(t) (1 + gamma/(gamma-1)*ff(t).^2 - b);
        %         t0 = 0.1; % starting point
        %         % tb_ = fzero(gg, t0);

        if j == 1
            lambda = 1e-4;
            myfun = @(Eta) zeroGas(Cp, Ca, lambda, gamma, Eta);
            eta_ss(i) = fzero(myfun, 1 - 0.01*L(i));
            xi_ss(i) = ((1 - lambda)./(eta_ss(i) - lambda)).^gamma;
        end
        
        %% Storing and other stuff
        
        Ca_(i,j) = Ca;
        Cp_(i,j) = Cp;
        alpha_(i,j) = alpha;
        beta_(i,j) = beta;
        %     plot(t, xi, t, tf(i,j)*cumtrapz(t, xi))
        xiMax(i,j) = 1/tf(i,j)*(trapz(t, xi));
        xiMean_(i,j) = xiMean(end);
        
        I = t(xiMidline >= b);
        if isempty(I)
            tbNum(i,j) = NaN;
            tbb_(i,j) = NaN;
            etabNum(i,j) = NaN;
            etabb_(i,j) = NaN;
        else
            aux1 = t(xiMidline >= b);
            aux2 = eta(xiMidline >= b);
            tbNum(i,j) = aux1(1);
            tbb_(i,j) = tbb;
            etabNum(i,j) = aux2(1);
            etabb_(i,j) = etabb;
        end
        
    end
end

%% Calculate AR values for oscillations smaller than 10% and choked exit flow
LambdaC = repmat(L', 1, length(d));
% deltaC = repmat(d', 1, length(L));
I = etaMax < 0.1;
LambdaMax = LambdaC.*I;
LambdaMax = max(LambdaMax);

I = xiMax >= b;
I = double(I);
I(I == 0) = NaN;
LambdaXiMin = LambdaC.*I;
LambdaXiMin = min(LambdaXiMin);

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
subplot(313), plot(t, xi, t, xi*0 + b, 'k--', ...
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
figure,
subplot(211), plot(t*tc(i,j), 3.6*gradient(Wt)./gradient(t))
ylabel('$P_t$ (kW)')

subplot(212), plot(t*tc(i,j), 3.6*Wt)
ylabel('$W_t$ (kWh)')
xlabel('$t$ (s)')

% saveas(gcf, ['figs/dim_Wt_lambda_' num2str(lambda)], 'fig');

figure,
subplot(211), plot(t, gradient(Wt)./gradient(t),...
    t_, gradient(Wt_)./gradient(t_), 'r--')
ylabel('$\dot{\mathcal{W}}_t$')

subplot(212), plot(t, Wt, t_, Wt_, 'r--')
ylabel('$\mathcal{W}_t$')
xlabel('$\tau$')

% saveas(gcf, ['figs/Wt_lambda_' num2str(lambda)], 'fig');

%% Energy summary plot
figure,
contourf(d/D, L/D, Wtt, 12), hold on
xlabel('$\delta$'), ylabel('$\Lambda$')
cb = colorbar;
% caxis([1 4])
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\mathcal{W}_t$';
cb.Label.FontSize = 20;
cb.TickLabelInterpreter = 'latex';
plot(d/D, LambdaMax, 'r--', 'linewidth', 2.5)
plot(d/D, LambdaXiMin, 'k--', 'linewidth', 2.5)
axis square

figure,
contourf(d/D, L/D, WtFF, 12), hold on
xlabel('$\delta$'), ylabel('$\Lambda$')
cb = colorbar;
% caxis([1 4])
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\mathcal{W}_m$';
cb.Label.FontSize = 20;
cb.TickLabelInterpreter = 'latex';
plot(d/D, LambdaMax, 'r--', 'linewidth', 2.5)
plot(d/D, LambdaXiMin, 'k--', 'linewidth', 2.5)
axis square

figure,
p = pcolor(d/D, L/D, WtFF); hold on
shading flat
set(p,'EdgeColor','none');
xlabel('$\delta$'), ylabel('$\Lambda$')
cb = colorbar;
% caxis([1 4])
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\mathcal{W}_m$';
cb.Label.FontSize = 20;
cb.TickLabelInterpreter = 'latex';
plot(d/D, LambdaMax, 'r--', 'linewidth', 2.5)
plot(d/D, LambdaXiMin, 'k--', 'linewidth', 2.5)
axis square



% figure,
% subplot(121),
% contourf(d/D, L/D, Wtt, 12), hold on
% xlabel('$\delta$'), ylabel('$\Lambda$')
% % cb = colorbar;
% caxis([1 3])
% cb.Label.Interpreter = 'latex';
% % cb.Label.String = '$\mathcal{W}_t$';
% cb.Label.FontSize = 20;
% cb.TickLabelInterpreter = 'latex';
% plot(d/D, LambdaMax, 'r--', 'linewidth', 2.5)
% plot(d/D, LambdaXiMin, 'k--', 'linewidth', 2.5)
% axis square
% 
% subplot(122),
% contourf(d/D, L/D, Wtt_, 12), hold on
% xlabel('$\delta$'), %ylabel('$\Lambda$')
% yticks([]);
% cb = colorbar;
% caxis([1 3])
% cb.Label.Interpreter = 'latex';
% cb.Label.String = '$\mathcal{W}_t$';
% cb.Label.FontSize = 20;
% cb.TickLabelInterpreter = 'latex';
% plot(d/D, LambdaMax, 'r--', 'linewidth', 2.5)
% plot(d/D, LambdaXiMin, 'k--', 'linewidth', 2.5)
% axis square

figure,
subplot(121),
contourf(d/D, L/D, Ca_, 12)
xlabel('$\delta$'), ylabel('$\Lambda$')
cb = colorbar;
% caxis([1 4])
cb.Label.Interpreter = 'latex';
cb.Label.String = '$C_A$';
cb.Label.FontSize = 20;
cb.TickLabelInterpreter = 'latex';
axis square

subplot(122),
contourf(d/D, L/D, Cp_, 12)
xlabel('$\delta$'), %ylabel('$\Lambda$')
yticks([]);
cb = colorbar;
% caxis([1 4])
cb.Label.Interpreter = 'latex';
cb.Label.String = '$C_P$';
cb.Label.FontSize = 20;
cb.TickLabelInterpreter = 'latex';
axis square

figure,
subplot(121),
contourf(d/D, L/D, alpha_, 12)
xlabel('$\delta$'), ylabel('$\Lambda$')
cb = colorbar;
% caxis([1 4])
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\alpha$';
cb.Label.FontSize = 20;
cb.TickLabelInterpreter = 'latex';
axis square

subplot(122),
contourf(d/D, L/D, beta_, 12)
xlabel('$\delta$'), %ylabel('$\Lambda$')
yticks([]);
cb = colorbar;
% caxis([1 4])
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\beta$';
cb.Label.FontSize = 20;
cb.TickLabelInterpreter = 'latex';
axis square

figure,
contourf(d/D, L/D, Ca_./Cp_)
xlabel('$\delta$'),
ylabel('$\Lambda$')
cb = colorbar;
cb.Label.Interpreter = 'latex';
cb.Label.String = '$C_A/C_P$';
cb.Label.FontSize = 20;
cb.TickLabelInterpreter = 'latex';

figure,
contourf(d/D, L/D, xiMean_)
xlabel('$\delta$'),
ylabel('$\Lambda$')
cb = colorbar;
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\overline{\xi}$';
cb.Label.FontSize = 24;
cb.TickLabelInterpreter = 'latex';

figure,
contourf(d/D, L/D, tf)
xlabel('$\delta$'),
ylabel('$\Lambda$')
cb = colorbar;
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\tau_f$';
cb.Label.FontSize = 20;
cb.TickLabelInterpreter = 'latex';

figure,
contourf(d/D, L/D, etaMax), hold on
xlabel('$\delta$'),
ylabel('$\Lambda$')
cb = colorbar;
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\max{\eta}$';
cb.Label.FontSize = 20;
cb.TickLabelInterpreter = 'latex';
plot(d/D, LambdaMax, 'r--', 'linewidth', 2.5)

figure,
contourf(d/D, L/D, xiMax), hold on
xlabel('$\delta$'),
ylabel('$\Lambda$')
cb = colorbar;
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\max{\xi}$';
cb.Label.FontSize = 20;
cb.TickLabelInterpreter = 'latex';
plot(d/D, LambdaXiMin, 'r--', 'linewidth', 2.5)

figure,
subplot(221),
contourf(d/D, L/D, real(etab0))
xlabel('$\delta$'), ylabel('$\Lambda$')
cb = colorbar;
% caxis([0 0.75])
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\eta_b$';
cb.Label.FontSize = 20;
cb.TickLabelInterpreter = 'latex';
axis square

subplot(222),
contourf(d/D, L/D, tb0)
xlabel('$\delta$'), %ylabel('$\Lambda$')
cb = colorbar;
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\tau_b$';
cb.Label.FontSize = 20;
cb.TickLabelInterpreter = 'latex';
hold on,
% plot(d, d*0 + 10, 'r--', d, d*0 + 9, 'r--')
axis square

figure,
subplot(221),
contourf(d/D, L/D, etabNum)
xlabel('$\delta$'), ylabel('$\Lambda$')
cb = colorbar;
caxis([0.1 0.6])
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\eta_b$';
cb.Label.FontSize = 20;
cb.TickLabelInterpreter = 'latex';
axis square

subplot(222),
contourf(d/D, L/D, etabb_)
xlabel('$\delta$'), ylabel('$\Lambda$')
cb = colorbar;
caxis([0.1 0.6])
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\eta_{b0}$';
cb.Label.FontSize = 20;
cb.TickLabelInterpreter = 'latex';
axis square

subplot(223),
contourf(d/D, L/D, tbNum)
xlabel('$\delta$'), %ylabel('$\Lambda$')
cb = colorbar;
caxis([0.1 0.7])
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\tau_b$';
cb.Label.FontSize = 20;
cb.TickLabelInterpreter = 'latex';
axis square

subplot(224),
contourf(d/D, L/D, tbb_)
xlabel('$\delta$'), ylabel('$\Lambda$')
cb = colorbar;
caxis([0.1 0.7])
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\tau_{b0}$';
cb.Label.FontSize = 20;
cb.TickLabelInterpreter = 'latex';
axis square



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
subplot(311), plot(t, eta, t_, eta_, 'r--',...
    tbb, etabb, 'ks'), hold on
% ylim([0 1])
% plot(t(end), etaf, 'ro'), hold on
% plot(t(end), etaft, 'gs')
legend('RK45', 'Theo.', 'Theo.', 'location', 'southwest')
ylabel('$\eta$')

subplot(312), plot(t, etaDot, t_, gradient(eta_)./gradient(t_), 'r--')
ylabel('$\dot{\eta}$')

% subplot(313), plot(t, etaDdot)
% ylabel('$\ddot{\eta}$')
% xlabel('$\tau$')
% saveas(gcf, 'eta', 'fig');
% saveas(gcf, 'eta','epsc')
% saveas(gcf, 'eta','png')

% figure,
subplot(313), plot(t, xi, t_, xi_, 'r--',...
    t, xi*0 + b, 'k--', tb, b, 'ko', tbb, xibb, 'ks'), hold on
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
subplot(221), plot(xi, eta, xi_, eta_, 'r--',...
    eta*0 + b, linspace(min(min(eta), min(eta_)),...
    max(max(eta), max(eta_)), length(eta)), 'k--')
ylabel('$\eta$')
% legend('RK45', 'Theo.', 'Theo.', 'location', 'southwest')
axis square

subplot(222), plot(etaDot, eta, gradient(eta_)./gradient(t_), eta_, 'r--')
xlabel('$\dot{\eta}$')
axis square

subplot(223), plot(xi, etaDot, xi_, gradient(eta_)./gradient(t_), 'r--',...
    etaDot*0 + b, linspace(min(min(etaDot), min(gradient(eta_)./gradient(t_))),...
    max(max(etaDot), max(gradient(eta_)./gradient(t_))), length(etaDot)), 'k--')
xlabel('$\xi$')
ylabel('$\dot{\eta}$')
axis square

print(gcf, 'phasePortraits','-dpdf','-fillpage')

% saveas(gcf, ['figs/eta_vs_xi_lambda_' num2str(lambda)], 'fig');

%% Closed-valve solution
% Theoretical solution
[eta_ss0, xi_ss0] = computeTheoSS(Ca_(:,1), Cp_(:,1), gamma);

figure,
% scatter(L, eta_ss, [], xi_ss, 'filled')
yyaxis left
p1 = plot(L, eta_ss, '-', L, eta_ss0, '--');
ylabel('$\eta_s$')
yyaxis right
plot(L, xi_ss, '-', L, xi_ss0, '--')
ylabel('$\xi_s$')
xlabel('$\Lambda$')
axis square
legend([p1(1) p1(2)], 'Num.', 'Theo.', 'location', 'northwest')

% cb = colorbar;
% cb.TickLabelInterpreter = 'latex';
% cb.Label.String = '$\xi_ss$';
% cb.Label.Interpreter = 'latex';




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
