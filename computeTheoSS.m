function [eta_ss0, xi_ss0] = computeTheoSS(Ca, Cp, gamma)

eps = (1/2).*Cp.^(-1).*(Ca+(-1).*gamma.*(Cp-(gamma.^(-2).*(4.*Cp.*gamma+ ...
  (Ca+(-1).*Cp.*gamma).^2)).^(1/2)));

% aa = 0.5*Cp*gamma*(gamma - 1) + Ca;
% bb = gamma*Cp - Ca;
% 
% eps = (-bb + sqrt(bb.^2 + 4*aa))./(2*aa);


xi_ss0 = 1 + eps;
eta_ss0 = 1./(1 + 1/gamma*eps);


end

