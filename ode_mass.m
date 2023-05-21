function dydt = ode_mass(t, y, Ca, Cp, gamma, Lambda, beta, alpha0)

% State variables
% y(1) -> eta
% y(2) -> eta'
% y(3) -> xi

dydt = zeros(3,1);

a = sqrt(2/(gamma - 1))*((gamma + 1)/2)^((gamma + 1)/(2*gamma - 2));
b = ((gamma + 1)/2)^(gamma/(gamma - 1));

dydt(1) = y(2);
dydt(2) = (-1 + Cp*(y(3) - 1) - Ca*(1 - y(1)) - beta*y(2))/alpha0;

if y(3) <= b
    dydt(3) = (-a*sqrt(y(3)^((gamma - 1)/gamma) - 1) - ...
        y(2)*y(3)^(1/gamma))/((y(1) - Lambda)/gamma*y(3)^(1/gamma - 1));
else
    dydt(3) = (-y(3)^((gamma + 1)/(2*gamma)) - y(2)*y(3)^(1/gamma))/...
        ((y(1) - Lambda)/gamma*y(3)^(1/gamma - 1));
    
end
           
end

