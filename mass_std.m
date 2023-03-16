function F = mass_std(x, y, Ca, Cp, gamma, Lambda, beta, alpha0)


a = sqrt(2/(gamma - 1))*((gamma + 1)/2)^((gamma + 1)/(2*gamma - 2));
b = ((gamma + 1)/2)^(gamma/(gamma - 1));

eta = 

F = (-Ca*(1 - y(1)) - 1  + Cp*(y(3) - 1) - beta*y(2))/alpha0;

if y(3) <= b
    dydt(3) = (-a*sqrt(y(3)^((gamma - 1)/gamma) - 1) - ...
        y(2)*y(3)^(1/gamma))/((y(1) - Lambda)/gamma*y(3)^(1/gamma - 1));
else
    dydt(3) = (-y(3)^((gamma + 1)/(2*gamma)) - y(2)*y(3)^(1/gamma))/...
        ((y(1) - Lambda)/gamma*y(3)^(1/gamma - 1));
    
end


F(1) = -Mh2o*(1 - x(1)) - 1 + (x(2) - 1)*Pa;

end

