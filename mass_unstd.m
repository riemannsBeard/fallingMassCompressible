function sol = mass_unstd(eta, etaDot, Mh2o, P0, Patm, Lambda, gamma, lambda, Fbrake)

% sol = -Mh2o.*(1 - eta) - 1 - Patm + P0.*(((1 - Lambda)./(eta(1) - Lambda)).^gamma - 1) -...
%     etaDot.*abs(etaDot)*Mh2o + etaDot*Mh2o*lambda + Fbrake;

sol = -Mh2o.*(1 - eta) - 1 + etaDot*Fbrake - ...
    Patm + P0.*(((1 - Lambda)./(eta - Lambda)).^gamma);

end

