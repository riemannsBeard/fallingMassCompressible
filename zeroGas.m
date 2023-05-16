function y = zeroGas(Cp, Ca, lambda, gamma, eta)

y = Cp*(((1 - lambda)./(eta - lambda)).^gamma - 1) - Ca*(1 - eta) - 1;

end

