function [value, isterminal, direction] = myEvent(t, y)
value      = (y(1) - 0.125 < 0);
isterminal = 1;   % Stop the integration
direction  = 0;
end

