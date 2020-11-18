function [value,isterminal,direction] = eventmoon(t,x)
% Stop at shortest approach to Moon
value = x(1);     %
isterminal = 1;   % stop the integration
direction = -1;   % negative direction


end