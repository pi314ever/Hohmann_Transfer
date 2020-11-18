function [value,isterminal,direction] = eventalpha(t,x)
% Stop at ideal burn location for transfer orbit
value = x(1);     %
isterminal = 1;   % stop the integration
direction = -1;   % negative direction


end