function [value,isterminal,direction] = eventalpha(t,x,al,xe,xm,Re,Rm)
% Stop at ideal burn location for transfer orbit
hitearth = sqrt((xe+x(1))^2+x(2)^2+x(3)^2)-Re;
hitmoon = sqrt((x(1)-xm)^2+x(2)^2+x(3)^2)-Rm;
value = (acos((xe+x(1))/sqrt((xe+x(1))^2+x(2)^2+x(3)^2))-al)*hitmoon*hitearth;
isterminal = 1;
direction = 0;
end