function [value,isterminal,direction] = eventmoon(t,x,xe,xm,Re,Rm)
% Stop at shortest approach to Moon
hitearth = sqrt((xe+x(1))^2+x(2)^2+x(3)^2)-Re;
hitmoon = sqrt((x(1)-xm)^2+x(2)^2+x(3)^2)-Rm;
value = ((x(1)-xm)*x(4)+x(2)*x(5)+x(3)*x(6))*hitmoon*hitearth+...
    0*(sqrt((x(1)-xm)^2+x(2)^2+x(3)^2)<(8*Rm));
isterminal = 1;
direction = 0;
end