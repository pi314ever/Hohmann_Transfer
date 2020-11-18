function [value,isterminal,direction] = eventhit(t,x,xe,xm,Re,Rm,Rem)
% Stop if hit Moon or Earth
hitearth = sqrt((xe+x(1))^2+x(2)^2+x(3)^2)-Re;
hitmoon = sqrt((x(1)-xm)^2+x(2)^2+x(3)^2)-Rm;
outofbounds = (sqrt((xe+x(1))^2+x(2)^2+x(3)^2)<3*Rem);
value = hitmoon*hitearth*outofbounds;
isterminal = 1;
direction = 0;
end