function [sol] = ABM4(func,tspan,y0,h)
%ABM4 Adams Bashforth Method 4th Order
%   
x0 = tspan(1);
xf = tspan(end);
n = floor((xf-x0)/h);
x = zeros(1,n);
y = zeros(length(y0),n);

% Initial conditions
y(:,1) = y0;
x(1) = x0;

if n < 5
    % Use RK4 only
    sol = RK4(func,tspan,y0,h);
else
    % RK4 for first 4 steps
    tspanRK = [x0 x0+3*h];
    solRK = RK4(func,tspanRK,y0,h);
    x(1:4) = solRK.x;
    y(:,1:4) = solRK.y;
    for ii = 4:n-1
        x(ii+1) = x(ii) + h;
        ypred = y(:,ii)+h/24*(55*func(x(ii),y(:,ii))...
            -59*func(x(ii-1),y(:,ii-1))+37*func(x(ii-2),y(:,ii-2))...
            -9*func(x(ii-3),y(:,ii-3)));
        y(:,ii+1) = y(:,ii)+h/24*(9*func(x(ii+1),ypred)...
            +19*func(x(ii),y(:,ii))-5*func(x(ii-1),y(:,ii-1))...
            +func(x(ii-2),y(:,ii-2)));
    end
    sol.x = x;
    sol.y = y;
end

end