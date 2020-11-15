function [sol] = RK4(func,tspan,y0,h)
% func = function of x and y
x0 = tspan(1);
xf = tspan(end);
n = floor((xf-x0)/h)+1;
x = zeros(1,n);
y = zeros(length(y0),n);

% Initial conditions
y(:,1) = y0;
x(1) = x0;

for ii = 1:n-1
    x(ii+1) = x0+ii*h;
    y(:,ii+1) = RK4_calc(x(ii),h,y(:,ii));
end

% Hit final time value straight on
if xf > x(end)
    h_end = xf-x(n);
    x(n+1) = xf;
    y(:,n+1) = RK4_calc(x(n),h_end,y(:,n));
end

% Save to solution 
sol.x = x;
sol.y = y;

function yi2 = RK4_calc(xi,dx,yi1)
    k1 = func(xi,yi1);
    k2 = func(xi+.5*dx,yi1+.5*k1*dx);
    k3 = func(xi+.5*dx,yi1+.5*k2*dx);
    k4 = func(xi+dx,yi1+k3*dx);
    yi2 = yi1 + (1/6)*(k1+2*k2+2*k3+k4)*dx; 
end
end