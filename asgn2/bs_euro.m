format long; clc;

% call

% initialise variables
%  parameters
Smax = 75; Smin = 25; K = 50; sigma = 0.03; rate = 0.06; div = 0;
qd = 2*(rate - div)/(sigma^2); q = 2*rate/sigma^2;
%  time
T = 1; tau = T*sigma^2/2; n = 100;
%  x-axis
a = log(Smin/K); b = log(Smax/K); m = 50;
%  solution
lambda = (tau/n)/((b-a)/m)^2;
w = zeros(m,n);


% initial conditions
for i=1:m
    x = a + i*(b-a)/m;
    w(i,1) = max(exp(x*(qd+1)/2) - exp(x*(qd-1)/2), 0);
end

% boundary conditions
for i=1:n
    t = i*tau/n;
    w(1,i) = 0;
    w(end,i) = exp(b*(qd+1)/2 + t*(qd+1)^2/4);
end

% iterations
G = full(gallery('tridiag',m-2,-1,2,-1));
I = eye(m-2);
for i=2:n
    c = (I - (lambda/2)*G)*w(2:end-1, i-1) + (lambda/2)*[w(1,i-1)+w(1,i); zeros(m-4,1); w(end,i-1)+w(end,i)];
    w(2:end-1,i) = (I + (lambda/2)*G)\c;
end

% inverse transform to get value of option
for i=1:m % x
    for j=1:n % t
        x = a + i*(b-a)/m; t = j*tau/n;
        w(i,j) = K * exp(-(qd-1)*x/2 - ((qd-1)^2/4 + q)*t);
    end
end