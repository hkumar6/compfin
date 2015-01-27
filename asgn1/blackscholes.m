format long; clear all; clc;

T=1; K=10; r=0.06; sig=0.3; d=0;

nx = 10; ny = 250;
U = zeros(nx, ny);

% initial conditions
for i=1:nx
    U(i,1) = max(1-2*i/nx, 0);
end

% boundary conditions
for i=1:ny
    U(1,i) = exp(-r*i/ny)*U(1,1);
    U(end, i) = exp(-d*i/ny)*U(end, 1);
end

k = T/ny; h = 1/nx;
% evaluating equation
for x=2:nx-1
    for y=2:ny
        U(x,y) = U(x+1,y-1)*(A(x/nx, k, h, sig)/2 + B(x/nx, r, d, k, h))...
            + U(x, y-1)*( A(x/nx, k, h, sig) + B(x/nx, r, d, k, h) + k*r*(1-x/nx) + d*k*x/nx + 1 )...
            + U(x-1, y-1)*A(x/nx, k, h, sig)/2;
    end
end

k/h^2
surf(U)