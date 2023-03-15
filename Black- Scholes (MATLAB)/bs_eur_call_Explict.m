function [u,U] = bs_eur_call_Explict(M,N,sigma,r,K,X,T)

dx = X / M;  %  Space steps
dt = T / N;      %  Time steps
x = (1:M)'*dx;

a = -(sigma^2*(1:M)'.^2 - r*(1:M)') / 2;
b = sigma^2*(1:M)'.^2 + r;
c = -(sigma^2*(1:M)'.^2 + r*(1:M)') / 2;


a(M) = a(M) + c(M);

b1 = 1/dt - b;
a1 = [a(2:M);0];
c1 = [0;c(1:M-1)];
B = spdiags([-a1 b1 -c1],-1:1,M,M);

U = zeros(M,N);

g = max(x-K,0); 
u = g;
f = zeros(M,1);
f(M) = -2*dx*c(M);
for j = N:-1:1    
    rhs = B*u +f;
    u = dt*rhs;
    U(:,j) = u;
end