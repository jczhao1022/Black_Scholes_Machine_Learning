function [u,U] = bs_eur_put_RungeKutta(M,N,sigma,r,K,X,T)

dx = X / (M+1);  %  Space steps
dt = T / N;      %  Time steps
x = (1:M)'*dx;

a = -(sigma^2*(1:M)'.^2 - r*(1:M)') / 2;
b = sigma^2*(1:M)'.^2 + r;
c = -(sigma^2*(1:M)'.^2 + r*(1:M)') / 2;
a1 = [a(2:M);0];
c1 = [0;c(1:M-1)];
A = spdiags([a1 b c1],-1:1,M,M);

theta =  1 - 1/sqrt(2);
B = 1/dt*speye(M) + theta*A;
C = 1/dt*speye(M) - (1-theta)*A;
D = 1/dt*speye(M) - 1/2*A;

U = zeros(M,N);
g = max(K-x,0); 
u = g;
f = zeros(M,1);

for j = N:-1:1
    f(1) = -a(1)*K*exp(-r*(T-(j-1)*dt));
    rhs = C*u+f;
    u_temp = B \rhs;
    rhs = D*u-(1/2-theta)*A*u_temp+f;
    u = B \ rhs;
    U(:,j) = u;
end







