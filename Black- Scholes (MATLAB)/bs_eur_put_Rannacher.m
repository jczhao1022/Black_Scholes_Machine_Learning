function [u,U] = bs_eur_put_Rannacher(M,N,sigma,r,K,X,T)

dx = X / (M+1);  %  Space steps
dt = T / N;      %  Time steps
x = (1:M)'*dx;

a = -(sigma^2*(1:M)'.^2 - r*(1:M)') / 2;
b = sigma^2*(1:M)'.^2 + r;
c = -(sigma^2*(1:M)'.^2 + r*(1:M)') / 2;

g = max(K-x,0); %put option payoff column vector
u = g;
f = zeros(M,1);

% 隐式欧拉格式
b = b + 1/dt;
a1 = [a(2:M);0];
c1 = [0;c(1:M-1)];
B = spdiags([a1 b c1],-1:1,M,M);
U = zeros(M,N);
for j = N:-1:N-3
    f(1) = -a(1)*K*exp(-r*(T-j*dt));
    rhs = u/dt +f;
    u = B \ rhs;
    U(:,j) = u;
end



% Crank-Nicolson 格式 
a = a/2;
c = c/2;
b1 = 1/dt + b/2;
b2 = 1/dt - b/2;

a1 = [a(2:M);  0];
c1 = [0;  c(1:M-1)];
C = spdiags([a1  b1  c1],-1:1,M,M);
B = spdiags([-a1  b2  -c1],-1:1,M,M);

for j = N-4:-1:1
    f(1) = -a(1)*K*( exp(-r*(T-(j-1)*dt)) + exp(r*(T-j*dt)) )/2;
    rhs = B*u + f;
    u = C \ rhs;
    U(: , j) = u;
end
