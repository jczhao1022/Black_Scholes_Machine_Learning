function [u,it,U] = bs_am_call_implict(M,N,sigma,r,K,X,T,para_omg)

% dx = X / (M);  %  Space steps
% dt = T / N;      %  Time steps
% x = (1:M)'*dx;
%
% a = -(sigma^2*(1:M)'.^2 - r*(1:M)') / 2;
% b = sigma^2*(1:M)'.^2 + r;
% c = -(sigma^2*(1:M)'.^2 + r*(1:M)') / 2;
% a(M) = a(M) + c(M);
%
% b = b + 1/dt;
% a1 = [a(2:M);0];
% c1 = [0;c(1:M-1)];
% B = spdiags([a1 b c1],-1:1,M,M);

dx = X / (M+1);  %  Space steps
dt = T / (N+1);      %  Time steps
x = (1:M)'*dx;

a = -(sigma^2*(1:M)'.^2 - r*(1:M)') / 2;
b = sigma^2*(1:M)'.^2 + r;
c = -(sigma^2*(1:M)'.^2 + r*(1:M)') / 2;
% -B
% a(M) = a(M) + c(M);

a = a/2;
c = c/2;
% b1 = 1/dt + b/2;
% 1/tau - B/2
b2 = 1/dt - b/2;
% 1/tau + B/2

a1 = [a(2:M);  0];
c1 = [0;  c(1:M-1)];
% C = spdiags([a1  b1  c1],-1:1,M,M);
% 1/tau - B/2
B = spdiags([-a1  b2  -c1],-1:1,M,M);
% 1/tau + B/2

maxit = 100;
tol = 1e-6;


U = zeros(M,N);

g = max(x-K,0);
u = g;
v = 0*g;
f = zeros(M,1);
f(M) = -2*c(M)*(X - K);
it = zeros(N,1);
for j = N:-1:1
    rhs = u/dt +f;
    q = B*g - rhs;
    [v,it(j)] = psor(B,q,v,para_omg,tol,maxit);
    u = v + g;
    U(:,j) = u;
end
it = mean(it);
U = [zeros(1, N +1);
     U, g;
     (X-K)*ones(1, N+1)];