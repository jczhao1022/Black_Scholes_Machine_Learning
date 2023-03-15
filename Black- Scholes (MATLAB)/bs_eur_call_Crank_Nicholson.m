function [u,U] = bs_eur_call_Crank_Nicholson(M,N,sigma,r,K,X,T)

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
    b1 = 1/dt + b/2;
    % 1/tau - B/2
    b2 = 1/dt - b/2;
    % 1/tau + B/2
    
    a1 = [a(2:M);  0];
    c1 = [0;  c(1:M-1)];
    C = spdiags([a1  b1  c1],-1:1,M,M);
    % 1/tau - B/2
    B = spdiags([-a1  b2  -c1],-1:1,M,M);
    % 1/tau + B/2
    
    U = zeros(M,N);
    g = max(x-K,0); 
    u = g;
    f = zeros(M,1);
    % f(M) = -2*dx*c(M);
    f(M) = -2*c(M)*(X - K);
    for j = N:-1:1    
        rhs = B*u + f;
        u = C \ rhs;
        U(: , j) = u;
    end
    U = [zeros(1, N +1);
        U, g;
        (X-K)*ones(1, N+1)];
