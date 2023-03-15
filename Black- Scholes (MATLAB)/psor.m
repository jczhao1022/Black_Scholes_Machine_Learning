function [x,it,norm_res] = psor(A,q,x,para_omega,tol,maxit)
%PSOR    Projected Successive Overrelaxation Method (PSOR) 
%   for linear complementarity problem (LCP)  
%   x >= 0, Ax+q >= 0, x'(Ax+q) = 0 
% 
%   [x,it,norm_res] = psor(A,q,x,para_omega,tol,maxit)
% 
%   Inputs:
%   ----------
%       A : m x m matrix   
%       q : given vector   
%       x : initial vector   
%       para_omega   : relaxation parameter (1.2)  
%       tol : tolerance   
%       maxit  : maximum iterations   
% 
%   Outputs:
%   ----------
%       x   : numerical solution 
%       it  : number of iterations   
%       norm_res   : norm of residual 
%
%
%   Copyright Dr. Ning Zheng 
%   September 1, 2022 


m = size(A,1);

r = -(A*x+q);
res = min(x,-r);
norm_res = zeros(maxit,1);
norm_res(1) = norm(res);
it = 0;
while norm_res(it+1)/norm_res(1) > tol   &&   it < maxit
    it = it + 1;
    for k = 1:m
        x(k) = x(k) + para_omega * (-q(k) - A(k, :) * x) / A(k,k); %SOR method
        x(k) = max(x(k),0);                               %Projection
    end
    r = -(A*x + q);
    res = min(x,-r);
    norm_res(it+1) = norm(res);
end
norm_res(it+2:end) = [];
