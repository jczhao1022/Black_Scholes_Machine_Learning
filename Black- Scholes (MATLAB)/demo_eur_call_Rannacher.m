% function demo_eur_call_Rannacher
clear;clc

%% Option model parameters 
sigma = 0.6;       % volatility
r = 0.25;        % interest rate 
K = 10;          % strike price
X = 50;               % asset value
T = 1;                % time


%% European put option calculation 
% 验证空间方向二阶精度
disp('验证空间方向二阶精度')
M = 10*2.^(0:8)-1;
N = 40;

close all

tic;
[u_star, U_star]  = bs_eur_call_Rannacher(M(4),N,sigma,r,K,X,T);
% u_star  = bs_eur_call_Rannacher(M(end),N,sigma,r,K,X,T);
%  u_star = bs_eur_call_fd(M(end),N,sigma,r,K,X,T,'Rannacher');
%  u_star = bs_eur_fd(M(end),N,sigma,r,K,X,T,'call','Rannacher');
toc;
surf(U_star)
fprintf('Fine discretization grids (M,N)=(%.0f,%.0f)  \n',[M(end);N]);

%%
err = zeros(length(M),1);
ratio = zeros(length(M),1);
for k = 1:length(M)
%     u = bs_eur_fd(M(k),N,sigma,r,K,X,T,'call','Rannacher');
%     u = bs_eur_call_fd(M(k),N,sigma,r,K,X,T,'Rannacher');
    u   = bs_eur_call_Rannacher(M(k),N,sigma,r,K,X,T);
    u_temp = u_star(2^8/2^(k-1):2^8/2^(k-1):end);
    err(k) = norm(u-u_temp) / norm(u_temp);
    if k>1
        ratio(k) = log2( err(k-1) ./ err(k) );
    end
    fprintf('Rannacher: (M,N)=(%.0f,%.0f) Error=%.2e Order=%.2f\n',[M(k);N;err(k);ratio(k)]);
end

%% 验证时间方向二阶精度
disp('验证时间方向二阶精度')
M = 10;
N =100*2.^(0:8);

tic;
u_star  = bs_eur_call_Rannacher(M,N(end),sigma,r,K,X,T);
%  u_star = bs_eur_call_fd(M,N(end),sigma,r,K,X,T,'Rannacher');
%  u_star = bs_eur_fd(M,N(end),sigma,r,K,X,T,'call','Rannacher');
toc;

fprintf('Fine discretization grids (M,N)=(%.0f,%.0f)  \n',[M;N(end)]);

%%
err = zeros(length(N),1);
ratio = zeros(length(N),1);
for k = 1:length(N)
%     u = bs_eur_fd(M,N(k),sigma,r,K,X,T,'call','Rannacher');
%     u = bs_eur_call_fd(M,N(k),sigma,r,K,X,T,'Rannacher');
    u   = bs_eur_call_Rannacher(M,N(k),sigma,r,K,X,T);
    u_temp = u_star(2^8/2^(k-1):2^8/2^(k-1):end);
    err(k) = norm(u-u_temp) / norm(u_temp);
    if k>1
        ratio(k) = log2( err(k-1) ./ err(k) );
    end
    fprintf('Implicit Euler: (M,N)=(%.0f,%.0f) Error=%.2e Order=%.2f\n',[M;N(k);err(k);ratio(k)]);
end
