function demo_eur_RungeKutta
%% 使用RungeKutta格式计算欧式期权

clear;clc

%% 模型参数
sigma = 0.6;       % volatility
r = 0.25;        % interest rate 
K = 10;          % strike price
X = 50;               % asset value
T = 1;                % time


%% 欧式看跌期权计算 
% 验证空间方向二阶精度
disp('验证空间方向二阶精度')
M = 10*2.^(0:8)-1;
N = 40000;

tic;
u_star  = bs_eur_put_RungeKutta(M(end),N,sigma,r,K,X,T);
toc;

fprintf('Fine discretization grids (M,N)=(%.0f,%.0f)  \n',[M(end);N]);

%%
err = zeros(length(M),1);
ratio = zeros(length(M),1);
for k = 1:length(M)
    u   = bs_eur_put_RungeKutta(M(k),N,sigma,r,K,X,T);
    u_temp = u_star(2^8/2^(k-1):2^8/2^(k-1):end);
    err(k) = norm(u-u_temp) / norm(u_temp);
    if k>1
        ratio(k) = log2( err(k-1) ./ err(k) );
    end
    fprintf('RungeKutta: (M,N)=(%.0f,%.0f) Error=%.2e Order=%.2f\n',[M(k);N;err(k);ratio(k)]);
end

% 验证时间方向二阶精度
disp('验证时间方向二阶精度')
M = 1000;
N =100*2.^(0:8);

tic;
u_star = bs_eur_put_RungeKutta(M,N(end),sigma,r,K,X,T);
toc;

fprintf('Fine discretization grids (M,N)=(%.0f,%.0f)  \n',[M;N(end)]);

%%
err = zeros(length(N),1);
ratio = zeros(length(N),1);
for k = 1:length(N)
    u = bs_eur_put_RungeKutta(M,N(k),sigma,r,K,X,T);
    err(k) = norm(u-u_star) / norm(u_star);
    if k>1
        ratio(k) = log2( err(k-1) ./ err(k) );
    end
    fprintf('RungeKutta: (M,N)=(%.0f,%.0f) Error=%.2e Order=%.2f\n',[M;N(k);err(k);ratio(k)]);
end

%% 欧式看涨期权计算 
% 验证空间方向二阶精度
disp('验证空间方向二阶精度')
M = 10*2.^(0:8);
N = 40000;

tic;
u_star  = bs_eur_call_RungeKutta(M(end),N,sigma,r,K,X,T);
toc;

fprintf('Fine discretization grids (M,N)=(%.0f,%.0f)  \n',[M(end);N]);

%%
err = zeros(length(M),1);
ratio = zeros(length(M),1);
for k = 1:length(M)
    u   = bs_eur_call_RungeKutta(M(k),N,sigma,r,K,X,T);
    u_temp = u_star(2^8/2^(k-1):2^8/2^(k-1):end);
    err(k) = norm(u-u_temp) / norm(u_temp);
    if k>1
        ratio(k) = log2( err(k-1) ./ err(k) );
    end
    fprintf('RungeKutta: (M,N)=(%.0f,%.0f) Error=%.2e Order=%.2f\n',[M(k);N;err(k);ratio(k)]);
end

% 验证时间方向二阶精度
disp('验证时间方向二阶精度')
M = 1000;
N =100*2.^(0:8)-1;

tic;
u_star = bs_eur_call_RungeKutta(M,N(end),sigma,r,K,X,T);
toc;

fprintf('Fine discretization grids (M,N)=(%.0f,%.0f)  \n',[M;N(end)]);

%%
err = zeros(length(N),1);
ratio = zeros(length(N),1);
for k = 1:length(N)
    u = bs_eur_call_RungeKutta(M,N(k),sigma,r,K,X,T);
    err(k) = norm(u-u_star) / norm(u_star);
    if k>1
        ratio(k) = log2( err(k-1) ./ err(k) );
    end
    fprintf('RungeKutta: (M,N)=(%.0f,%.0f) Error=%.2e Order=%.2f\n',[M;N(k);err(k);ratio(k)]);
end

