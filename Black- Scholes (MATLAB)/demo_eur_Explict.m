function demo_eur_Explict
clear;clc

%% 模型参数 
sigma = 0.6;       % volatility
r = 0.25;        % interest rate 
K = 10;          % strike price
X = 50;               % asset value
T = 1;                % time


%% 计算欧式看跌期权

disp('验证空间方向二阶精度')
M = 10*2.^(0:8);
N = 10;

disp('按下任意键开始参考解的计算 ');
pause;

tic;
u_star  = bs_eur_put_Explict(M(end),N,sigma,r,K,X,T);
toc;

fprintf('Fine discretization grids (M,N)=(%.0f,%.0f)  \n',[M(end);N]);

%%
err = zeros(5,1);
ratio = zeros(5,1);
for k = 1:5
    u   = bs_eur_put_Explict(M(k),N,sigma,r,K,X,T);
    u_temp = u_star(2^8/2^(k-1):2^8/2^(k-1):end);
    err(k) = norm(u-u_temp) / norm(u_temp);
    if k>1
        ratio(k) = log2( err(k-1) ./ err(k) );
    end
    fprintf('Explict Euler: (M,N)=(%.0f,%.0f) Error=%.2e Order=%.2f\n',[M(k);N;err(k);ratio(k)]);
end

%% 验证时间方向一阶精度
disp('验证时间方向一阶精度')
M = 100;
N =100*2.^(0:8)-1;

tic;
u_star  = bs_eur_put_Explict(M,N(end),sigma,r,K,X,T);
toc;

fprintf('Fine discretization grids (M,N)=(%.0f,%.0f)  \n',[M;N(end)]);

%%
err = zeros(length(N),1);
ratio = zeros(length(N),1);
for k = 1:length(N)
    u   = bs_eur_put_Explict(M,N(k),sigma,r,K,X,T);    
    err(k) = norm(u-u_star) / norm(u_star);
    if k>1
        ratio(k) = log2( err(k-1) ./ err(k) );
    end
    fprintf('Implicit Euler: (M,N)=(%.0f,%.0f) Error=%.2e Order=%.2f\n',[M;N(k);err(k);ratio(k)]);
end



%%%%%%%     
%% 计算欧式看涨期权

disp('验证空间方向二阶精度')
M = 10*2.^(0:8);
N = 10;

disp('按下任意键开始参考解的计算 ');
pause;

tic;
u_star  = bs_eur_call_Explict(M(end),N,sigma,r,K,X,T);
toc;

fprintf('Fine discretization grids (M,N)=(%.0f,%.0f)  \n',[M(end);N]);

%%
err = zeros(5,1);
ratio = zeros(5,1);
for k = 1:5
    u   = bs_eur_call_Explict(M(k),N,sigma,r,K,X,T);
    u_temp = u_star(2^8/2^(k-1):2^8/2^(k-1):end);
    err(k) = norm(u-u_temp) / norm(u_temp);
    if k>1
        ratio(k) = log2( err(k-1) ./ err(k) );
    end
    fprintf('Explict Euler: (M,N)=(%.0f,%.0f) Error=%.2e Order=%.2f\n',[M(k);N;err(k);ratio(k)]);
end

%% 验证时间方向一阶精度
disp('验证时间方向一阶精度')
M = 100;
N =100*2.^(0:8)-1;

tic;
u_star  = bs_eur_call_Explict(M,N(end),sigma,r,K,X,T);
toc;

fprintf('Fine discretization grids (M,N)=(%.0f,%.0f)  \n',[M;N(end)]);

%%
err = zeros(length(N),1);
ratio = zeros(length(N),1);
for k = 1:length(N)
    u   = bs_eur_call_Explict(M,N(k),sigma,r,K,X,T);    
    err(k) = norm(u-u_star) / norm(u_star);
    if k>1
        ratio(k) = log2( err(k-1) ./ err(k) );
    end
    fprintf('Implicit Euler: (M,N)=(%.0f,%.0f) Error=%.2e Order=%.2f\n',[M;N(k);err(k);ratio(k)]);
end


