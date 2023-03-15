function demo_am_Implict
clear;clc

%% 模型参数 
sigma = 0.6;       % volatility
r = 0.25;        % interest rate 
K = 10;          % strike price
X = 50;               % asset value
T = 1;                % time

para_omg = 1.2;
%%
mer = 1e4;
sam = 1e2;
M = mer - 1;
N = mer;
[u,it,U] = bs_am_call_implict(M,N,sigma,r,K,X,T,para_omg);
U_sample = U(1:mer/sam:end, 1:mer/sam:end);
U_sample = U_sample(2:end-1, 1:end-1);
close all
surf(U_sample')
writematrix(U_sample,'bs_ame_call_implict.csv')
%% 计算欧式期权
% 验证空间方向二阶精度
fprintf(' \n');


M = 100*2.^(0:8)-1;
N = 40000;


disp('按下任意键开始参考解的计算 ');
pause;
tic;
[u_star,it]  = bs_eur_put_implict(M(end),N,sigma,r,K,X,T);
toc;

fprintf('Fine discretization grids (M,N)=(%.0f,%.0f) Iteration=%.2f  \n',[M(end);N;it]);

%%
disp('验证空间方向二阶精度')
err = zeros(length(M),1);
ratio = zeros(length(M),1);
for k = 1:length(M)
    u   = bs_eur_put_implict(M(k),N,sigma,r,K,X,T);
    u_temp = u_star(2^8/2^(k-1):2^8/2^(k-1):end);
    err(k) = norm(u-u_temp) / norm(u_temp);
    if k>1
        ratio(k) = log2( err(k-1) ./ err(k) );
    end
    fprintf('Implicit Euler: (M,N)=(%.0f,%.0f) Error=%.2e Order=%.2f\n',[M(k);N;err(k);ratio(k)]);
end

%% 验证时间方向一阶精度
fprintf(' \n');

M = 10000;
N =100*2.^(0:8)-1;

tic;
[u_star,it]  = bs_eur_put_implict(M,N(end),sigma,r,K,X,T);
toc;

fprintf('Fine discretization grids (M,N)=(%.0f,%.0f) Iteration=%.2f  \n',[M;N(end);it]);

%%
disp('验证时间方向一阶精度')
err = zeros(length(N),1);
ratio = zeros(length(N),1);
for k = 1:length(N)
    u   = bs_eur_put_implict(M,N(k),sigma,r,K,X,T);
    err(k) = norm(u-u_star) / norm(u_star);
    if k>1
        ratio(k) = log2( err(k-1) ./ err(k) );
    end
    fprintf('Implicit Euler: (M,N)=(%.0f,%.0f) Error=%.2e Order=%.2f\n',[M;N(k);err(k);ratio(k)]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 欧式看涨期权的计算

fprintf(' \n');
disp('开始欧式看涨期权的计算...');

%% 计算欧式看涨期权
% 验证空间方向二阶精度



M = 100*2.^(0:8);
N = 10000;


disp('按下任意键开始参考解的计算 ');
pause;
tic;
[u_star,it]  = bs_am_call_implict(M(end),N,sigma,r,K,X,T,para_omg);
toc;
fprintf('Fine discretization grids (M,N)=(%.0f,%.0f) Iteration=%.2f  \n',[M(end);N;it]);

%%
disp('验证空间方向二阶精度')
err = zeros(length(M),1);
ratio = zeros(length(M),1);
for k = 1:length(M)
    u   = bs_am_call_implict(M(k),N,sigma,r,K,X,T,para_omg);
    u_temp = u_star(2^8/2^(k-1):2^8/2^(k-1):end);
    err(k) = norm(u-u_temp) / norm(u_temp);
    if k>1
        ratio(k) = log2( err(k-1) ./ err(k) );
    end
    fprintf('Implicit Euler: (M,N)=(%.0f,%.0f) Error=%.2e Order=%.2f\n',[M(k);N;err(k);ratio(k)]);
end

%% 验证时间方向一阶精度
fprintf(' \n');

M = 10000;
N =100*2.^(0:8);

tic;
[u_star,it]  = bs_am_call_implict(M,N(end),sigma,r,K,X,T,para_omg);
toc;
fprintf('Fine discretization grids (M,N)=(%.0f,%.0f) Iteration=%.2f  \n',[M;N(end);it]);
%%
disp('验证时间方向一阶精度')
err = zeros(length(N),1);
ratio = zeros(length(N),1);
for k = 1:length(N)
    u   = bs_am_call_implict(M,N(k),sigma,r,K,X,T,para_omg);
    err(k) = norm(u-u_star) / norm(u_star);
    if k>1
        ratio(k) = log2( err(k-1) ./ err(k) );
    end
    fprintf('Implicit Euler: (M,N)=(%.0f,%.0f) Error=%.2e Order=%.2f\n',[M;N(k);err(k);ratio(k)]);
end

