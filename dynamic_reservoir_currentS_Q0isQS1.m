function [QS3] = dynamic_reservoir_currentS_Q0isQS1(QSS3,U,dt,Q0,c_min,c_max,a,b,c1_init,Wref)
% dynamic_reservoir_currentS_Q0isQS1
%
% Q0 = QS3(1)
% 即第一时段出流量已知，从第二时段开始递推。
%
% 动态判别量：
%   W(i) = dt * Qin(i) + S3(i)
%   x(i) = W(i) / (W(i) + Wref)
%
% 动态消退系数：
%   c(i) = c_min + (c_max-c_min)/(1+exp(a*(x(i)-b)))
%
% 水量平衡递推：
%   QS3(i) = CS(i-1)/(CS(i)+dt)*QS3(i-1)
%          + dt/(CS(i)+dt)*Qin(i)
%
% 输入：
%   QSS3     - 入流序列
%   U        - 转换系数
%   dt       - 时间步长
%   Q0       - 初始出流量，即 QS3(1)
%   c_min    - 最小无量纲消退系数
%   c_max    - 最大无量纲消退系数
%   a        - Logistic 陡峭系数
%   b        - Logistic 阈值
%   c1_init  - 第一时段初始猜测消退系数，可省略
%   Wref     - 参考水量，可省略
%
% 输出：
%   QS3           - 模拟出流过程
%   S3            - 蓄水量过程
%   c             - 动态无量纲消退系数
%   CS            - 动态蓄量系数
%   x             - 当前水量强度指标
%   balance_error - 水量平衡误差，第1项为 NaN

    %% ===============================
    % 1. 输入处理
    %% ===============================

    QSS3 = QSS3(:);
    Qin  = QSS3 * U;

    n = length(Qin);

    QS3 = zeros(n,1);
    S3  = zeros(n,1);
    c   = zeros(n,1);
    CS  = zeros(n,1);
    x   = zeros(n,1);

    balance_error = nan(n,1);

    eps0 = 1e-12;

    if nargin < 9 || isempty(c1_init)
        c1_init = c_max;
    end

  
  



    %% ===============================
    % 2. 自动估计 Wref
    %% ===============================

    CS_init = c1_init * dt / (1 - c1_init);

    if nargin < 10 || isempty(Wref)
        Wref = max([max(dt * Qin), CS_init * abs(Q0), eps0]);

        if Wref <= eps0
            Wref = 1.0;
        end
    end

    %% ===============================
    % 3. 迭代参数
    %% ===============================

    maxIter = 50;
    tol     = 1e-8;
    omega   = 0.7;

    %% ===============================
    % 4. 第一时段：已知 QS3(1)=Q0
    %% ===============================

    QS3(1) = Q0;

    c_guess = c1_init;

    for iter = 1:maxIter

        CS_guess = c_guess * dt / (1 - c_guess);

        S_guess = CS_guess * QS3(1);

        W_current = dt * Qin(1) + S_guess;

        x_new = W_current / (W_current + Wref + eps0);

        arg = a * (x_new - b);
        arg = max(-60, min(60, arg));

        c_update = c_min + ...
            (c_max - c_min) / (1 + exp(arg));

        c_next = (1 - omega) * c_guess + omega * c_update;

        c_next = max(1e-6, min(1 - 1e-6, c_next));

        if abs(c_next - c_guess) < tol
            c_guess = c_next;
            break;
        end

        c_guess = c_next;
    end

    c(1)  = c_guess;
    CS(1) = c(1) * dt / (1 - c(1));
    S3(1) = CS(1) * QS3(1);

    W_current = dt * Qin(1) + S3(1);
    x(1) = W_current / (W_current + Wref + eps0);

    % 注意：因为 QS3(1)=Q0 是外部给定值，
    % 第1时段没有用水量平衡反推，所以 balance_error(1)=NaN

    %% ===============================
    % 5. 第二时段以后递推
    %% ===============================

    for i = 2:n

        c_guess = c(i-1);

        for iter = 1:maxIter

            CS_guess = c_guess * dt / (1 - c_guess);

            Q_guess = CS(i-1) / (CS_guess + dt) * QS3(i-1) ...
                    + dt / (CS_guess + dt) * Qin(i);

            S_guess = CS_guess * Q_guess;

            W_current = dt * Qin(i) + S_guess;

            x_new = W_current / (W_current + Wref + eps0);

            arg = a * (x_new - b);
            arg = max(-60, min(60, arg));

            c_update = c_min + ...
                (c_max - c_min) / (1 + exp(arg));

            c_next = (1 - omega) * c_guess + omega * c_update;

            c_next = max(1e-6, min(1 - 1e-6, c_next));

            if abs(c_next - c_guess) < tol
                c_guess = c_next;
                break;
            end

            c_guess = c_next;
        end

        c(i)  = c_guess;
        CS(i) = c(i) * dt / (1 - c(i));

        QS3(i) = CS(i-1) / (CS(i) + dt) * QS3(i-1) ...
               + dt / (CS(i) + dt) * Qin(i);

        S3(i) = CS(i) * QS3(i);

        W_current = dt * Qin(i) + S3(i);
        x(i) = W_current / (W_current + Wref + eps0);

        balance_error(i) = S3(i) - S3(i-1) ...
                         - dt * (Qin(i) - QS3(i));
    end

end
