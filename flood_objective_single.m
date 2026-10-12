function F = flood_objective_single(Qsim, Qobs)

% 只考虑最高峰洪峰误差、峰现时间误差和 NSE
% Qsim : 模拟流量，单列数组
% Qobs : 实测流量，单列数组
% F    : 目标函数值，越小越好

Qsim = Qsim(:);
Qobs = Qobs(:);

if length(Qsim) ~= length(Qobs)
    error('Qsim 和 Qobs 长度必须一致');
end

valid = ~isnan(Qsim) & ~isnan(Qobs);
Qsim = Qsim(valid);
Qobs = Qobs(valid);

n = length(Qobs);
if n < 1
    error('Qsim 和 Qobs 没有有效数据');
end

eps_val = 1e-6;

%% 权重设置
wQ   = 0.60;   % 最高洪峰相对误差权重
wT   = 0.30;   % 峰现时间误差权重
wNSE = 0.10;   % NSE 权重

%% 最高峰及峰现时间
[Qp_sim, kp_sim] = max(Qsim);
[Qp_obs, kp_obs] = max(Qobs);

%% 洪峰相对误差
EQ = (Qp_sim - Qp_obs) / (Qp_obs + eps_val);

%% 峰现时间误差，按总历时归一化
ET = (kp_sim - kp_obs) / n;

%% NSE 损失项
FNSE = sum((Qobs - Qsim).^2) / ...
       (sum((Qobs - mean(Qobs)).^2) + eps_val);

%% 目标函数
F = wQ * EQ^2 + wT * ET^2 + wNSE * FNSE;

end
