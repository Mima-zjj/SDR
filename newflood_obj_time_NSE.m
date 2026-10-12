function F = newflood_obj_time_NSE(Qsim, Qobs)

% 目标函数：
% 实测洪峰峰现时间误差 + 整体过程 NSE + 洪峰相对误差惩罚
%
% Qsim : 模拟流量过程，单列数组
% Qobs : 实测流量过程，单列数组
% F    : 目标函数值，越小越好

%% 数据整理
Qsim = Qsim(:);
Qobs = Qobs(:);

if length(Qsim) ~= length(Qobs)
    error('Qsim 和 Qobs 长度必须一致');
end

valid = isfinite(Qsim) & isfinite(Qobs);
Qsim = Qsim(valid);
Qobs = Qobs(valid);

n = length(Qobs);

if n < 3
    F = 100;
    return;
end

eps_val = 1e-6;

%% 权重设置
wT   = 0.70;   % 峰现时间误差权重
wNSE = 0.30;   % NSE 约束权重

% 新增：洪峰误差惩罚权重
% 这个权重不一定要和 wT、wNSE 加起来等于 1
% 如果洪峰仍然偏差大，可以把 wQ 调大，例如 5、10、20
wQ = 5.00;

NSE0 = 0.75;       % NSE 最低要求
QerrLimit = 0.15;  % 洪峰相对误差上限，15%

%% 计算 NSE
NSE = 1 - sum((Qobs - Qsim).^2) / ...
          (sum((Qobs - mean(Qobs)).^2) + eps_val);

% NSE 低于 NSE0 才惩罚
FNSE = max(0, (NSE0 - NSE) / (1 - NSE0 + eps_val))^2;

%% 识别实测主峰
minPromFactor = 0.08;
minPeakDistance = max(1, round(0.05 * n));
matchWindow     = max(1, round(0.10 * n));

promObs = max(eps_val, minPromFactor * (max(Qobs) - min(Qobs)));
promSim = max(eps_val, minPromFactor * (max(Qsim) - min(Qsim)));

[pObs, locObs] = findpeaks(Qobs, ...
    'MinPeakProminence', promObs, ...
    'MinPeakDistance', minPeakDistance);

% 如果没有识别到局部峰，则取最大洪峰
if isempty(pObs)
    [pObs, locObs] = max(Qobs);
end

pObs = pObs(:);
locObs = locObs(:);

%% 识别模拟峰
[pSim, locSim] = findpeaks(Qsim, ...
    'MinPeakProminence', promSim, ...
    'MinPeakDistance', minPeakDistance);

if isempty(pSim)
    [pSim, locSim] = max(Qsim);
end

pSim = pSim(:);
locSim = locSim(:);

Mobs = length(pObs);
Msim = length(pSim);

%% 实测洪峰权重
% 大洪峰权重大，小洪峰权重小
alpha = abs(pObs) ./ (sum(abs(pObs)) + eps_val);
alpha = alpha(:);

%% 峰现时间误差 + 洪峰误差匹配
usedSim = false(Msim, 1);

% 峰现时间误差
FT_each = ones(Mobs, 1);

% 洪峰相对误差惩罚
% 默认给 1，表示如果某个实测峰没有匹配到模拟峰，也要惩罚
FQ_each = ones(Mobs, 1);

for i = 1:Mobs

    dt = abs(locSim - locObs(i));
    dt(usedSim) = inf;

    [minDt, idx] = min(dt);

    if minDt <= matchWindow

        usedSim(idx) = true;

        %% 峰现时间误差
        ET = (locSim(idx) - locObs(i)) / matchWindow;
        FT_each(i) = ET^2;

        %% 洪峰相对误差惩罚
        EQ = abs(pSim(idx) - pObs(i)) / (abs(pObs(i)) + eps_val);

        % 15% 以内不惩罚，超过 15% 才惩罚
        FQ_each(i) = max(0, (EQ - QerrLimit) / QerrLimit)^2;

    end
end

%% 总峰现时间误差
FT = sum(alpha .* FT_each);

%% 总洪峰误差惩罚
FQ = sum(alpha .* FQ_each);

%% 总目标函数
F = wT * FT + wNSE * FNSE + wQ * FQ;

end
