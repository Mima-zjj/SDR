function F = flood_obj_time_NSE(Qsim, Qobs)

% 目标函数：
% 只考虑实测主峰的峰现时间误差 + 整体过程 NSE
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

valid = ~isnan(Qsim) & ~isnan(Qobs);
Qsim = Qsim(valid);
Qobs = Qobs(valid);

n = length(Qobs);

if n < 3
    error('数据长度太短，无法识别洪峰');
end

eps_val = 1e-6;

%% 权重设置
wT   = 0.70;   % 峰现时间误差权重
wNSE = 0.30;   % NSE 约束权重

NSE0 = 0.75;   % NSE 最低要求，可改为 0.80 或 0.85

%% 计算 NSE
NSE = 1 - sum((Qobs - Qsim).^2) / ...
          (sum((Qobs - mean(Qobs)).^2) + eps_val);

% NSE 低于 NSE0 才惩罚
FNSE = max(0, (NSE0 - NSE) / (1 - NSE0 + eps_val))^2;

%% 识别实测主峰
minPromFactor = 0.08;              % 主峰识别强度，越大识别峰越少
minPeakDistance = max(1, round(0.05 * n));
matchWindow     = max(1, round(0.10 * n));

[pObs, locObs] = findpeaks(Qobs, ...
    'MinPeakProminence', minPromFactor * (max(Qobs) - min(Qobs)), ...
    'MinPeakDistance', minPeakDistance);

% 如果没有识别到局部峰，则取最大洪峰
if isempty(pObs)
    [pObs, locObs] = max(Qobs);
end

%% 识别模拟峰
[pSim, locSim] = findpeaks(Qsim, ...
    'MinPeakProminence', minPromFactor * (max(Qsim) - min(Qsim)), ...
    'MinPeakDistance', minPeakDistance);

if isempty(pSim)
    [pSim, locSim] = max(Qsim);
end

Mobs = length(pObs);
Msim = length(pSim);

%% 实测主峰权重
% 大峰权重大，小峰权重小
alpha = pObs ./ (sum(pObs) + eps_val);

%% 峰现时间匹配
usedSim = false(Msim, 1);
FT_each = ones(Mobs, 1);   % 如果某个实测峰没有匹配到模拟峰，给最大时间惩罚

for i = 1:Mobs

    dt = abs(locSim - locObs(i));
    dt(usedSim) = inf;

    [minDt, idx] = min(dt);

    if minDt <= matchWindow

        usedSim(idx) = true;

        ET = (locSim(idx) - locObs(i)) / matchWindow;

        FT_each(i) = ET^2;
    end
end

%% 所有实测主峰的峰现时间误差
FT = sum(alpha .* FT_each);

%% 总目标函数
F = wT * FT + wNSE * FNSE;

end
