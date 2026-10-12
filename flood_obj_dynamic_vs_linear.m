function [F] = flood_obj_dynamic_vs_linear(Qdyn, Qobs, Qlin, opts)

% 目标函数：
% 以 1 - NSE 为主目标；
% 同时要求动态水库相对于线性水库逐场改进：
% 1）洪峰相对误差不大于线性水库，且尽量不超过 15%
% 2）最大洪峰峰现时间误差不大于线性水库，且尽量不超过 3 小时
% 3）NSE 不低于线性水库
%
% Qdyn    : 动态水库模拟流量过程
% Qobs    : 实测流量过程
% Qlin    : 线性水库模拟流量过程
% dt_hour : 时间步长，单位小时
%           逐小时数据：dt_hour = 1
%           30分钟数据：dt_hour = 0.5
%           15分钟数据：dt_hour = 0.25
% opts    : 可选参数结构体
%
% F       : 目标函数值，越小越好
% detail  : 输出各项指标，方便检查

%% 默认设置
dt_hour=1;
if nargin < 4 || isempty(dt_hour)
    dt_hour = 1;
end

if nargin < 5
    opts = struct();
end

% 主目标权重
wMain = get_opt(opts, 'wMain', 1.0);

% 惩罚权重
% 这些权重不要求加起来等于 1
% 因为这是“主目标 + 约束惩罚”的形式
wQ        = get_opt(opts, 'wQ', 50);    % 洪峰误差惩罚
wT        = get_opt(opts, 'wT', 30);    % 峰现时间惩罚
wNSEfloor = get_opt(opts, 'wNSEfloor', 100); % NSE 不低于线性结果的惩罚

% 工程约束
Q_limit_abs  = get_opt(opts, 'Q_limit_abs', 0.15); % 洪峰相对误差 15%
T_limit_hour = get_opt(opts, 'T_limit_hour', 3);   % 峰现时间误差 3 小时

% 容许误差
% 默认严格要求动态结果不差于线性结果
Q_tol   = get_opt(opts, 'Q_tol', 0.000);  % 洪峰误差容许值，例如 0.005 表示允许 0.5%
T_tol   = get_opt(opts, 'T_tol', 0.000);  % 峰现时间容许值，小时
NSE_tol = get_opt(opts, 'NSE_tol', 0.000);% NSE 容许下降值

% NSE 下降惩罚尺度
NSE_scale = get_opt(opts, 'NSE_scale', 0.02);

penalty_big = 1e6;
eps_val = 1e-6;

%% 数据整理
Qdyn = Qdyn(:);
Qobs = Qobs(:);
Qlin = Qlin(:);

if length(Qdyn) ~= length(Qobs) || length(Qlin) ~= length(Qobs)
    F = penalty_big;
    detail = [];
    return;
end

valid = isfinite(Qdyn) & isfinite(Qobs) & isfinite(Qlin);

Qdyn = Qdyn(valid);
Qobs = Qobs(valid);
Qlin = Qlin(valid);

n = length(Qobs);

if n < 3
    F = penalty_big;
    detail = [];
    return;
end

%% 计算动态水库和线性水库指标
mDyn = calc_metrics(Qdyn, Qobs, dt_hour);
mLin = calc_metrics(Qlin, Qobs, dt_hour);

%% 构造逐场约束目标

% 洪峰相对误差允许值：
% 要求动态水库 <= 线性水库；
% 如果线性水库本身超过 15%，则进一步要求向 15% 靠近
Q_allow = min(mLin.peakRelErrAbs, Q_limit_abs) + Q_tol;

% 峰现时间误差允许值：
% 要求动态水库 <= 线性水库；
% 如果线性水库本身超过 3 小时，则进一步要求向 3 小时靠近
T_allow = min(mLin.peakTimeErrAbsHour, T_limit_hour) + T_tol;

% NSE 下限：
% 要求动态水库 NSE 不低于线性水库 NSE
NSE_allow = mLin.NSE - NSE_tol;

%% 主目标：1 - NSE
F_main = 1 - mDyn.NSE;

%% 惩罚项 1：洪峰相对误差惩罚
% 动态洪峰误差小于允许值，不惩罚；
% 超过允许值，二次惩罚
P_Q = max(0, (mDyn.peakRelErrAbs - Q_allow) / max(Q_allow, 0.05))^2;

%% 惩罚项 2：峰现时间误差惩罚
% 动态峰现时间误差小于允许值，不惩罚；
% 超过允许值，二次惩罚
P_T = max(0, (mDyn.peakTimeErrAbsHour - T_allow) / max(T_allow, 1.0))^2;

%% 惩罚项 3：NSE 不得低于线性水库
P_NSEfloor = max(0, (NSE_allow - mDyn.NSE) / NSE_scale)^2;

%% 总目标函数
F = wMain * F_main + wQ * P_Q + wT * P_T + wNSEfloor * P_NSEfloor;

if ~isfinite(F)
    F = penalty_big;
end



end


