function [IRW, PSLR, ISLR, IRW_null] = analyze_IR(s, show_flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   analyze_IR 脉冲冲激响应分析
%
%   输入：
%       s 脉冲响应信号序列（dB）
%       show_flag 是否展示寻找到的主瓣、最大旁瓣峰值点、3dB位置（默认false不展示）
%   输出：
%       IRW 脉冲3dB主瓣宽度（样点间隔数），注意部分对此有详细说明
%       PSLR 峰值旁瓣比（dB）
%       ISLR 积分旁瓣比
%   注意：
%       严格来说对于sinc函数：半功率点应该是-3.010299956639812dB；
%       而0.886/2点对应的是-3.011085043320321dB.
%       本程序按半功率点计算，因此此程序测得的分辨率要略大于按3dB宽度计算的分辨率
%       但又略小于按0.886/2点计算得到得分辨率.
%       **可以通过调整V_3dB的值来调整IRW测量标准**
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author: 
% create: 2021-09-04
% update: 2021-09-04
% version: v0.0.1
% Note: 未来可以对s先进行sinc插值，再测量相应参数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    show_flag = false;
end

% 参数预备
V_3dB = -3.010299956639812;  % 3dB值定义
L = length(s);               % 脉冲序列的长度
w = max(round(L/32), 1);     % 搜索窗长度
w_half = floor(w/2);
is_odd = w - 2*w_half;
ind_pattern = -w_half:(w_half-1+is_odd);

% 搜寻主瓣峰值点（带平滑搜索）
v_max = -Inf;
i_main = 1; % 主瓣位置索引
for n = 1:L
    ind = ind_pattern + n;
    if n <= w_half
        ind = ind(ind>0);
    end
    if n > L-w_half
        ind = ind(ind<=L);
    end
    
    v = sum(s(ind), 'all') / length(ind);
    if v > v_max
        v_max = v;
        i_main = n;
    end
end

% 搜寻主瓣左零点和3dB点
i_left = 1;
i_left_3dB = L;     % 令其初始值为L主要是为了标记左3dB点是否被找到
prev_v = s(i_main); % 上一个值
for n = i_main-1:-1:1
    if i_left_3dB == L && s(n) < V_3dB
        i_left_3dB = n + 0.5;
    end
    if s(n) > prev_v
        i_left = n+1;
        % 主瓣零点至少低于主瓣3dB
        if s(i_left) < s(i_main) - 3
            break
        end
    end
    prev_v = s(n);
end

% 搜寻主瓣右零点
i_right = L;
i_right_3dB = 1;     % 令其初始值为1主要是为了标记右3dB点是否被找到
prev_v = s(i_main); % 上一个值
for n = i_main+1:L
    if i_right_3dB == 1 && s(n) < V_3dB
        i_right_3dB = n - 0.5;
    end
    if s(n) > prev_v
        i_right = n-1;
        % 主瓣零点至少低于主瓣3dB
        if s(i_right) < s(i_main) - 3
            break
        end
    end
    prev_v = s(n);
end

% 1. 计算IRW
IRW = i_right_3dB - i_left_3dB; % 3dB主瓣宽度
IRW_null = i_right - i_left;    % 零点主瓣宽度

% 2. 计算PSLR
x = s([1:i_left,i_right:L]);
PSLR = max(x, [], 'all');

% 3. 计算ISLR
A_2 = 10.^(s/10);    % 幅值的平方，因为s的单位是dB
P_main = sum(A_2(i_left:i_right), 'all');
P_total = sum(A_2, 'all');
P_sidelobe = P_total - P_main;
ISLR = 10 * log10(P_sidelobe/P_main);

if (show_flag)
figure;
    plot(s);
    hold on;
    scatter(i_main, s(i_main));
    scatter(i_left, s(i_left));
    scatter(i_right, s(i_right));
    title('脉冲响应分析结果图');
    legend('脉冲响应归一化dB曲线', '主瓣峰值点','主瓣左零点','主瓣右零点');
    hold off;
end
end

