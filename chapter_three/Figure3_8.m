clc
clear
close all

% 参数设置
TBP = 42;              % 时间带宽积
T = 7.2e-6;            % 脉冲持续时间
tc = 1e-6;             % 时间偏移量
% 参数计算
B = TBP/T;             % 信号带宽
K = B/T;               % 线性调频频率
alpha_os = 50;         % 过采样率，使用较高的过采样率是为了提高采样频率
F = alpha_os*B;        % 采样频率
N = 2*ceil(F*T/2);     % 采样点数
dt = T/N;              % 采样时间间隔
% 变量设置
t = -T/2:dt:T/2-dt;    % 时间变量
t_out = linspace(2*t(1),2*t(end),2*length(t)-1);    % 循环卷积后的信号长度    
% 信号表达
st = exp(1j*pi*K*(t-tc).^2);          % Chirp信号复数表达式
ht = conj(fliplr(st));                % 时域匹配滤波器
sout = conv(st,ht);                   % 匹配滤波器输出
% 信号变换
sout_nor = sout/max(sout);                          % 单位化
sout_log = 20*log10(abs(sout)./max(abs(sout))+eps); % 归一化
% 绘图
figure
subplot(221),plot(t*1e+6,real(st))
axis([-4 4,-1.2 1.2])
title('(a)原始信号实部'),ylabel('幅度')

subplot(222),plot(t_out*1e+6,sout_log)
axis([-1 1,-30 5])
title('(c)压缩后信号(经扩展)'),ylabel('幅度')
[pslr] = get_pslr(sout_log);
text(0,3,['PSLR= ',num2str(pslr),'dB'],'HorizontalAlignment','center')

subplot(223),plot(t_out*1e+6,real(sout_nor))
axis([-4 4,-0.3 1.3])
title('(b)压缩后信号'),xlabel('相对于t_0时间(\mus)'),ylabel('弧度(dB)')
[hw,~,~] = get_hw(sout_log);
hw = hw*dt;
text(0,1.2,['HW= ',num2str(hw*1e+6),'\mus'],'HorizontalAlignment','center')

subplot(224),plot(t_out*1e+6,angle(sout_nor))
axis([-1 1,-5 5])
title('(d)压缩后信号相位(经扩展)'),xlabel('相对于t_0时间(\mus)'),ylabel('相位(弧度)')
suptitle('图3.8 非基带线性调频信号的匹配滤波')

%% HW函数
function [hw,locleft,locright] = get_hw(Af)
    % 找到Af的最大位置
    [~,locmax] = max(Af);
    % 找到locmax左边最接近-3dB的位置
    [~,locleft] = min(abs(Af(1:locmax)+3));
    % 找到locmax右边最接近-3dB的位置
    [~,locright] = min(abs(Af(locmax:end)+3));
    locright = locright + locmax - 1;
    % 得到3dB波束宽度
    hw = locright-locleft;
end
%% PSLR函数
function [pslr] = get_pslr(Af)
    % 找到所有的pesks
    peaks = findpeaks(Af);
    % 对peaks进行降序排列
    peaks = sort(peaks,'descend');
    % 得到第一旁瓣
    pslr = peaks(2);
end