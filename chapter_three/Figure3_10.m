clc
clear
close all

% 参数设置
TBP = 100;             % 时间带宽积
T = 7.2e-6;            % 脉冲持续时间
% 参数计算
B = TBP/T;             % 信号带宽
K = B/T;               % 线性调频频率
alpha_os = 1.25;       % 过采样率，使用较高的过采样率是为了提高采样频率
F = alpha_os*B;        % 采样频率
N = 2*ceil(F*T/2);     % 采样点数
dt = T/N;              % 采样时间间隔
df = F/N;              % 采样频率间隔
% 变量设置
t = -T/2:dt:T/2-dt;    % 时间变量
f = -F/2:df:F/2-df;    % 频率变量
% 信号表达
st = exp(1j*pi*K*t.^2);           % Chirp信号复数表达式
Sf = fft((st));                   % Chirp信号频谱表达式
Hf = exp(1j*pi*f.^2/K);           % 频域匹配滤波器
Soutf = Sf.*Hf;                   % 匹配滤波器输出
% 窗函数
window = kaiser(N,2.5)';          % 时域窗
Window = fftshift(window);        % 频域窗
% 信号变换
st_window = window.*exp(1j*pi*K*t.^2);          % 加窗后的Chirp信号
Hf_Window = Window.*Hf;                         % 加窗后的频域频谱滤波器
Soutf_Window = Hf_Window.*Sf;                   % 加窗后的匹配滤波器输出
% 绘图
H = figure;
set(H,'position',[500,500,600,300]);
subplot(221),plot(t*1e+6,window)
axis([-4 4,0 1.2])
title('时域窗函数'),ylabel('幅度')
subplot(222),plot(f*1e-6,Window)
axis([-10 10,0 1.2])
title('频域窗函数')
subplot(223),plot(t*1e+6,real(st_window))
axis([-4 4,-1.2 1.2])
title('加窗后的信号实部'),xlabel('时间(\mus)'),ylabel('幅度')
subplot(224),plot(f*1e-6,real(Soutf_Window))
axis([-10 10,-15 15])
title('加窗后的频谱实部'),xlabel('频率(MHz)')
suptitle('图3.10 Kaiser窗在时域和频域中的实现形式')