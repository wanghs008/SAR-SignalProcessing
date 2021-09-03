clc
clear
close all

% 参数设置
TBP = 100;             % 时间带宽积
T = 7.2e-6;            % 脉冲持续时间
t_0 = 1e-6;            % 脉冲回波时延
% 参数计算
B = TBP/T;             % 信号带宽
K = B/T;               % 线性调频频率
alpha_os = 1.25;       % 过采样率，使用较高的过采样率是为了提高采样频率
F = alpha_os*B;        % 采样频率
N = 2*ceil(F*T/2);     % 采样点数
dt = T/N;              % 采样时间间隔
df = F/N;              % 采样频率间隔
% 变量设置
t = -T/2+t_0:dt:T/2+t_0-dt;             % 时间变量
f = -F/2:df:F/2-df;                     % 频率变量
% 信号表达
st = exp(1j*pi*K*t.^2);                 % Chirp信号复数表达式
srt = exp(1j*pi*K*(t-t_0).^2);          % Chirp信号时延表达式
Srf = fft(srt);                         % Chirp信号频谱表达式
Hf = exp(1j*pi*f.^2/K);                 % 频域匹配滤波器
Soutf = Srf.*Hf;                        % 匹配滤波器输出
% 绘图
H = figure;
set(H,'position',[500,500,600,150]);
subplot(121),plot(f*1e-6,real(Soutf))
axis([-10 10,-15 15])
title('频谱实部'),xlabel('频率(MHz)'),ylabel('幅度')
subplot(122),plot(f*1e-6,imag(Soutf))
axis([-10 10,-15 15])
title('频谱虚部'),xlabel('频率(MHz)'),ylabel('幅度')
% suptitle('图3.9 匹配滤波后的信号频谱')