clc
clear
close all

% 参数设置
TBP = 720;             % 时间带宽积
T = 10e-6;             % 脉冲持续时间
% 参数计算
B = TBP/T;             % 信号带宽
K = B/T;               % 线性调频频率
alpha = 1.25;          % 过采样率
F = alpha*B;           % 采样频率
N = 2*ceil(F*T/2);     % 采样点数
dt = T/N;              % 采样时间间隔
df = F/N;              % 采样频率间隔
% 变量设置
t = -T/2:dt:T/2-dt;    % 时间变量
f = -F/2:df:F/2-df;    % 频率变量
% 信号表达
st = exp(1j*pi*K*t.^2);               % Chirp信号复数表达式
Sf1 = exp(-1j*pi*f.^2/K);             % Chirp信号频谱表达式
Sf2 = fftshift(fft(fftshift(st)));    % Chirp信号频谱表达式
% 绘图
% figure
% subplot(221),plot(f*1e-6,real(Sf1))
% axis([-10 10,-1.2 1.2])
% title('(a)频谱实部'),xlabel('归一化后频率(Hz)'),ylabel('幅度')
% subplot(222),plot(f*1e-6,abs(Sf1))
% axis([-10 10,-1.2 1.2])
% title('(b)频谱相位'),xlabel('归一化后频率(Hz)'),ylabel('弧度')
% subplot(223),plot(f*1e-6,imag(Sf1))
% axis([-10 10,-1.2 1.2])
% title('(c)频谱虚部'),xlabel('归一化后频率(Hz)'),ylabel('幅度')
% subplot(224),plot(f*1e-6,unwrap(angle(Sf1)))
% axis([-28 28,0 400])
% title('(d)频谱相位'),xlabel('归一化后频率(Hz)'),ylabel('相位(弧度)')
% suptitle('图3.2 线性调频脉冲的复频谱')

figure
subplot(221),plot(f*1e-6,real(Sf2))
axis([-10 10,-40 40])
title('(a)频谱实部'),xlabel('归一化后频率(Hz)'),ylabel('幅度')
subplot(222),plot(f*1e-6,abs(Sf2))
axis([-50 50,0 40])
title('(b)频谱相位'),xlabel('归一化后频率(Hz)'),ylabel('弧度')
subplot(223),plot(f*1e-6,imag(Sf2))
axis([-10 10,-40 40])
title('(c)频谱虚部'),xlabel('归一化后频率(Hz)'),ylabel('幅度')
subplot(224),plot(f*1e-6,unwrap(angle(Sf2)))
axis([-50 50,0 900])
title('(d)频谱相位'),xlabel('归一化后频率(Hz)'),ylabel('相位(弧度)')
suptitle('图3.2 线性调频脉冲的复频谱')