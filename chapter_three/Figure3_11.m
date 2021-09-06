clc
clear
close all

% 参数设置
TBP = 724;             % 时间带宽积
T = 42e-6;             % 脉冲持续时间
Nfft = 2^11;           % FFT长度
% 参数计算
B = TBP/T;             % 信号带宽
K = B/T;               % 线性调频频率
alpha_os = 1.07;       % 过采样率
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
% 窗函数
window = kaiser(N,2.5)';          % 时域窗
Window = fftshift(window);        % 频域窗
% 信号变换
ht = st;                          % 复制信号
ht_window = window.*ht;           % 加窗
Hf_2 = conj(fft(ht_window,Nfft));       % 计算补零离散傅里叶变换
% 绘图
H = figure;
set(H,'position',[500,500,600,400]);
subplot(211),plot(abs(Hf_2),'k')
axis([-50 2100,0 35])
title('(a)加权后频域匹配滤波器的幅度谱'),ylabel('幅度')
line([ 860, 860],[0,35],'Color','k','LineStyle','--')
line([1190,1190],[0,35],'Color','k','LineStyle','--')
line([990 ,990 ],[0,35],'Color','r','LineStyle','--')
line([1060,1060],[0,35],'Color','r','LineStyle','--')
subplot(212),plot(1:1:990,unwrap(-angle(Hf_2(1:1:990))),'k'),hold on
plot(1060:1:2048,unwrap(-angle(Hf_2(1060:1:2048))),'k')
axis([-50 2100,-2000 600])
title('(b)加权后频域匹配滤波器的相位谱'),xlabel('频率(FFT采样点)'),ylabel('弧度')
line([990 ,990 ],[-2000,600],'Color','r','LineStyle','--')
line([1060,1060],[-2000,600],'Color','r','LineStyle','--')
suptitle('图3.11 方式二生成的频域匹配滤波器频率响应函数的幅度和相位')