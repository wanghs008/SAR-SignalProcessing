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
Hf_3 = Window.*exp(1j*pi*f.^2/K); % 计算补零离散傅里叶变换
HF_phi = pi*f.^2/K;               % 频域匹配滤波器的相位
HF_fre = -F/2:df:F/2-df;          % 频域匹配滤波器的频率
% 绘图
H = figure;
set(H,'position',[300,300,600,600]);
subplot(311),plot(abs(Hf_3),'k')
axis([-15 791,0 1.2])
title('(a)加权后频域匹配滤波器的幅度'),ylabel('幅度')
line([362,362],[0,1.2],'Color','k','LineStyle','--')
line([414,414],[0,1.2],'Color','k','LineStyle','--')
arrow([130,0.3],[0  ,0.3],'Color','k','Linewidth',1);
arrow([232,0.3],[362,0.3],'Color','k','Linewidth',1);
text(181,0.3,'B/2','Color','red','HorizontalAlignment','center')
arrow([544,0.3],[414,0.3],'Color','k','Linewidth',1);
arrow([646,0.3],[776,0.3],'Color','k','Linewidth',1);
text(595,0.3,'B/2','Color','red','HorizontalAlignment','center')
line([388,388],[0,1.2],'Color','r','LineStyle','--')
subplot(312),plot(1:1:388,HF_fre(388:1:775)*1e-6,'k'),hold on
plot(388:1:775,HF_fre(1:1:388)*1e-6,'k')
title('(b)加权后频域匹配滤波器的频率'),ylabel('频率(MHz)')
axis([-15 791,-10 10])
line([388,388],[-10 10],'Color','r','LineStyle','--')
subplot(313),plot(fftshift(HF_phi),'k')
title('(c)加权后频域匹配滤波器的相位'),xlabel('频率(采样点)'),ylabel('弧度')
axis([-15 791,-50 700])
line([388,388],[-500,700],'Color','r','LineStyle','--')
suptitle('图3.12 方式三生成的匹配滤波器')