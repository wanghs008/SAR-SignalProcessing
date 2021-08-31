clc
clear
close all

% 参数设置
B = 5.80e+6;           % 信号带宽
T = 7.26e-6;           % 脉冲持续时间
% 参数计算
K = B/T;               % 线性调频频率
alpha = 5;             % 过采样率
F = alpha*B;           % 采样频率
N = F*T;               % 采样点数
dt = T/N;              % 采样时间间隔
% 变量设置
t = -T/2:dt:T/2-dt;    % 时间变量

% % 设置自变量的几种方式
% t1 = linspace(-T/2,T/2-T/N,N)
% t2 = -T/2:dt:T/2-dt
% t3 = (-N/2:N/2-1)/N*T

% 信号表达
st = exp(1j*pi*K*t.^2);    % Chirp信号复数表达式
% 其他参数
f = K*t;                   % 瞬时频率
phi = pi*K*t.^2;           % 瞬时相位
% 绘图
figure
subplot(221),plot(t*1e+6,real(st))
title('(a)信号实部'),xlabel('相对于t_0时间(\mus)'),ylabel('幅度')
axis([-4 4,-1.2 1.2])
subplot(222),plot(t*1e+6,phi)
title('(b)信号相位'),xlabel('相对于t_0时间(\mus)'),ylabel('相位(弧度)')
axis([-4 4,-5 40])
subplot(223),plot(t*1e+6,imag(st))
title('(c)信号虚部'),xlabel('相对于t_0时间(\mus)'),ylabel('幅度')
axis([-4 4,-1.2 1.2])
subplot(224),plot(t*1e+6,f*1e-6)
title('(da)信号频率'),xlabel('相对于t_0时间(\mus)'),ylabel('MHz')
axis([-4 4,-3.2 3.2])
suptitle('图3.1 线性调频脉冲的相位和频率')