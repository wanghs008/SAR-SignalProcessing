clc
clear
close all

% 参数设置
TBP = 100;             % 时间带宽积
T = 10e-6;             % 脉冲持续时间
% 参数计算
B = TBP/T;             % 信号带宽
K = B/T;               % 线性调频频率
alpha_os = 50;         % 过采样率，使用较高的过采样率是为了提高采样频率
F = alpha_os*B;        % 采样频率
N = 2*ceil(F*T/2);     % 采样点数
dt = T/N;              % 采样时间间隔
df = F/N;              % 采样频率间隔
% 变量设置
t = -T/2:dt:T/2-dt;    % 时间变量
f = -F/2:df:F/2-df;    % 频率变量
t_out = linspace(2*t(1),2*t(end),2*length(t)-1);    % 循环卷积后的信号长度    
% 信号表达
st = exp(1j*pi*K*t.^2);               % Chirp信号复数表达式
ht = conj(fliplr(st));                % 时域匹配滤波器
sout = conv(st,ht);                   % 匹配滤波器输出
sout = sout/max(sout);                % 归一化
% 绘图
figure
plot(t_out*1e+6,real(sout))
axis([-1 1,-0.4 1.2])
xlabel('时间(\mus)'),ylabel('幅度')
line([-1,1],[ 0,  0],'Color','k')
line([ 0,0],[-0.4,1.2],'Color','k')
line([-1,-0.05],[0.707,0.707],'Color','k','LineStyle','--')
line([ 0.05, 1],[0.707,0.707],'Color','k','LineStyle','--')
arrow([-0.3,0.707],[-0.05,0.707]);
arrow([ 0.3,0.707],[ 0.05,0.707]);
suptitle('图3.5 匹配滤波器输出的3dB分辨率的测量')