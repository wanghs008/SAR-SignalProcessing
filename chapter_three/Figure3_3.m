clc
clear
close all

% 参数设置
TBP = [25,50,100,200,400]; % 时间带宽积
T = 1e-6;                  % 脉冲持续时间

H = figure;
for i = 1:5
    % 参数计算
    B = TBP(i)/T;          % 信号带宽
    K = B/T;               % 线性调频频率
    alpha_os = 1.25;       % 过采样率
    F = alpha_os*B;        % 采样频率
    N = 2*ceil(F*T/2);     % 采样点数
    dt = T/N;              % 采样时间间隔
    df = F/N;              % 采样频率间隔
    % 变量设置
    t = -T/2:dt:T/2-dt;    % 时间变量
    f = -F/2:df:F/2-df;    % 频率变量
    % 信号表达
    st = exp(1j*pi*K*t.^2);               % Chirp信号复数表达式
    Sf = fftshift(fft(fftshift(st)));     % Chirp信号频谱表达式
    % 绘图
    figure(H);
    % 频谱幅度
    subplot(5,2,2*i-1)
    plot(f*1e-6,abs(Sf))
    axis([-F*1e-6/2-F*1e-6/100 F*1e-6/2+F*1e-6/100,0 8+6*(i-1)])
    if(i==5)
        xlabel('频率(MHz)')
    end
    ylabel('幅度')
    line([-B*1e-6/2,-B*1e-6/2],[0,sqrt(1/K)*1e+6*N],'color','k','linestyle','--')
    line([ B*1e-6/2, B*1e-6/2],[0,sqrt(1/K)*1e+6*N],'color','k','linestyle','--')
    line([-B*1e-6/2, B*1e-6/2],[sqrt(1/K)*1e+6*N,sqrt(1/K)*1e+6*N],'color','k','linestyle','--')
    % 频谱相位
    subplot(5,2,2*i)
    plot(f*1e-6,unwrap(angle(Sf))-max(unwrap(angle(Sf)))),hold on
    plot(f*1e-6,(-pi*f.^2/K)-max(-pi*f.^2/K),'k--');
    % set(gca,'YDir','reverse')    % 设置坐标轴翻转
    axis([-F*1e-6/2-F*1e-6/10 F*1e-6/2+F*1e-6/10,-32*2^(i-1) 0])
    if(i==5)
        xlabel('频率(MHz)')
    end
    ylabel('相位(弧度)')
    text(0,-TBP(i)/2,['TBP= ',num2str(TBP(i))],'HorizontalAlignment','center')
end
suptitle('图3.3 不同TBP值的离散傅里叶变换频谱变化')