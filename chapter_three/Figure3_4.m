clc
clear
close all

% 参数设置
TBP = 100;             % 时间带宽积
T = 10e-6;             % 脉冲持续时间
alpha_os = [1.4,1.2,1.0,0.8];             % 过采样率

% G = figure;
H = figure;
for i=1:length(alpha_os)
    % 参数计算
    B = TBP/T;             % 信号带宽
    K = B/T;               % 线性调频频率
    F = alpha_os(i)*B;     % 采样频率
    N = 2*ceil(F*T/2);     % 采样点数
    dt = T/N;              % 采样时间间隔
    df = F/N;              % 采样频率间隔
    % 变量设置
    t = -T/2:dt:T/2-dt;    % 时间变量
    f = -F/2:df:F/2-df;    % 频率变量
    f_zero = -F/2:F/(2*N):F/2-F/(2*N);    % 补零后的频率变零
    % 信号表达
    st = exp(1j*pi*K*t.^2);               % Chirp信号复数表达式
    Sf1 = fft(fftshift(st));              % Chirp信号频谱表达式
    st_zero = [st,zeros(1,N)];            % Chirp信号补零表达式
    Sf2 = fft(fftshift(st_zero));         % Chirp信号补零后的频谱表达式
%     % 绘图
%     figure(G);
%     subplot(4,2,2*i-1),plot(t*1e+6,real(st))
%     axis([-5 5,-1.2 1.2])
%     if(i==1)
%         title('信号实部')
%     end
%     if(i==4)
%         xlabel('时间(\mus)')
%     end
%     
%     % n1 = 0:1:N-1;
%     % subplot(4,2,2*i),plot(n1,abs(Sf1))
%     % axis([0 140,0 17])
%     
%     subplot(4,2,2*i),plot(f*1e-6,abs(Sf1))
%     axis([-4 4,0 22])
%     if(i==1)
%         title('频谱幅度')
%     end
%     if(i==4)
%         xlabel('频率单元(MHz)')
%     end
%     text(2.7,18,['\alpha_{os}= ',num2str(alpha_os(i))],'HorizontalAlignment','center')
    
    % 绘图
    figure(H);
    subplot(4,2,2*i-1),plot(t*1e+6,real(st))
    axis([-5 5,-1.2 1.2])
    if(i==1)
        title('信号实部')
    end
    if(i==4)
        xlabel('时间(\mus)')
    end
    
    % n2 = 0:1:2*N-1;
    % subplot(4,2,2*i),plot(n2,abs(Sf2))
    % axis([0 280,0 17])
    
    subplot(4,2,2*i),plot(f_zero*1e-6,abs(Sf2))
    axis([-4 4,0 22])
    if(i==1)
        title('频谱幅度')
    end
    if(i==4)
        xlabel('频率单元(MHz)')
    end
    text(2.7,18,['\alpha_{os}= ',num2str(alpha_os(i))],'HorizontalAlignment','center')
end
% figure(G);suptitle('图3.4 过采样率\alpha_{os}在频谱中引起的能量间隙')
figure(H);suptitle('图3.4 过采样率\alpha_{os}在频谱中引起的能量间隙')