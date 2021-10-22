clc
close all
clear all

% LFM信号 实信号 s(t) = cos(2*pi*(f0.*t + 0.5*K*t.^2))
% LFM信号 复信号 s(t) = exp(1i*2*pi*(f0.*t + 0.5*K*t.^2))

% 参数设置
T = 1;
B = 150;
K = B/T;

% figure('Name','...','NumberTitle','off')
% figure属性设置
% set(gca,'XTick',[xmin:间距:xmax])
% set(gca,'YTick',[ymin:间距:ymax])
% grid on
% set(gca,'Xgrid','off')
% set(gca,'Ygrid','off')
% box off
% legend('lenname1','lenname2'),set(legend,'Location','NorthEastOutside'),set(legend,'Position',[xmin,xmax,ymin.ymax])

set(figure,'position',[100,100,900,800]);

for i = 1:1:5
    
    f0 = -75+(i-1)*100;

    % 原始信号
    dt = 1/2000;
    t1 = 0:dt:T-dt;
    
    At = t1+2;
    At = At./max(max(At));  % 信号幅度修正
    
    st = exp(1i*2*pi*(f0*t1+0.5*K*t1.^2)).*At;
    
    % 傅里叶变换
    N1 = 2^20;
    df = 1/dt;
    f1 = -df/2:df/N1:df/2-df/N1;
    
    Sf = fft(st,N1);
    Sf = fftshift(abs(Sf)/max(abs(Sf)));
    
    % 采样信号
    Fs = 400;
    Dt = 1/Fs;
    t2 = 0:Dt:T-Dt;
    
    Bt = t2+2;
    Bt = Bt./max(max(Bt));
    
    st_c = exp(1i*2*pi*(f0*t2+0.5*K*t2.^2)).*Bt;
    
    % 傅里叶变换
    N2 = 2^20;
    f2 = -Fs/2:Fs/N2:Fs/2-Fs/N2;
    
    Sf_c = fft(st_c,N2);
    Sf_c = fftshift(abs(Sf_c)/max(abs(Sf_c)));
    
    % figure 
    % plot(t1,st,'b')
    % xlabel('t'),ylabel('st')
    %grid on
    
    % figure
    % plot(f1,Sf,'b')
    % xlabel('f'),ylabel('Sf')
    % grid on
    
    % figure
    % plot(f2,Sf_c,'b')
    % xlabel('f'),ylabel('Sf_c')
    % grid on
    
    subplot(5,3,[3*(i-1)+1,3*(i-1)+2]),plot(f1,Sf,'b')
    xlabel('频率'),ylabel('幅度')
    axis([-200 600, 0 1])
    grid on
    subplot(5,3,3*i)
    plot(f2,Sf_c,'b')
    xlabel('频率'),ylabel('幅度')
    axis([-200 200, 0 1])
    grid on
end

suptitle('图2.8 采样引起的频谱平移(复信号)')