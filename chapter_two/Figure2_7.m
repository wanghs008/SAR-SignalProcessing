clc
close all
clear all

T = 10e-3;
F = 300;
t = 0:T/10000:T-T/10000;

A = sin(2*pi*F*t); % 原始信号

% figure
% plot(t,A,'k'),hold on
% title('不同采样率对300Hz的正弦波采样来说明混叠现象'),xlabel('t','Interpreter','latex'),ylabel('A','Interpreter','latex')
% grid on

% 不满足奈奎斯特采样定理
fs1 = 400;   % 设置采样频率
dt1 = 1/fs1; % 定义采样间隔

t1 = 0:dt1:T;

A1 = sin(2*pi*F*t1); % 采样信号

cs1 = spline(t1,A1);
xx1 = linspace(t1(1),t1(end),1000); % 插值点
yy1 = ppval(cs1,xx1); % 插值

% plot(t1,A1,'rp'),hold on
% plot(xx1,yy1,'r'),hold on

% 满足奈奎斯特采样定理
fs2 = 1600;  % 设置采样频率
dt2 = 1/fs2; % 定义采样间隔

t2 = 0:dt2:T;

A2 = sin(2*pi*F*t2);    % 采样信号

cs2 = spline(t2,A2);
xx2 = linspace(t2(1),t2(end),1000);     % 插值点
yy2 = ppval(cs2,xx2);   % 插值

% plot(t2,A2,'bh'),hold on
% plot(xx2,yy2,'b'),hold on

% set(legend,'Location','NorthEastOutside')
% legend('原始信号','','400Hz采样','','1600Hz采样')

%% 绘图
figure
plot(t,A,'k'),hold on
title('不同采样率对300Hz的正弦波采样来说明混叠现象'),xlabel('t','Interpreter','latex'),ylabel('A','Interpreter','latex')
plot(t1,A1,'rp'),hold on
plot(xx1,yy1,'r'),hold on
plot(t2,A2,'bh'),hold on
plot(xx2,yy2,'b'),hold on
grid on
%% 图例
set(legend,'Location','NorthEastOutside')
legend('原始信号','','400Hz采样','','1600Hz采样')
