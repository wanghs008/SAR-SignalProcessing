clc
clear
close all
% 已知参数--》距离向参数
R_eta_c_1 = 450;                            % 景中心斜距
R_eta_c_2 = 850;                            % 景中心斜距
R_eta_c_3 = 1250;                           % 景中心斜距
% 已知参数--》方位向参数
c  = 3e+8;                                  % 电磁传播速度
f0 = 5.3e+6;                                % 雷达工作频率
Vr = 7100;                                  % 等效雷达速度
Ta = 0.64;                                  % 目标照射时间
Ka = 2095;                                  % 方位向调频率
theta_r_c = 0*pi/180;                       % 斜视角
% 参数计算
lambda = c/f0;                              % 雷达工作波长
R0_1 = R_eta_c_1*cos(theta_r_c);            % 最短斜距
R0_2 = R_eta_c_2*cos(theta_r_c);            % 最短斜距
R0_3 = R_eta_c_3*cos(theta_r_c);            % 最短斜距
Delta_f_dop = abs(Ka*Ta);                   % 方位信号带宽
t_eta_c_1 = -R_eta_c_1*sin(theta_r_c)/Vr;   % 波束中心穿越时刻
t_eta_c_2 = -R_eta_c_2*sin(theta_r_c)/Vr;   % 波束中心穿越时刻
t_eta_c_3 = -R_eta_c_3*sin(theta_r_c)/Vr;   % 波束中心穿越时刻
% 参数设置
alpha_a_s = 1;                              % 方位过采样率
Fa = alpha_a_s*Delta_f_dop;                 % 方位采样频率PRF
Na = 2*ceil(Fa*Ta/2);                       % 方位采样点数
dt = Ta/Na;                                 % 采样时间间隔
df = Fa/Na;                                 % 采样频率间隔 
% 变量设置
t_eta =  (-Ta/16:dt:Ta/16-dt);              % 方位时间变量
f_eta =  (-Fa/16:df:Fa/16-df);              % 方位频率变量
% 信号表达
R_eta_1 = R0_1 + Vr^2*t_eta.^2/(2*R0_1);    % 目标轨迹
R_eta_2 = R0_2 + Vr^2*t_eta.^2/(2*R0_2);    % 目标轨迹
R_eta_3 = R0_3 + Vr^2*t_eta.^2/(2*R0_3);    % 目标轨迹
% 信号表达                                          
R_rd_1 = R0_1 + lambda^2*R0_1/(8*Vr^2)*f_eta.^2;    % 目标轨迹
R_rd_2 = R0_2 + lambda^2*R0_2/(8*Vr^2)*f_eta.^2;    % 目标轨迹
R_rd_3 = R0_3 + lambda^2*R0_3/(8*Vr^2)*f_eta.^2;    % 目标轨迹
% 绘图                                                                
G = figure();                                                           
set(G,'position',[100,100,800,600]);
plot(R_eta_1,t_eta,'k','LineWidth',2),hold on
plot(R_eta_2,t_eta,'k','LineWidth',2),hold on
plot(R_eta_3,t_eta,'k','LineWidth',2),set(gca,'ydir','reverse')
axis([400 1350,-0.05 +0.05])
xlabel('距离(m)'),ylabel('方位时间(s)')
sgtitle('时域中的目标轨迹','Fontsize',20,'color','k')
% 绘图                                                                
G = figure();                                                           
set(G,'position',[100,100,800,600]);
plot(R_rd_1,f_eta,'k','LineWidth',2),hold on
plot(R_rd_2,f_eta,'k','LineWidth',2),hold on
plot(R_rd_3,f_eta,'k','LineWidth',2),set(gca,'ydir','reverse')
axis([400 1350,-110 +110])
xlabel('距离(m)'),ylabel('方位频率(Hz)')
sgtitle('距离多普勒域中的目标轨迹','Fontsize',20,'color','k')