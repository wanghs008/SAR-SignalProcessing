clc
clear
close all
%% 参数设置
%  已知参数--》距离向参数
R_eta_c = 20e+3;                % 景中心斜距
Tr = 25e-6;                     % 发射脉冲时宽
Kr = 0.25e+12;                  % 距离向调频率
alpha_os_r = 1.2;               % 距离过采样率
Nrg = 256;                      % 距离线采样点数
%  计算参数--》距离向参数
Bw = abs(Kr)*Tr;                % 距离信号带宽
Fr = alpha_os_r*Bw;             % 距离向采样率
%  已知参数--》方位向参数
c = 3e+8;                       % 电磁传播速度
Vr = 150;                       % 等效雷达速度
Vs = Vr;                        % 卫星平台速度
Vg = Vr;                        % 波束扫描速度
f0 = 5.3e+9;                    % 雷达工作频率
Delta_f_dop = 80;               % 多普勒带宽
alpha_os_a = 1.3;               % 方位过采样率
Naz = 256;                      % 距离线数
theta_r_c = [0,-22.8]*pi/180;   % 波束斜视角
t_eta_c = [0,-51.7];            % 波束中心穿越时刻
f_eta_c = [0,+2055];            % 多普勒中心频率
%  计算参数--》方位向参数
lambda = c/f0;                  % 雷达工作波长
Fa = alpha_os_a*Delta_f_dop;    % 方位向采样率
%  参数计算
R0 = R_eta_c*cos(theta_r_c(2));                         % 最短斜距
La = 0.886*2*Vs*cos(theta_r_c(2))/Delta_f_dop;          % 实际天线长度
theta_bw = 0.886*lambda/La;                             % 方位向3dB波束宽度
Trr = Nrg/Fr;                   % 发射脉冲时宽
Taa = Naz/Fa;                   % 目标照射时间
Ka = 2*Vr^2*cos(theta_r_c(2))^2/lambda/R0;              % 方位向调频率
d_t_tau = 1/Fr;                 % 距离采样时间间隔
d_t_eta = 1/Fa;                 % 方位采样时间间隔
d_f_tau = Fa/Nrg;               % 距离采样频率间隔    
d_f_eta = Fa/Naz;               % 方位采样频率间隔
%% 变量设置
%  时间变量                                                    
t_tau = (-Trr/2:d_t_tau:Trr/2-d_t_tau) + 2*R_eta_c/c;   % 距离时间变量
t_eta = (-Taa/2:d_t_eta:Taa/2-d_t_eta) + t_eta_c(2);    % 方位时间变量
%  坐标设置                                                                                                           
[t_tauX,t_eta_Y] = meshgrid(t_tau,t_eta);               % 设置二维网络坐标
%% 信号设置
R_eta = sqrt(R0^2 + Vr^2*t_eta_Y.^2);                   % 瞬时斜率
A0 = 1;                                                 % 后向散射系数幅度
wr = (abs(t_tauX-2*R_eta/c) <= Tr/2);                   % 距离向包络
wa = sinc(0.886*atan(Vg*(t_eta_Y-t_eta_c(2))/R0)/theta_bw).^2;      % 方位向包络
%  接收信号
srt = A0*wr.*wa.*exp(-1j*4*pi*f0*R_eta/c)...
               .*exp(+1j*pi*Kr*(t_tauX-2*R_eta/c).^2);                                                           
srt_z = A0*wr.*wa.*exp(-1j*4*pi*f0*R_eta/c)...              
                 .*exp(+1j*pi*Kr*(t_tauX-2*R_eta/c).^2);% 正扫频
srt_f = A0*wr.*wa.*exp(-1j*4*pi*f0*R_eta/c)... 
                 .*exp(-1j*pi*Kr*(t_tauX-2*R_eta/c).^2);% 负扫频
%  距离时域-方位频域
Srf_rd = fft(srt);
%  距离频域-方位频域
SrF_2d = fft2(srt);
%% 绘图
%  距离时域-方位时域
H = figure();
set(H,'position',[100,100,900,300]);                    
subplot(131),imagesc(abs(srt_z)),colorbar               
xlabel('距离向(采样点)'),ylabel('方位向(采样点)'),title('(a)幅度');
subplot(132),imagesc(angle(srt_z))
xlabel('距离向(采样点)'),title('(b)相位(正扫频)');
subplot(133),imagesc(angle(srt_f))
xlabel('距离向(采样点)'),title('(c)相位(负扫频)');
sgtitle('图5.18 非零斜视角情况下单个点目标的时域性质','Fontsize',20,'color','k')