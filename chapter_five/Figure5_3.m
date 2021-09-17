clc
clear
close all
% 已知参数
Ta = 128;                      % 脉冲持续时间
Ka = -0.085;                   % 方位向调频率
% 参数计算
Delta_f_dop = abs(Ka*Ta);      % 方位信号带宽
% 参数设置
alpha_as = [5,0.25];           % 方位过采样率
Fa1 = alpha_as(1)*Delta_f_dop; % 方位采样频率PRF
N1 = 2*ceil(Fa1*Ta/2);         % 方位采样点数
dt1 = Ta/N1;                   % 采样时间间隔
df1 = Fa1/N1;                  % 采样频率间隔
Fa2 = alpha_as(2)*Delta_f_dop; % 方位采样频率PRF
N2 = 2*ceil(Fa2*Ta/2) ;        % 方位采样点数
dt2 = Ta/N2;                   % 采样时间间隔
df2 = Fa2/N2;                  % 采样频率间隔
% 变量设置
t1 = -Ta/2:dt1:Ta/2-dt1;       % 时间变量
t2 = -Ta/2:dt2:Ta/2-dt2;       % 时间变量
f1 = -Fa1/2:df1:Fa1/2-df1;     % 频率变量
f2 = -Fa2/2:df2:Fa2/2-df2;     % 频率变量
% 信号表达
st1 = exp(1j*pi*Ka*t1.^2);     % Chirp信号复数表达式
st2 = exp(1j*pi*Ka*t2.^2);     % Chirp信号复数表达式
% 参数计算
F1 = Ka*t1/Fa2;
F2 = (Ka*t2+floor((Fa2/2-Ka*t2)/Fa2)*Fa2)/Fa2;
% 绘图
H = figure();
set(H,'position',[100,100,800,600]);
subplot(311),plot(t1,real(st1),'k')
axis([-Ta/2-5 Ta/2+5,-1.2 1.2])
title('(a)混叠前方位Chirp信号实部'),ylabel('幅度')
subplot(312),plot(t2,real(st2),'k')
axis([-Ta/2-5 Ta/2+5,-1.2 1.2])
title('(b)混叠后方位Chirp信号实部'),ylabel('幅度')
subplot(313),plot(t1,F1,'K--',t2,F2,'K')
axis([-Ta/2-5 Ta/2+5,-2.2 2.2])
title('(c)信号瞬时频率'),xlabel('方位时间'),ylabel('频率(PRF)')
arrow([-16,1],[+16,1],'Color','k','Linewidth',1);
arrow([+16,1],[-16,1],'Color','k','Linewidth',1);
text(0,+1.5,'PRF时间','FontSize',14,'Color','red','HorizontalAlignment','center')
line([-16,-16],[-2,+2],'Color','r','LineStyle','--')
line([+16,+16],[-2,+2],'Color','r','LineStyle','--')
text(0,-1.5,'未混叠频率','FontSize',14,'Color','r','HorizontalAlignment','center')
text(-40,-1.5,'混叠频率','FontSize',14,'Color','r','HorizontalAlignment','center')
text(+40,-1.5,'混叠频率','FontSize',14,'Color','r','HorizontalAlignment','center')
sgtitle('图5.3 离散脉冲对方位信号的采样造成的方位混叠','Fontsize',20,'color','k')