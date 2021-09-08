clc
clear
close all

% 参数设置
TBP = 42;                % 时间带宽乘积
T = 7.2e-6;              % 脉冲持续时间
N_st = 2^11;                         
% 参数计算
B = TBP/T;               % 信号带宽
K = B/T;                 % 调频频率
alpha_os = 200;          % 过采样率，使用过高的过采样率是为了方便地实现升采样
F = alpha_os*B;          % 采样频率
N = 2*ceil(F*T/2);       % 采样点数
dt = T/N;                % 采样时间间隔
df = F/N;                % 采样频率间隔
% 参数设置
t_c = 0e-6;              % 时间偏移
% 参数计算
f_c = -K*t_c;            % 中心频点
% 变量设置
t1 = -T/2:dt:T/2-dt;     % 时间变量
f1 = -F/2:df:F/2-df;     % 频率变量
% 信号表达                             
st = exp(1j*pi*K*(t1-t_c).^2);          % 发射信号
% % 绘图
% H1 = figure;
% set(H1,'position',[100,100,600,300]);
% subplot(211),plot(t1,real(st),'k')
% subplot(212),plot(t1,imag(st),'k')
% suptitle('发射信号')
% 参数设置
t_0 = 0e-6;              % 回波时延
% 变量设置
t2 = -T/2+t_0:dt:T/2+t_0-dt;            % 时间变量
f2 = -F/2+f_c:df:F/2+f_c-df;            % 频率变量
% 信号表达                                                                 
srt = exp(1j*pi*K*(t2-t_c-t_0).^2);     % 回波信号
% % 绘图
% H2 = figure;
% set(H2,'position',[100,100,600,300]);
% subplot(211),plot(t2,real(srt),'k')
% subplot(212),plot(t2,imag(srt),'k')
% suptitle('回波信号')
%%
% % 窗函数
% window_1 = kaiser(N,2.5)';              % 时域窗
% Window_1 = fftshift(window_1);          % 频域窗
% % 信号变换-->方式一
% ht_1 = conj(fliplr(srt));               % 将时间反褶后的复制脉冲取复共轭
% ht_window_1 = window_1.*ht_1;           % 加窗
% Hf_1 = fftshift(fft(ht_1,N));           % 计算补零离散傅里叶变换
% % 绘图
% H3 = figure;
% set(H3,'position',[100,100,600,450]);
% subplot(311),plot(real(Hf_1),'k')
% axis([4000 4400,-1200 1600])
% subplot(312),plot(imag(Hf_1),'k')
% axis([4000 4400,-1200 1600])
% subplot(313),plot( abs(Hf_1),'k')
% suptitle('方式一生成的频域匹配滤波器')
% axis([4000 4400,    0 1600])
% % 窗函数 
% window_2 = kaiser(N,2.5)';              % 时域窗
% Window_2 = fftshift(window_2);          % 频域窗
% % 信号变换-->方式二
% ht_2 = srt;                             % 复制信号
% ht_window_2 = window_2.*ht_2;           % 加窗
% Hf_2 = fftshift(conj(fft(ht_2,N)));     % 计算补零离散傅里叶变换
% % 绘图
% H3 = figure;
% set(H3,'position',[100,100,600,450]);
% subplot(311),plot(real(Hf_2),'k')
% axis([4000 4400,-1200 1600])
% subplot(312),plot(imag(Hf_2),'k')
% axis([4000 4400,-1200 1600])
% subplot(313),plot( abs(Hf_2),'k')
% axis([4000 4400,    0 1600])
% suptitle('方式二生成的频域匹配滤波器')
% % 窗函数
% window_3 = kaiser(N,2.5)';              % 时域窗
% Window_3 = fftshift(window_3);          % 频域窗
% % 信号变换-->方式三
% Hf_3 = Window_3.*exp(1j*pi*f2.^2/K);    % 计算补零离散傅里叶变换
% % 绘图
% H3 = figure;
% set(H3,'position',[100,100,600,450]);
% subplot(311),plot(real(Hf_3),'k')
% axis([2000 6400,-1.2 1.2])
% subplot(312),plot(imag(Hf_3),'k')
% axis([2000 6400,-1.2 1.2])
% subplot(313),plot( abs(Hf_3),'k')
% axis([2000 6400,-0.2 1.2])
% suptitle('方式三生成的频域匹配滤波器')
% 信号表达
Srf = fftshift(fft(srt));
% % 绘图
% H4 = figure;
% set(H4,'position',[100,100,600,300]);
% subplot(211),plot(real(srt),'k')
% axis([0 N,-1.2 1.2])
% subplot(212),plot(abs(Srf),'k')
% axis([4000 4400,-50 1600])
% suptitle('回波信号的频谱分析')
% % 信号表达
% Soutf_1 = Srf.*Hf_1;
% soutt_1 = ifft(ifftshift(Soutf_1));     % 方式一匹配滤波结果
% soutt_1_nor = abs(soutt_1)./max(abs(soutt_1));               % 归一化
% soutt_1_log = 20*log10(abs(soutt_1)./max(abs(soutt_1))+eps); % 对数化
% Soutf_2 = Srf.*Hf_2;
% soutt_2 = ifft(ifftshift(Soutf_2));     % 方式二匹配滤波结果
% soutt_2_nor = abs(soutt_2)./max(abs(soutt_2));               % 归一化
% soutt_2_log = 20*log10(abs(soutt_2)./max(abs(soutt_2))+eps); % 对数化
% Soutf_3 = Srf.*Hf_3;
% soutt_3 = ifft(ifftshift(Soutf_3));     % 方式三匹配滤波结果 
% soutt_3_nor = abs(soutt_3)./max(abs(soutt_3));               % 归一化
% soutt_3_log = 20*log10(abs(soutt_3)./max(abs(soutt_3))+eps); % 对数化
% % 参数计算-->IRW
% [irw1,a1,b1] = get_irw(fftshift(soutt_1_nor));irw11 = irw1*dt;
% [irw2,a2,b2] = get_irw(fftshift(soutt_2_nor));irw12 = irw2*dt;
% [irw3,a3,b3] = get_irw(soutt_3_nor);irw13 = irw3*dt;
% % 参数计算-->ISLR
% [islr1] = get_islr(fftshift(soutt_1_nor),5);
% [islr2] = get_islr(fftshift(soutt_2_nor),5);
% [islr3] = get_islr(soutt_3_nor,5);
% % 绘图                                        
% H5 = figure;
% set(H5,'position',[100,100,600,600]);
% subplot(411),plot(real(srt),'k')
% axis([0 N,-1.2 1.2])
% subplot(412),plot(fftshift(soutt_1_nor),'k')
% axis([0 N,-0.1,1.2])
% line([4110,4110],[-0.1,1.2],'Color','r','LineStyle','--')
% line([4290,4290],[-0.1,1.2],'Color','r','LineStyle','--')
% text(1500,0.7,['IRW= ',num2str(irw11*1e+6),'\mus'],'HorizontalAlignment','center')
% text(7000,0.7,['ISLR= ',num2str(islr1),'dB'],'HorizontalAlignment','center')
% subplot(413),plot(fftshift(soutt_2_nor),'k')
% axis([0 N,-0.1,1.2])
% line([4111,4111],[-0.1,1.2],'Color','r','LineStyle','--')
% line([4291,4291],[-0.1,1.2],'Color','r','LineStyle','--')
% text(1500,0.7,['IRW= ',num2str(irw12*1e+6),'\mus'],'HorizontalAlignment','center')
% text(7000,0.7,['ISLR= ',num2str(islr2),'dB'],'HorizontalAlignment','center')
% subplot(414),plot(soutt_3_nor,'k')
% axis([0 N,-0.1,1.2])
% line([4113,4113],[-0.1,1.2],'Color','r','LineStyle','--')
% line([4289,4289],[-0.1,1.2],'Color','r','LineStyle','--')
% text(1500,0.7,['IRW= ',num2str(irw13*1e+6),'\mus'],'HorizontalAlignment','center')
% text(7000,0.7,['ISLR= ',num2str(islr3),'dB'],'HorizontalAlignment','center')
% % 参数计算-->PSLR
% [pslr1] = get_pslr(fftshift(soutt_1_log));
% [pslr2] = get_pslr(fftshift(soutt_2_log));
% [pslr3] = get_pslr(soutt_3_log);
% % 绘图                                        
% H6 = figure;
% set(H6,'position',[100,100,600,600]);
% subplot(411),plot(real(srt),'k')
% axis([0 N,-1.2 1.2])
% subplot(412),plot(fftshift(soutt_1_log),'k')
% axis([0 N,-35 0])
% text(7000,-10,['PSLR= ',num2str(pslr1),'dB'],'HorizontalAlignment','center')
% subplot(413),plot(fftshift(soutt_2_log),'k')
% axis([0 N,-35 0])
% text(7000,-10,['PSLR= ',num2str(pslr2),'dB'],'HorizontalAlignment','center')
% subplot(414),plot(soutt_3_log,'k')
% axis([0 N,-35 0])
% text(7000,-10,['PSLR= ',num2str(pslr3),'dB'],'HorizontalAlignment','center')
%%
% 参数设置
QPE = linspace(0,0.8*pi,N);            % 二次相位误差
dk = QPE/(pi*(T/2)^2);                 % 调频率误差
% 参数设置
IRW1  = zeros(1,N);                    % 初始化冲激响应宽度
PSLR1 = zeros(1,N);                    % 初始化峰值旁瓣比
ISLR1 = zeros(1,N);                    % 初始化积分旁瓣比
IRW2  = zeros(1,N);                    % 初始化冲激响应宽度
PSLR2 = zeros(1,N);                    % 初始化峰值旁瓣比
ISLR2 = zeros(1,N);                    % 初始化积分旁瓣比
IRW3  = zeros(1,N);                    % 初始化冲激响应宽度
PSLR3 = zeros(1,N);                    % 初始化峰值旁瓣比
ISLR3 = zeros(1,N);                    % 初始化积分旁瓣比
% 显示运行时间
tic
% 显示进度条框
wait_title = waitbar(0,'Program Initializing ...');                                            
% 绘图
H8 = figure;
set(H8,'position',[100,100,600,350]);
% 循环计算
pause(1);
for i = 1:N
    % 变量设置
    B_dk = (K+dk(i))*T;
    F_dk = alpha_os*B_dk;
    df_dk = F_dk/N;
    f3 = -F_dk/2+f_c:df_dk:F_dk/2+f_c-df_dk;           % 频率变量
    % 信号表达                                                                 
    st_dk = exp(1j*pi*(K+dk(i))*t1.^2);                  
    Sf_dk = fftshift(fft(st_dk));
    % 信号变换-->频域方式一
    window_1 = kaiser(N,2.5)';                         % 时域窗
    Window_1 = fftshift(window_1);                     % 频域窗
    ht_dk_1 = conj(fliplr(st_dk));                     % 将时间反褶后的复制脉冲取复共轭
    ht_window_dk_1 = window_1.*ht_dk_1;                % 加窗
    Hf_dk_1 = fftshift(fft(ht_window_dk_1,N));         % 计算补零离散傅里叶变换
    % 信号变换-->频域方式二
    window_2 = kaiser(N,2.5)';                         % 时域窗
    Window_2 = fftshift(window_2);                     % 频域窗
    ht_dk_2 = st_dk;                                   % 复制信号
    ht_window_dk_2 = window_2.*ht_dk_2;                % 加窗
    Hf_dk_2 = fftshift(conj(fft(ht_window_dk_2,N)));   % 计算补零离散傅里叶变换
    % 信号变换-->频域方式三
    window_3 = kaiser(N,2.5)';                         % 时域窗
    Window_3 = fftshift(window_3);                     % 频域窗
    Hf_dk_3 = Window_3.*exp(1j*pi*f3.^2/(K+dk(i)));    % 计算补零离散傅里叶变换
% 参数计算-->方式一                                         
    Soutf_dk_1 = Srf.*Hf_dk_1;
    soutt_dk_1 = ifft(ifftshift(Soutf_dk_1));          % 方式一匹配滤波结果 
    soutt_dk_1_nor = abs(soutt_dk_1)./max(abs(soutt_dk_1));               % 归一化
    soutt_dk_1_log = 20*log10(abs(soutt_dk_1)./max(abs(soutt_dk_1))+eps); % 对数化
    % 参数计算-->IRW
    [irw_dk_1,~,~] = get_irw(fftshift(soutt_dk_1_nor));
    IRW1(i) = irw_dk_1;
    % 参数计算-->PSLR
    [pslr_dk_1] = get_pslr(fftshift(soutt_dk_1_log));
    PSLR1(i) = pslr_dk_1;
    % 参数计算-->ISLR
    [islr_dk_1] = get_islr(fftshift(soutt_dk_1_nor),5);
    ISLR1(i) = islr_dk_1;
     % 绘图
    if i == 2990
       figure(H8);
       subplot(121),plot(fftshift(soutt_dk_1_log),'k')
       axis([3600 4800,-25 -15])
       title('(a)|QPE|略小于0.28\pi弧度'),xlabel('采样点序列'),ylabel('幅度/dB')
    end
    if i == 3100
       figure(H8);
       subplot(122),plot(fftshift(soutt_dk_1_log),'k')
       axis([3600 4800,-25 -15])
       title('(b)|QPE|略大于0.28\pi弧度'),xlabel('采样点序列'),ylabel('幅度/dB')
       suptitle('图3.15 最大旁瓣位置不同，而脉冲响应相似时的情况')
    end
% % 参数计算-->方式二                                        
%     Soutf_dk_2 = Srf.*Hf_dk_2;
%     soutt_dk_2 = ifft(ifftshift(Soutf_dk_2));          % 方式二匹配滤波结果 
%     soutt_dk_2_nor = abs(soutt_dk_2)./max(abs(soutt_dk_2));               % 归一化
%     soutt_dk_2_log = 20*log10(abs(soutt_dk_2)./max(abs(soutt_dk_2))+eps); % 对数化
%     % 参数计算-->IRW
%     [irw_dk_2,~,~] = get_irw(fftshift(soutt_dk_2_nor));
%     IRW2(i) = irw_dk_2;
%     % 参数计算-->PSLR
%     [pslr_dk_2] = get_pslr(fftshift(soutt_dk_2_log));
%     PSLR2(i) = pslr_dk_2;
%     % 参数计算-->ISLR
%     [islr_dk_2] = get_islr(fftshift(soutt_dk_2_nor),5);
%     ISLR2(i) = islr_dk_2;
% % 参数计算-->方式三                                         
%     Soutf_dk_3 = Srf.*Hf_dk_3;
%     soutt_dk_3 = ifft(ifftshift(Soutf_dk_3));          % 方式三匹配滤波结果 
%     soutt_dk_3_nor = abs(soutt_dk_3)./max(abs(soutt_dk_3));               % 归一化
%     soutt_dk_3_log = 20*log10(abs(soutt_dk_3)./max(abs(soutt_dk_3))+eps); % 对数化
%     % 参数计算-->IRW
%     [irw_dk_3,~,~] = get_irw(soutt_dk_3_nor);
%     IRW3(i) = irw_dk_3;
%     % 参数计算-->PSLR
%     [pslr_dk_3] = get_pslr(soutt_dk_3_log);
%     PSLR3(i) = pslr_dk_3;
%     % 参数计算-->ISLR
%     [islr_dk_3] = get_islr(soutt_dk_3_nor,5);
%     ISLR3(i) = islr_dk_3;           
    % 进度内容
    pause(0);
    Time_Trans   = Time_Transform(toc);
    Time_Disp    = Time_Display(Time_Trans);
    Display_Data = num2str(roundn(i/N*100,-1));
    Display_Str  = ['Computation Progress ... ',Display_Data,'%',' --- ',...
                    'Using Time: ',Time_Disp];

    waitbar(i/N,wait_title,Display_Str)
end
pause(1);
close(wait_title);
toc
% % 绘图
% H7 = figure;
% set(H7,'position',[50,50,900,900]);
% subplot(331),plot(QPE/pi,(IRW1-IRW1(1))/IRW1(1)*100,'k')
% title('(a)IRW'),xlabel('|QPE|(\pi弧度)'),ylabel('展宽百分比')
% subplot(332),plot(QPE/pi,PSLR1,'k')
% title('(b)PSLR'),xlabel('|QPE|(\pi弧度)'),ylabel('PSLR/dB')
% subplot(333),plot(QPE/pi,ISLR1,'k')
% title('(c)ISLR'),xlabel('|QPE|(\pi弧度)'),ylabel('ISLR/dB')
% subplot(334),plot(QPE/pi,(IRW2-IRW2(1))/IRW2(1)*100,'k')
% title('(a)IRW'),xlabel('|QPE|(\pi弧度)'),ylabel('展宽百分比')
% subplot(335),plot(QPE/pi,PSLR2,'k')
% title('(b)PSLR'),xlabel('|QPE|(\pi弧度)'),ylabel('PSLR/dB')
% subplot(336),plot(QPE/pi,ISLR2,'k')
% title('(c)ISLR'),xlabel('|QPE|(\pi弧度)'),ylabel('ISLR/dB')
% subplot(337),plot(QPE/pi,(IRW3-IRW3(1))/IRW3(1)*100,'k')
% title('(a)IRW'),xlabel('|QPE|(\pi弧度)'),ylabel('展宽百分比')
% subplot(338),plot(QPE/pi,PSLR3,'k')
% title('(b)PSLR'),xlabel('|QPE|(\pi弧度)'),ylabel('PSLR/dB')
% subplot(339),plot(QPE/pi,ISLR3,'k')
% title('(c)ISLR'),xlabel('|QPE|(\pi弧度)'),ylabel('ISLR/dB')
% suptitle('当\beta=2.5时的IRW、PSLR、ISLR与QPE之间的关系')
%% 提取冲击响应宽度
function [irw,locleft,locright] = get_irw(Af)
    % 找到Af的最大位置
    [~,locmax] = max(Af);
    % 找到locmax左边最接近-3dB的位置
    [~,locleft] = min(abs(Af(1:locmax-1)/max(abs(Af(1:locmax-1)))-0.707));
    % 找到locmax右边最接近-3dB的位置
    [~,locright] = min(abs(Af(locmax+1:end)/max(abs(Af(locmax+1:end)))-0.707));
    locright = locright + locmax;
    % 得到3dB波束宽度
    irw = locright-locleft;
end
%% 提取峰值旁瓣比
function [pslr] = get_pslr(Af)
    % 找到所有的pesks
    peaks = findpeaks(Af);
    % 对peaks进行排序
    peaks = sort(peaks,'descend');
    % 得到第一旁瓣
    pslr = peaks(2);
end
%% 提取积分旁瓣比
function [islr] = get_islr(Af,Nr)
    % 找到Af的最大位置
    [~,locmax] = max(Af);
    % 找到locmax左边最接近-3dB的位置
    [~,locleft] = min(abs(Af(1:locmax-1)/max(abs(Af(1:locmax-1)))-0.707));
    % 找到locmax右边最接近-3dB的位置
    [~,locright] = min(abs(Af(locmax+1:end)/max(abs(Af(locmax+1:end)))-0.707));
    locright = locright + locmax;
    % 计算总功率
    P_total = sum(Af(locleft-Nr:locright+Nr).^2);
    % 计算主瓣功率
    P_main = sum(Af(locleft:locright).^2);
    % 一维积分旁瓣比
    islr = 10*log10((P_total-P_main)./P_main);
end
%% 时间转换函数
function y = Time_Transform(u)
    Time_in = u(1);
    Hours   = fix(Time_in/3600);
    Minutes = fix((Time_in-Hours*3600)/60);
    Seconds = fix(Time_in-Hours*3600-Minutes*60);
    Time_out = [Hours Minutes Seconds];
    y = Time_out;
end
%% 时间显示函数
function y = Time_Display(u)
    Hours   = u(1);
    Minutes = u(2);
    Seconds = u(3);
    
    if Hours == 0
        if Minutes == 0
            Time_out = [num2str(Seconds),'','s'];
        else
            Time_out = [num2str(Minutes),'','m','',...
                            num2str(Seconds),'','s'];
        end 
    else
        Time_out = [num2str(  Hours),'','h','',...
                        num2str(Minutes),'','m','',...
                        num2str(Seconds),'','s'];
    end
    y = Time_out;
end