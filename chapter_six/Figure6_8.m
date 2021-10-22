clc
clear
close all
%% 参数设置
%  已知参数--》距离向参数
R_eta_c = 20e+3;                % 景中心斜距
Tr = 2.5e-6;                    % 发射脉冲时宽
Kr = 20e+12;                    % 距离向调频率
alpha_os_r = 1.2;               % 距离过采样率
Nrg = 320;                      % 距离线采样点数
%  计算参数--》距离向参数
Bw = abs(Kr)*Tr;                % 距离信号带宽
Fr = alpha_os_r*Bw;             % 距离向采样率
Nr = round(Fr*Tr);              % 距离采样点数(脉冲序列长度)
%  已知参数--》方位向参数
c = 3e+8;                       % 电磁传播速度
Vr = 150;                       % 等效雷达速度
Vs = Vr;                        % 卫星平台速度
Vg = Vr;                        % 波束扫描速度
f0 = 5.3e+9;                    % 雷达工作频率
Delta_f_dop = 80;               % 多普勒带宽
alpha_os_a = 1.25;              % 方位过采样率
Naz = 256;                      % 距离线数
theta_r_c = 1.6*pi/180;           % 波束斜视角
%  计算参数--》方位向参数
lambda = c/f0;                  % 雷达工作波长
t_eta_c = -R_eta_c*sin(theta_r_c)/Vr;
                                % 波束中心穿越时刻
f_eta_c = 2*Vr*sin(theta_r_c)/lambda;
                                % 多普勒中心频率
La = 0.886*2*Vs*cos(theta_r_c)/Delta_f_dop;               
                                % 实际天线长度
Fa = alpha_os_a*Delta_f_dop;    % 方位向采样率
Ta = 0.886*lambda*R_eta_c/(La*Vg*cos(theta_r_c));
                                % 目标照射时间
R0 = R_eta_c*cos(theta_r_c);    % 最短斜距
Ka = 2*Vr^2*cos(theta_r_c)^2/lambda/R0;              
                                % 方位向调频率
theta_bw = 0.886*lambda/La;     % 方位向3dB波束宽度
theta_syn = Vs/Vg*theta_bw;     % 合成角宽度
Ls = R_eta_c*theta_syn;         % 合成孔径长度
%  参数计算
rho_r = c/(2*Fr);               % 距离向分辨率
rho_a = La/2;                   % 距离向分辨率
Trg = Nrg/Fr;                   % 发射脉冲时宽
Taz = Naz/Fa;                   % 目标照射时间
d_t_tau = 1/Fr;                 % 距离采样时间间隔
d_t_eta = 1/Fa;                 % 方位采样时间间隔
d_f_tau = Fr/Nrg;               % 距离采样频率间隔    
d_f_eta = Fa/Naz;               % 方位采样频率间隔
%% 目标设置
%  设置目标点相对于景中心之间的距离
A_r =   0; A_a =   0;                                   % A点位置
B_r = -50; B_a = -50;                                   % B点位置
C_r = -50; C_a = +50;                                   % C点位置
D_r = +50; D_a = C_a + (D_r-C_r)*tan(theta_r_c);        % D点位置
%  得到目标点相对于景中心的位置坐标
A_x = R0 + A_r; A_Y = A_a;                              % A点坐标
B_x = R0 + B_r; B_Y = B_a;                              % B点坐标
C_x = R0 + C_r; C_Y = C_a;                              % C点坐标
D_x = R0 + D_r; D_Y = D_a;                              % D点坐标
NPosition = [A_x,A_Y;
             B_x,B_Y;
             C_x,C_Y;
             D_x,D_Y;];                                 % 设置数组
fprintf( 'A点坐标为[%+3.3f，%+3.3f]km\n', NPosition(1,1)/1e3, NPosition(1,2)/1e3 );
fprintf( 'B点坐标为[%+3.3f，%+3.3f]km\n', NPosition(2,1)/1e3, NPosition(2,2)/1e3 );
fprintf( 'C点坐标为[%+3.3f，%+3.3f]km\n', NPosition(3,1)/1e3, NPosition(3,2)/1e3 );
fprintf( 'D点坐标为[%+3.3f，%+3.3f]km\n', NPosition(4,1)/1e3, NPosition(4,2)/1e3 );
%  得到目标点的波束中心穿越时刻
Ntarget = 4;
Tar_t_eta_c = zeros(1,Ntarget);
for i = 1 : Ntarget
    DeltaX = NPosition(i,2) - NPosition(i,1)*tan(theta_r_c);
    Tar_t_eta_c(i) = DeltaX/Vs;
end
%  得到目标点的绝对零多普勒时刻
Tar_t_eta_0 = zeros(1,Ntarget);
for i = 1 : Ntarget
    Tar_t_eta_0(i) = NPosition(i,2)/Vr;
end
%% 变量设置
%  时间变量 以景中心的零多普勒时刻作为方位向零点
t_tau = (-Trg/2:d_t_tau:Trg/2-d_t_tau) + 2*R_eta_c/c;   % 距离时间变量
t_eta = (-Taz/2:d_t_eta:Taz/2-d_t_eta) + t_eta_c;       % 方位时间变量
%  长度变量
r_tau = (t_tau*c/2)*cos(theta_r_c);                     % 距离长度变量
%  频率变量 
f_tau = fftshift(-Fr/2:d_f_tau:Fr/2-d_f_tau);           % 距离频率变量
f_tau = f_tau - round((f_tau-0)/Fr)*Fr;                 % 距离频率变量(可观测频率)  
f_eta = fftshift(-Fa/2:d_f_eta:Fa/2-d_f_eta);           % 方位频率变量
f_eta = f_eta - round((f_eta-f_eta_c)/Fa)*Fa;           % 方位频率变量(可观测频率)
%% 坐标设置     
%  以距离时间为X轴，方位时间为Y轴
[t_tauX,t_etaY] = meshgrid(t_tau,t_eta);                % 设置距离时域-方位时域二维网络坐标
%  以距离长度为X轴，方位频率为Y轴                                                                                                            
[r_tauX,f_etaY] = meshgrid(r_tau,f_eta);                % 设置距离时域-方位频域二维网络坐标
%  以距离频率为X轴，方位频率为Y轴                                                                                                            
[f_tau_X,f_eta_Y] = meshgrid(f_tau,f_eta);              % 设置频率时域-方位频域二维网络坐标
%% 信号设置--》原始回波信号                                                        
tic
wait_title = waitbar(0,'开始生成雷达原始回波数据 ...');  
pause(1);                                                     
st_tt = zeros(Naz,Nrg);
for i = 1 : 1
    %  计算目标点的瞬时斜距
    R_eta = sqrt( NPosition(i,1)^2 +...
                  Vr^2*(t_etaY-Tar_t_eta_0(i)).^2 ); 
    %{
    R_eta = NPosition(i,1) + Vr^2*t_etaY.^2/(2*NPosition(i,1));   
    %}
    %  后向散射系数幅度
    A0 = [1,1,1,1]*exp(+1j*0);   
    %  距离向包络
    wr = (abs(t_tauX-2*R_eta/c) <= Tr/2);                               
    %  方位向包络
    wa = sinc(0.886*atan(Vg*(t_etaY-Tar_t_eta_c(i))/NPosition(i,1))/theta_bw).^2;      
    %  接收信号叠加
    st_tt_tar = A0(i)*wr.*wa.*exp(-1j*4*pi*f0*R_eta/c)...
                            .*exp(+1j*pi*Kr*(t_tauX-2*R_eta/c).^2); 
    %{
    st_tt_tar = A0(i)*wr.*wa.*exp(-1j*4*pi*R0/lambda)...
                            .*exp(-1j*pi*Ka*t_eta_Y.^2)...
                            .*exp(+1j*pi*Kr*(t_tauX-2*R_eta/c).^2);
    %}                                                          
    st_tt = st_tt + st_tt_tar;  
    
    pause(0.001);
    Time_Trans   = Time_Transform(toc);
    Time_Disp    = Time_Display(Time_Trans);
    Display_Data = num2str(roundn(i/Ntarget*100,-1));
    Display_Str  = ['Computation Progress ... ',Display_Data,'%',' --- ',...
                    'Using Time: ',Time_Disp];
    waitbar(i/Ntarget,wait_title,Display_Str)
end
pause(1);
close(wait_title);
toc
%% 信号设置--》一次距离压缩                  
%  信号变换-->方式三：根据脉冲频谱特性直接在频域生成频域匹配滤波器
%  加窗函数
window = kaiser(Nrg,2.5)';                              % 时域窗
Window = fftshift(window);                          	% 频域窗
Hrf = (abs(f_tau_X)<=Bw/2).*Window.*exp(+1j*pi*f_tau_X.^2/Kr);  
%  匹配滤波
Sf_ft = fft(st_tt,Nrg,2);
Srf_tf = Sf_ft.*Hrf;
srt_tt = ifft(Srf_tf,Nrg,2);
%% 信号设置--》方位向傅里叶变换
Saf_tf = fft(srt_tt,Naz,1);
%% 信号设置--》距离徙动校正
RCM = lambda^2*r_tauX.*f_etaY.^2/(8*Vr^2);              % 需要校正的距离徙动量
RCM = R0 + RCM - R_eta_c;                               % 将距离徙动量转换到原图像坐标系中
offset = RCM/rho_r;                                     % 将距离徙动量转换为距离单元偏移量
%  计算插值系数表(八点插值)
x_tmp = repmat(-4:3,[16,1]);                            % 插值长度                          
x_tmp = x_tmp + repmat(((1:16)/16).',[1,8]);            % 量化位移
hx = sinc(x_tmp);                                       % 生成插值核
kwin = repmat(kaiser(8,2.5).',[16,1]);                  % 加窗
hx = kwin.*hx;
hx = hx./repmat(sum(hx,2),[1,8]);                       % 核的归一化
%  插值表校正
tic
wait_title = waitbar(0,'开始进行距离徙动校正(八点插值) ...');  
pause(1);
Srcmf_tf_8 = zeros(Naz,Nrg);
for a_tmp = 1 : Naz
    for r_tmp = 1 : Nrg
        offset_ceil = ceil(offset(a_tmp,r_tmp));
        offset_frac = round((offset_ceil - offset(a_tmp,r_tmp)) * 16);
        if offset_frac == 0
           Srcmf_tf_8(a_tmp,r_tmp) = Saf_tf(a_tmp,ceil(mod(r_tmp+offset_ceil-0.1,Nrg))); 
        else
           Srcmf_tf_8(a_tmp,r_tmp) = Saf_tf(a_tmp,ceil(mod((r_tmp+offset_ceil-4:r_tmp+offset_ceil+3)-0.1,Nrg)))*hx(offset_frac,:).';
        end
    end
    
    pause(0.001);
    Time_Trans   = Time_Transform(toc);
    Time_Disp    = Time_Display(Time_Trans);
    Display_Data = num2str(roundn(a_tmp/Naz*100,-1));
    Display_Str  = ['Computation Progress ... ',Display_Data,'%',' --- ',...
                    'Using Time: ',Time_Disp];
    waitbar(a_tmp/Naz,wait_title,Display_Str)
end
pause(1);
close(wait_title);
toc
%% 信号设置--》距离徙动校正
RCM = lambda^2*r_tauX.*f_etaY.^2/(8*Vr^2);              % 需要校正的距离徙动量
RCM = R0 + RCM - R_eta_c;                               % 将距离徙动量转换到原图像坐标系中
offset = RCM/rho_r;                                     % 将距离徙动量转换为距离单元偏移量
%  计算插值系数表(最近邻域插值)
x_tmp = repmat(-4:3,[16,1]);                            % 插值长度                          
x_tmp = x_tmp + repmat(((1:16)/16).',[1,8]);            % 量化位移
hx = sinc(x_tmp);                                       % 生成插值核
kwin = repmat(kaiser(8,2.5).',[16,1]);                  % 加窗
hx = kwin.*hx;
hx = hx./repmat(sum(hx,2),[1,8]);                       % 核的归一化
%  插值表校正
tic
wait_title = waitbar(0,'开始进行距离徙动校正(最近邻域插值) ...');  
pause(1);
Srcmf_tf_0 = zeros(Naz,Nrg);
for a_tmp = 1 : Naz
    for r_tmp = 1 : Nrg
        offset_ceil = ceil(offset(a_tmp,r_tmp));
        offset_frac = round((offset_ceil - offset(a_tmp,r_tmp)) * 16);
        if offset_frac == 0
           Srcmf_tf_0(a_tmp,r_tmp) = Saf_tf(a_tmp,ceil(mod(r_tmp+offset_ceil-0.1,Nrg))); 
        else
           Srcmf_tf_0(a_tmp,r_tmp) = Saf_tf(a_tmp,ceil(mod((r_tmp+offset_ceil-4:r_tmp+offset_ceil+3)-0.1,Nrg)))*hx(offset_frac,:).';
           % 最近邻域
           point = r_tmp + ceil(offset(a_tmp,r_tmp));
           if (point < 1)
               point = 1;
           elseif(point > Nrg)
               point = Nrg;
           end
               Srcmf_tf_0(a_tmp,r_tmp) = Saf_tf(a_tmp,point);
        end
    end
    
    pause(0.001);
    Time_Trans   = Time_Transform(toc);
    Time_Disp    = Time_Display(Time_Trans);
    Display_Data = num2str(roundn(a_tmp/Naz*100,-1));
    Display_Str  = ['Computation Progress ... ',Display_Data,'%',' --- ',...
                    'Using Time: ',Time_Disp];
    waitbar(a_tmp/Naz,wait_title,Display_Str)
end
pause(1);
close(wait_title);
toc
%% 信号设置--》距离徙动校正
RCM = lambda^2*r_tauX.*f_etaY.^2/(8*Vr^2);              % 需要校正的距离徙动量
RCM = R0 + RCM - R_eta_c;                               % 将距离徙动量转换到原图像坐标系中
offset = RCM/rho_r;                                     % 将距离徙动量转换为距离单元偏移量
%  计算插值系数表(四点插值)
x_tmp = repmat(-2:1,[16,1]);                            % 插值长度                          
x_tmp = x_tmp + repmat(((1:16)/16).',[1,4]);            % 量化位移
hx = sinc(x_tmp);                                       % 生成插值核
kwin = repmat(kaiser(4,2.5).',[16,1]);                  % 加窗
hx = kwin.*hx;
hx = hx./repmat(sum(hx,2),[1,4]);                       % 核的归一化
%  插值表校正
tic
wait_title = waitbar(0,'开始进行距离徙动校正(四点插值) ...');  
pause(1);
Srcmf_tf_4 = zeros(Naz,Nrg);
for a_tmp = 1 : Naz
    for r_tmp = 1 : Nrg
        offset_ceil = ceil(offset(a_tmp,r_tmp));
        offset_frac = round((offset_ceil - offset(a_tmp,r_tmp)) * 16);
        if offset_frac == 0
           Srcmf_tf_4(a_tmp,r_tmp) = Saf_tf(a_tmp,ceil(mod(r_tmp+offset_ceil-0.1,Nrg))); 
        else
           Srcmf_tf_4(a_tmp,r_tmp) = Saf_tf(a_tmp,ceil(mod((r_tmp+offset_ceil-2:r_tmp+offset_ceil+1)-0.1,Nrg)))*hx(offset_frac,:).';
        end
    end
    
    pause(0.001);
    Time_Trans   = Time_Transform(toc);
    Time_Disp    = Time_Display(Time_Trans);
    Display_Data = num2str(roundn(a_tmp/Naz*100,-1));
    Display_Str  = ['Computation Progress ... ',Display_Data,'%',' --- ',...
                    'Using Time: ',Time_Disp];
    waitbar(a_tmp/Naz,wait_title,Display_Str)
end
pause(1);
close(wait_title);
toc
%% 信号设置--》方位压缩
%  变量设置
Ka = 2*Vr^2*cos(theta_r_c)^2./(lambda*r_tauX); 
%  计算滤波器
Haf = exp(-1j*pi*f_etaY.^2./Ka);
Haf_offset = exp(-1j*2*pi*f_etaY.*t_eta_c);
%  匹配滤波
Soutf_tf_8 = Srcmf_tf_8.*Haf.*Haf_offset;
soutt_tt_8 = ifft(Soutf_tf_8,Naz,1);
Soutf_tf_0 = Srcmf_tf_0.*Haf.*Haf_offset;
soutt_tt_0 = ifft(Soutf_tf_0,Naz,1);
Soutf_tf_4 = Srcmf_tf_4.*Haf.*Haf_offset;
soutt_tt_4 = ifft(Soutf_tf_4,Naz,1);
% %% 绘图
% H1 = figure();
% set(H1,'position',[100,100,600,300]); 
% subplot(121),imagesc(real(soutt_tt_8))
% %  axis([0 Naz,0 Nrg])
% xlabel('距离时间(采样点)→'),ylabel('←方位时间(采样点)'),title('(a)实部')
% subplot(122),imagesc( abs(soutt_tt_8))
% %  axis([0 Naz,0 Nrg])
% xlabel('距离时间(采样点)→'),ylabel('←方位时间(采样点)'),title('(b)幅度')
% sgtitle('图6.12 方位压缩后的仿真结果','Fontsize',16,'color','k')
% %% 绘图
% H2 = figure();
% set(H2,'position',[100,100,600,300]); 
% subplot(121),imagesc(real(soutt_tt_0))
% %  axis([0 Naz,0 Nrg])
% xlabel('距离时间(采样点)→'),ylabel('←方位时间(采样点)'),title('(a)实部')
% subplot(122),imagesc( abs(soutt_tt_0))
% %  axis([0 Naz,0 Nrg])
% xlabel('距离时间(采样点)→'),ylabel('←方位时间(采样点)'),title('(b)幅度')
% sgtitle('图6.12 方位压缩后的仿真结果','Fontsize',16,'color','k')
% %% 绘图
% H3 = figure();
% set(H3,'position',[100,100,600,300]); 
% subplot(121),imagesc(real(soutt_tt_4))
% %  axis([0 Naz,0 Nrg])
% xlabel('距离时间(采样点)→'),ylabel('←方位时间(采样点)'),title('(a)实部')
% subplot(122),imagesc( abs(soutt_tt_4))
% %  axis([0 Naz,0 Nrg])
% xlabel('距离时间(采样点)→'),ylabel('←方位时间(采样点)'),title('(b)幅度')
% sgtitle('图6.12 方位压缩后的仿真结果','Fontsize',16,'color','k')
%% 信号设置--》点目标分析
%  8点插值
srcmt_tt_8 = ifft(Srcmf_tf_8,Naz,1);
Arcm_8 = srcmt_tt_8(:,round(Nrg/2+1+2*(NPosition(1,1)-R0)/c*Fr));
Arcm_8 = abs(Arcm_8)/max(abs(Arcm_8));
%  最近邻域插值
srcmt_tt_0 = ifft(Srcmf_tf_0,Naz,1);
Arcm_0 = srcmt_tt_0(:,round(Nrg/2+1+2*(NPosition(1,1)-R0)/c*Fr));
Arcm_0 = abs(Arcm_0)/max(abs(Arcm_0));
%  4点插值
srcmt_tt_4 = ifft(Srcmf_tf_4,Naz,1);
Arcm_4 = srcmt_tt_4(:,round(Nrg/2+1+2*(NPosition(1,1)-R0)/c*Fr));
Arcm_4 = abs(Arcm_4)/max(abs(Arcm_4));
%  方位切片
%  8点插值
len_az = 16;
cut_az = -len_az/2:len_az/2-1;
out_az_8 = soutt_tt_8(round(Naz/2+1+NPosition(1,2)/Vr*Fa)+cut_az,...
                      round(Nrg/2+1+2*(NPosition(1,1)-R0)/c*Fr)).';
spec_az_8 = fft(out_az_8);
spec_up_az_8 = ifft(spec_az_8,len_az*16);
spec_up_az_8 = 20*log10(abs(spec_up_az_8)/max(abs(spec_up_az_8)));
%  方位切片
%  最近邻域插值
len_az = 16;
cut_az = -len_az/2:len_az/2-1;

out_az_0 = soutt_tt_0(round(Naz/2+1+NPosition(1,2)/Vr*Fa)+cut_az,...
                      round(Nrg/2+1+2*(NPosition(1,1)-R0)/c*Fr)).';
spec_az_0 = fft(out_az_0);
spec_up_az_0 = ifft(spec_az_0,len_az*16);
spec_up_az_0 = 20*log10(abs(spec_up_az_0)/max(abs(spec_up_az_0)));
%  方位切片
%  4点插值
len_az = 16;
cut_az = -len_az/2:len_az/2-1;

out_az_4 = soutt_tt_4(round(Naz/2+1+NPosition(1,2)/Vr*Fa)+cut_az,...
                      round(Nrg/2+1+2*(NPosition(1,1)-R0)/c*Fr)).';
spec_az_4 = fft(out_az_4);
spec_up_az_4 = ifft(spec_az_4,len_az*16);
spec_up_az_4 = 20*log10(abs(spec_up_az_4)/max(abs(spec_up_az_4)));
%  绘图
H4 = figure();
set(H4,'position',[100,100,900,600]); 
subplot(231),plot(abs(wa(:,1))),axis([0 Naz,0 1])
xlabel('方位向(采样点)'),ylabel('幅度'),title('(a)理想的距离徙动校正(8点插值)')
subplot(232),plot(Arcm_0),axis([0 Naz,0 1])
xlabel('方位向(采样点)'),ylabel('幅度'),title('(b)最近邻域插值的距离徙动校正')
subplot(233),plot(Arcm_4),axis([0 Naz,0 1])
xlabel('方位向(采样点)'),ylabel('幅度'),title('(c)4点插值的距离徙动校正')
subplot(234),plot(spec_up_az_8),axis([0 Naz,-30 0])
xlabel('方位向(采样点)'),ylabel('幅度(dB)')
subplot(235),plot(spec_up_az_0),axis([0 Naz,-30 0])
xlabel('方位向(采样点)'),ylabel('幅度(dB)')
subplot(236),plot(spec_up_az_4),axis([0 Naz,-30 0])
xlabel('方位向(采样点)'),ylabel('幅度(dB)')
sgtitle('图6.8 距离徙动校正不精确时，由调制引入的成对回波','Fontsize',16,'color','k')