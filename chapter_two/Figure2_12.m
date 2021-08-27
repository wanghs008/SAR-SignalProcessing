clc
close all
clear all

N = 100;
beta = [0,1,2,3,4,5,6];

hw = [];
pslr = [];
% tempn = [];

% set(figure,'position',[100,100,1200,600]);
for i = 1:1:7
    % 生成kauser窗
    A = kaiser(N,beta(i));
    % 绘制Kaiser窗
    % plot(A),hold on
    % 对kaiser窗进行傅里叶变换
    temp = fft(A,2^14);
    temp = 20*log10(abs(temp)./max(abs(temp)));
    % tempn = [tempn,temp];
    % 提取kaiser窗的3dB波束宽度
    [hwn,locleft,locright] = get_hw(temp);
    hw = [hw,hwn];
    % 提取Kaiser窗的峰值旁瓣比
    [pslrn] = get_pslr(temp);
    pslr = [pslr,pslrn];
end

% title('图2.11 不同\beta值的Kaiser窗形状'),xlabel('采样数'),ylabel('幅度')
% axis([-11 110,-0.2 1.2])
% grid on

% text('Interpreter','latex','String','$\beta=0.0$','Position',[-9 1.00],'FontSize',16);
% text('Interpreter','latex','String','$\beta=1.0$','Position',[-9 0.80],'FontSize',16);
% text('Interpreter','latex','String','$\beta=2.0$','Position',[-9 0.45],'FontSize',16);
% text('Interpreter','latex','String','$\beta=3.0$','Position',[-9 0.21],'FontSize',16);
% text('Interpreter','latex','String','$\beta=4.0$','Position',[-9 0.10],'FontSize',16);
% text('Interpreter','latex','String','$\beta=5.0$','Position',[-9 0.05],'FontSize',16);
% text('Interpreter','latex','String','$\beta=6.0$','Position',[-9 0.00],'FontSize',16);

% tempn
% tempn(:,1)
% figure,set(figure,'position',[100,100,1200,600]);
% plot(fftshift(tempn))
% title('图2.11 不同\beta值的Kaiser窗频谱'),xlabel('采样数'),ylabel('幅度/dB')
% axis([5000 11500,-150 10])
% grid on
% set(legend,'Location','NorthEastOutside')
% legend('\beta=0','\beta=1','\beta=2','\beta=3','\beta=4','\beta=5','\beta=6')

% figure,set(figure,'position',[100,100,1200,600]);
% plot(fftshift(tempn(:,1)))
% title('图2.11 \beta=0时的kaiser窗频谱'),xlabel('采样数'),ylabel('幅度/dB')
% axis([5000 11500,-150 10])
% grid on

% figure,set(figure,'position',[100,100,1200,600]);
% plot(fftshift(tempn(:,7)))
% title('图2.11 \beta=6时的kaiser窗频谱'),xlabel('采样数'),ylabel('幅度/dB')
% axis([5000 11500,-150 10])
% grid on

set(figure,'position',[100,100,1200,600]);
% 3dB宽度展宽比
% hw
subplot(121),plot((hw./hw(1)-1)*100)
title('(a)3dB宽度展宽比'),xlabel('Kaiser窗\beta'),ylabel('展宽比')
grid on
% 峰值旁瓣比
% pslr
subplot(122),plot(pslr)
title('(b)峰值旁瓣比'),xlabel('Kaiser窗\beta'),ylabel('PSLR/dB')
grid on
suptitle('图2.12 不同Kaiser窗的展宽和峰值旁瓣比')

function [hw,locleft,locright] = get_hw(Af)
    % 找到Af的最大位置
    [~,locmax] = max(Af);
    % 找到locmax左边最接近-3dB的位置
    [~,locleft] = min(abs(Af(1:locmax)+3));
    % 找到locmax右边最接近-3dB的位置
    [~,locright] = min(abs(Af(locmax:end)+3));
    locright = locright + locmax - 1;
    % 得到3dB波束宽度
    hw = locright-locleft;
end

function [pslr] = get_pslr(Af)
    % 找到所有的pesks
    peaks = findpeaks(Af);
    % 对peaks进行排序
    peaks = sort(peaks,'descend');
    % 得到第一旁瓣
    pslr = peaks(2);
end