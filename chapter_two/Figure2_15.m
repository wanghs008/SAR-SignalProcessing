clc
close all
clear all

N = 500;
t = -8:1/N:8-1/N;

st = sin(pi*t)./(pi*t);

figure
p1 = plot(t,st,'k');hold on
p2 = plot(t,st.*kaiser(length(t),2.5)','r-.');  % kaiser窗加权
title('Kaiser窗加权后的sinc函数，\beta=2.5'),xlabel('时间(采样点)'),ylabel('幅度')
axis([-8 8,-0.4 1.2])
grid on

legend('boxoff');
lgd = legend([p1,p2],{'未经加权','经过加权'},'FontSize',8,'TextColor','k','Location','northeast','NumColumns',1);
% title(lgd,'My Legend Title')