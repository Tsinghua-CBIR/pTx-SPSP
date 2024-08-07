function [] = showPulse(rf,grad,T,dt)
figure;
plot((0:dt:T)*1e3,abs(rf));
set(gca,'FontSize',12);
ax = gca;
ax.LineWidth = 2; % 坐标轴线条磅数
ax.FontName = 'Times New Roman'; % 字体格式
ax.FontSize = 14; % 字体大小

xlabel('T/ms', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('rf/|V|', 'FontName', 'Times New Roman', 'FontSize', 14);
grid on;

figure;
hold on
for i = 1:3
plot((0:dt:T)*1e3,(grad(i,:)));
end
ax = gca;
ax.LineWidth = 2; % 坐标轴线条磅数
ax.FontName = 'Times New Roman'; % 字体格式
ax.FontSize = 14; % 字体大小
legend('Gx','Gy','Gz','Location','Northeast');
set(gca,'Fontname','Times New Roman','FontSize',16);
xlabel('T/ms', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('grad/mTm^{-1}', 'FontName', 'Times New Roman', 'FontSize', 14);
end