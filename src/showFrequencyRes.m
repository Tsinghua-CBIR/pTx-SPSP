function [] = showFrequencyRes(ID,rf,grad,dt,RFA)
% We want to offer a template which is compatible with UP design.
% So here you can fill in the ID of calibration, forming testing-set.
% default to trainID(1) if num == 1.
myrf = rf*1e-6;
mygrad = grad*1e-3;
disp('Validation')
filename = ['data/s',num2str(ID),'/calibdata_sag'];
load(filename);
soi = 30:5:70;%%check for the direction
foxkt = [0.24,0.306,0.12];%% pay attention to the FOV! With fewer slices, the FOV would be different.
gamma = 2.675e8;
poffset = [0 0 0]; %offset
maskMS = mask_sag(:,:,soi);
b1mapsMSn = rfmap_sag(:,:,soi,:);
b0mapMS = b0map_sag(:,:,soi)/1e6; 


TFA = RFA;
sf = -1200;
ef = 200;
freq = sf:100:ef;
af = zeros(1,length(freq));
RMSEfat = zeros(1,length(freq));
RMSEwater = zeros(1,length(freq));
COV = zeros(1,length(freq));
for i = 1:length(freq)
    mxypat = run_bloch_sim (myrf,mygrad,b1mapsMSn,maskMS,foxkt,b0mapMS-freq(i)/(gamma/2/pi),0,[],dt,poffset);
    vec = abs(mxypat(maskMS));
    tar = sin(TFA/180*pi);
    af(i) = mean(vec/tar);
    RMSEfat(i) = sqrt(sum((vec-0).^2)/length(vec));
    RMSEwater(i) = sqrt(sum((vec-tar).^2)/length(vec));
    COV(i) = std(vec)/mean(vec);
end

figure;

plot(sf:100:ef,af,'o-', 'LineWidth', 2);  


% 设置坐标轴属性
ax = gca;
ax.LineWidth = 1.5; % 坐标轴线条磅数
ax.FontName = 'Times New Roman'; % 字体格式
ax.FontSize = 14; % 字体大小

% 设置x轴标签
xlabel('Frequency offset (Hz)', 'FontName', 'Times New Roman', 'FontSize', 14);

% 设置y轴标签为空
ylabel('', 'FontName', 'Times New Roman', 'FontSize', 14);


% 添加网格线
grid on;
ax.GridLineStyle = '--'; % 网格线样式
ax.GridAlpha = 0.5; % 网格线透明度

% 调整图形边距
xlim([sf,ef]);
ylim([0,1.1]);
set(gca,'FontSize',14);
end