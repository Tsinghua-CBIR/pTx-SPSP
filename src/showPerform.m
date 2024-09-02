function [] = showPerform(ID,rf,grad,dt,offset,RFA)
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

mxypat = run_bloch_sim (myrf,mygrad,b1mapsMSn,maskMS,foxkt,b0mapMS-(0+offset)/(gamma/2/pi),...
    0,[],dt,poffset);

vec = asin(abs(mxypat(maskMS)));
Cov = std(vec)/mean(vec); %%% COV output here
disp('%--------------------------------%');
fprintf('Bloch simulated CoV:%d \n',Cov);
disp('%--------------------------------%');


Img = cell(3,3);
for i = 1:3
    for j = 1:3
        temp = abs(mxypat(:,:,i*3-3+j));
        Img{i,j} = asin(temp)/pi*180;
    end
end
Img = cell2mat(Img);

figure;
imagesc(Img);
% maxm = max(max(abs(Img)));
caxis([0,ceil(RFA/10)*10]);
t = 0:0:0;
set(gca,'xtick',t);
set(gca,'ytick',t);
colorbar;
colormap hot;
axis image;

mxypat = run_bloch_sim (myrf,mygrad,b1mapsMSn,maskMS,foxkt,b0mapMS+(1050+offset)/(gamma/2/pi),...
    0,[],dt,poffset);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PLUS!!! left hand system (xyz2PRS) for Siemens!


Img = cell(3,3);
for i = 1:3
    for j = 1:3
        temp = abs(mxypat(:,:,i*3-3+j));
        Img{i,j} = asin(temp)/pi*180;
    end
end
Img = cell2mat(Img);
% 
figure;
imagesc(Img);
% maxm = max(max(abs(Img)));
caxis([0,5]);
t = 0:0:0;
set(gca,'xtick',t);
set(gca,'ytick',t);
colorbar;
colormap hot;
axis image;

end