function [hpKT,wt] = mainKTdopt(ID,RF_duration,X0, dt, TR, RFA)

addpath('src');
trainID = ID;
num = length(trainID);
UA = cell(num,1);
Ub = cell(num,1);
testID = ID;
global A maskMS b1arr posarr nKT gamma b0mapMS phasetrack b off trel sysmat myfox wt passband stopband point_P point_S

ID = 1; % 可以改成for循环用作UP训练

path = ['data/s',num2str(trainID(ID)),'/'];
filename = [path,'calibdata_sag'];
load(filename);
foxkt = [0.24,0.306,0.3]; %fov的计算，后期根据实际激发的区域来调整
poffset = [0 0 0];
soi = 1:1:100;
filename = [path,'mask_SS'];
load(filename);
temp = false(80,102,100);
for ss = 25:10:65
    temp(ss,:,:) = mask_SS(ss,:,:);
end
maskMS = temp;%进一步执行降采样，降低计算量
b1mapsMSn = rfmap_sag(:,:,soi,:);
b0mapMS = b0map_sag(:,:,soi)/1e6;
gamma = 2.675e8;
[b1arr,posarr] = create_array(b1mapsMSn,maskMS,foxkt,1e-3*poffset);
T = RF_duration;%单位s，设置波形总时间
nKT = 24;
dt = 10e-6;%dwell time，单位s
len = round(T/dt);

ngrps = 8;
mytheta = 0:180/ngrps:180;
mytheta = mytheta(1:end-1);
myphi = mytheta;

kp0 = [0 0 0;
       -1 0 0;
      1 0 0].';
kp= [];
for ind=1:length(mytheta)
    ikp = angle2dcm(deg2rad(mytheta(ind)),deg2rad(myphi(ind)),0,'ZYX') * kp0;
    kp = [kp ikp];
end
myfox = mean(foxkt);
kp = 2*pi./myfox.* kp;




[nspa,nchs] = size(b1arr);
sysmat = complex(zeros(nspa,nchs*nKT)); % this prealloc results in much faster construct

passband = -250:125:250;
stopband = -1300:125:-800;
% passband = [0];
% stopband = [-1050];

point_P = length(passband);
point_S = length(stopband);

%%%Prr
A = complex(zeros((point_P+point_S)*nspa,nchs*nKT));
b = [(1/180)*pi*ones(point_P*nspa,1);zeros(point_S*nspa,1)];

wt = zeros(384/2,1);
sysmat = complex(zeros(nspa,nchs*nKT));





temp = zeros(3,8);
for i = 1:8
    temp(:,i) = kp(:,i*3);
end
kp_ini = reshape(temp,[24,1]);
kp_ini = X0;

OPTS=cmaes; 
% OPTS.DiagonalOnly = 0;
OPTS.CMA.active = 2;
OPTS.StopFitness=1e-5;
OPTS.StopIter = 1000;
OPTS.StopFunEvals = 1e-1;
OPTS.LogFilenamePrefix = 'process/outcmaes';
OPTS.SaveFilename = 'process/variablescmaes.mat';
% opts.DispModulo = 10;
% OPTS.PopSize = 20;
% XP = cmaes('fun',enlargeX0,50,OPTS);
hpKT = cmaes('foptKT',kp_ini,10,OPTS);

%% Get Start-point and execute ADMM
hpKT = X0;%%%just for test, you can delete it.
foptKT(hpKT);
vx2 = [real(wt);imag(wt)];
disp(' ');
disp('ADMM, please wait for a while');
while (1) 
    %%% local SAR update
    vx0 = vx2;
    vx1 = localSARupKT(vx0,TR,RFA,dt);

    %%% FA update
    fun = @(x)FAIupdate(x,vx1);
    nonlcon = @holds;
    iter = 10;
    options = optimoptions('fmincon','MaxFunEvals',10^18,'TolFun',1e-10,'TolCon',1e-10,'TolX',1e-8,'MaxIter',iter,'Algorithm','sqp');
    vx2 = fmincon(fun,vx1,[],[],[],[],[],[],nonlcon,options);
    lenv = round(length(vx2)/2);
    wt = vx2(1:lenv)+1i*vx2(lenv+1:end);
    
    tempx = wt;
    fun = @(theta)foptKTfixX(theta,tempx);
    nonlcon = @holds;
    iter = 1;
    options = optimoptions('fmincon','MaxFunEvals',10^18,'TolFun',1e-10,'TolCon',1e-10,'TolX',1e-8,'MaxIter',iter,'Algorithm','sqp');
    hpKT = fmincon(fun,hpKT,[],[],[],[],[],[],nonlcon,options);
    %%% we've got a good startpoint, so we don't want to go far away from it.
    %%% Thus we used a large weight on |vx1-vx2| to keep the parameters
    %%% from changing too far
    %%% 
    if sum(abs(vx2-vx1))^2<1e-5 %%% quit loop 
        break
    end
end


kp = zeros(3,24);

for group = 1:8
    kp(:,group*3-2) = zeros(3,1);
    kp(:,group*3-1) = -hpKT(group*3-2:group*3);
    kp(:,group*3) = hpKT(group*3-2:group*3);
end
hpKT = kp;
end