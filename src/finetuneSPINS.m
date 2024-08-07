function [hpSPINS,wt] = finetuneSPINS(ID,RF_duration,dt, x0)

addpath('src/');
trainID = [ID];
num = length(trainID);
UA = cell(num,1);
Ub = cell(num,1);
testID = ID;

global A maskMS b1arr posarr nKT gamma b0mapMS sysmat phasetrack b off trel myfox wt passband stopband point_P point_S

ID = 1; % Could be changed for UP k-space optimization

path = ['data/s',num2str(trainID(ID)),'/'];
filename = [path,'calibdata_sag'];
load(filename);
foxkt = [0.24,0.306,0.3]; %FOV
poffset = [0 0 0];
soi = 1:1:100;
filename = [path,'mask_SS'];
load(filename);
temp = false(80,102,100);
for ss = 25:10:65
    temp(ss,:,:) = mask_SS(ss,:,:);
end
maskMS = temp;%downsample
b1mapsMSn = rfmap_sag(:,:,soi,:);
b0mapMS = b0map_sag(:,:,soi)/1e6;
gamma = 2.675e8;
[b1arr,posarr] = create_array(b1mapsMSn,maskMS,foxkt,1e-3*poffset);
T = RF_duration;%total time,sec
nKT = T/dt;
dt = 10e-6;%dwell time,sec
len = round(T/dt);
phasetrack = 1:len;
[nspa,nchs] = size(b1arr);

passband = -250:125:250;
stopband = -1300:125:-800;
% passband = [0];
% stopband = [-1050];
%%% For chasing a faster computing time, you could change the passband and
%%% stopband to make a more rough calculation. As we just want to find a
%%% good start point, it is resonable to abandon some precise pursuits.

point_P = length(passband);
point_S = length(stopband);

%%%
A = complex(zeros((point_P+point_S)*nspa,nchs*nKT));
b = [(1/180)*pi*ones(point_P*nspa,1);zeros(point_S*nspa,1)];
gObj = gradPulse(ones(1,phasetrack(end)),dt,'unitgrad','Tx');
kb0 = gObj.KspaceTrajectory;
off = b0mapMS(maskMS)*kb0(phasetrack);
trel = 1/(gamma/2/pi)*kb0(phasetrack);
wt = zeros(nKT*nchs,1);
sysmat = complex(zeros(nspa,nchs*nKT)); % this prealloc results in much faster construct

Xstart = x0;
fun = @(x)FAupdate(x);
nonlcon = @holds;
iter = 100;
options = optimoptions('fmincon','Display','iter','MaxFunEvals',10^18,'TolFun',1e-10,'TolCon',1e-10,'TolX',1e-8,'MaxIter',iter,'Algorithm','active-set');
theta0 = fmincon(fun,x0,[],[],[],[],[],[],nonlcon,options);