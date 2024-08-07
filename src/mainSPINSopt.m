function [hpSPINS,wt] = mainSPINSopt(ID,RF_duration,X0, dt, TR, RFA)

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
nKT = 24;
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
wt = zeros(384/2,1);
sysmat = complex(zeros(nspa,nchs*nKT)); 
%%% The prealloc results in much faster construct

%% CMA-ES
OPTS=cmaes; 
% OPTS.DiagonalOnly = 1;
% OPTS.CMA.active = 2;
OPTS.StopFitness=1e-5;
OPTS.StopIter = 50;
OPTS.StopFunEvals = 6e-2;
OPTS.LogFilenamePrefix = 'process/outcmaes';
OPTS.SaveFilename = 'process/variablescmaes.mat';
% opts.DispModulo = 10;
% OPTS.PopSize = 32;
% XP = cmaes('fun',enlargeX0,50,OPTS);
tp = cmaes('foptSPINS',X0,10,OPTS);
% tp = cmaes('foptSPINS',[20,8/2.88*6,2/2.88*30,20,20],10,OPTS);
addpath('process/');
load('process/variablescmaes.mat');
hpSPINS = xmin;


%% Get Start-point and execute ADMM
wt0 = repelem(wt,12,1);
nKT = round(T/dt);
A = complex(zeros((point_P+point_S)*nspa,nchs*nKT));
b = [(1/180)*pi*ones(point_P*nspa,1);zeros(point_S*nspa,1)];
gObj = gradPulse(ones(1,phasetrack(end)),dt,'unitgrad','Tx');
kb0 = gObj.KspaceTrajectory;
off = b0mapMS(maskMS)*kb0(phasetrack);
trel = 1/(gamma/2/pi)*kb0(phasetrack);
sysmat = complex(zeros(nspa,nchs*nKT)); 
%%% The prealloc results in much faster construct
hpSPINS = X0;%%%just for test, you can delete it.
foptSPINSos(hpSPINS,wt0);
%%% Here you could optimize hpSPINS together, but as it will slow down the
%%% caculation, we do not recommend you to do that in the demo.

vx0 = [real(wt);imag(wt)];
disp('ADMM, please wait for a while');
while (1) 
    %%% local SAR update
    vx2 = vx0;
    fun = @(x)SARupdate(x,vx2);
    nonlcon = @(x)SARcon(x,TR,RFA,dt);
    iter = 10;
    options = optimoptions('fmincon','MaxFunEvals',10^18,'TolFun',1e-10,'TolCon',1e-10,'TolX',1e-8,'MaxIter',iter,'Algorithm','sqp');
    vx1 = fmincon(fun,vx2,[],[],[],[],[],[],nonlcon,options);

    %%% FA update
    fun = @(x)FAIupdate(x,vx1);
    nonlcon = @holds;
    iter = 10;
    options = optimoptions('fmincon','MaxFunEvals',10^18,'TolFun',1e-10,'TolCon',1e-10,'TolX',1e-8,'MaxIter',iter,'Algorithm','sqp');
    vx2 = fmincon(fun,vx1,[],[],[],[],[],[],nonlcon,options);
    lenv = round(length(vx2)/2);
    wt = vx2(1:lenv)+1i*vx2(lenv+1:end);
    %%% we've got a good startpoint, so we don't want to go far away from it.
    %%% Thus we used a large weight on |vx1-vx2| to keep the parameters
    %%% from changing too far
    %%% 
    if sum(abs(vx2-vx1))<1e-7 %%% quit loop 
        break
    end
end
end