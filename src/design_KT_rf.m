function [rf,grad,localSAR] = design_KT_rf(ID, RF_duration, dt, TR, RFA)
% rf_design_function - MATLAB function for designing a non-selective pTx-SPSP KT-points RF pulse
%
% Inputs:
%   ID (array) - serial number of calibration, set here so the function
%   could be easily extended to UP training.
%   RF_duration (double) - Duration of the RF pulse, default to 2.88e-3s
%   dt (double) - Dwell time, default to 10e-6s
%   TR (double) - Repetition time, default to 50e-3s
%   RFA - Round flip angle , default ernstAngle at 7T
%
% Outputs:
%   rf, grad, localSAR - Resulting RF & grad pulse based on the input parameters, and
%   local SAR according to the vendor-provided 'SARDataUser'
%
% Example:
%   result = design_KT_rf(); % Uses default parameters
%   result = design_KT_rf(5,2.88,10e-6,50e-3,ernstAngle(50e-3)); % Uses custom parameters. Pay
%   attention the RF_duration can not be changed in the version 1.0.

% Set default values
if nargin < 1
    ID = [1];
end
if nargin < 2
    RF_duration = 2.88e-3;
end
if nargin < 3
    dt = 10e-6;
end
if nargin < 4
    TR = 50e-3;
end
if nargin < 5
    RFA = 13;
end

load KT_X0 xmin
if length(ID) == 1 %%%subject specific (need to be quick), premalloc optimization & pre calculation used
    [theta,vx0] = mainKTopt(ID(1),RF_duration, xmin, dt, TR, RFA); %%% step1: reduced problem solve
else %%% universal pulse training, no acceleration here
end
kp = reshape(theta,[3,24]);
nchs = 8;
nKT = size(kp,2);
maxamp = 50e-3;
maxsr = 160;
kpback = [kp zeros(3,1)];
dkval = diff(kpback,1,2);
sgrad = {};

nblips = size(dkval,2);
for ind= 1: nblips     
    dkmax = max(abs(dkval(:,ind)));
    %blipmax = design_toptgrad1D(0,0,dkmax,maxamp,maxsr,dt);
    blipmax = design_grad_trapz(dkmax,maxamp,maxsr,dt);
    blips = [blipmax;blipmax;blipmax];
    sf = diag(dkval(:,ind)./ dkmax);
    sgrad{ind} = sf*blips;  
end

npts = length(sgrad);
subrfLen = zeros(1,npts);

bw= 2*1050; %ideal bandwidth
dt0 = 1/bw;
ntimepts = round(dt0/dt);

npts0= 3;
for ind=1:npts0:npts,
    g0=[];
    for jdx = 1:npts0
        g0= [g0 sgrad{ind+jdx-1}];
    end

    nt = round((ntimepts - size(g0,2))/npts0); 
    nt = 5;
    if nt>0
        subrfLen(ind:ind+npts0-1) = nt; else
        subrfLen(ind:ind+npts0-1) = 1;
    end
end

phasetrack = subrfLen;
for ind= 1: length(sgrad)-1
    phasetrack(ind+1) = phasetrack(ind)+ size(sgrad{ind},2)+ subrfLen(ind+1);
end


wt11 = reshape(vx0,[],nchs);
wt11 = wt11.';
srf = cell(size(subrfLen));
for ind=1:length(subrfLen)
    nt = subrfLen(ind);
    srf{ind} = 1/nt.* wt11(:,ind) * ones(1,nt);
end

mygrad = [];
myrf = [];
nchs = size(srf{1},1);
npts = length(srf);
for ind=1:npts
    mygrad = [mygrad zeros(3,size(srf{ind},2)),sgrad{ind}];
    myrf = [myrf srf{ind} complex(zeros(nchs,size(sgrad{ind},2)))];
end

grad = [zeros(3,1),mygrad];
rf = [zeros(8,1),myrf];

rf = RFA*rf*1e6; grad = grad*1e3; %%% uV->V, T->mT
localSAR = localSARcom(rf,TR,dt);
end