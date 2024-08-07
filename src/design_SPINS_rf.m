function [rf,grad,localSAR] = design_SPINS_rf(ID, RF_duration, dt, TR, RFA)
% rf_design_function - MATLAB function for designing a non-selective pTx-SPSP SPINS RF pulse
%
% Inputs:
%   ID (int) - serial number of calibration, set here so the function
%   could be easily extended to UP training.
%   RF_duration (double) - Duration of the RF pulse, default to 2.88 ms
%
% Outputs:
%   rf, grad, localSAR - Resulting RF & grad pulse based on the input parameters, and
%   local SAR according to the vendor-provided 'SARDataUser'
%
% Example:
%   result = design_SPINS_rf(); % Uses default parameters
%   result = design_SPINS_rf(5,2.88,10e-6,50e-3,ernstAngle(50e-3)); % Uses custom parameters. Pay
%   attention the RF_dutation can not be changed in the version 1.0.

% Set default values
if nargin < 1
    ID = 5;
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
nKT = 12;%%% accelarating facotr for SPINS
nchs = 8;%%% number of channel, considering commonly used Nova 8-ch transmit
%%% if you want to use another tranmit coil, just remember to change "nchs"
%%% in all the files.
load SPINS_X0.mat %%% optimal startpoint
[theta,vx0] = mainSPINSopt(ID,RF_duration, X0, dt, TR, RFA); %%% step1: reduced problem solve

rf = vx0;
% vx0 = repelem(vx0,nKT,1);
% Hp = [theta;real(vx0);imag(vx0)];
% Wp = finetuneSPINS(5, RF_duration, dt, Hp);
% theta = Wp(1:5);
% rf = Wp(6:2309)+1i*Wp(2310:4613);


rf = reshape(rf,[],nchs);
rf = rf.';
gamma = 2.675e8;
T = RF_duration;
len = round(T/dt);
kp = zeros(3,len);
Kpp = theta(1);
Ktt = theta(2)/6;
Kphph = theta(3)/30;
alpha = theta(4)/2;
beta = theta(5)/40;
%%% Weighting, as the parameters of SPINS have been normalized

for i = 1:len
    ttt = i*dt;
    Kp = Kpp/(1+exp(alpha*(ttt/T-beta)));
    Ktheta = Ktt*pi*1000*ttt;
    Kphi = Kphph*pi*1000*ttt;
    Kp_temp = [Kp*sin(Ktheta)*cos(Kphi);Kp*sin(Ktheta)*sin(Kphi);Kp*cos(Ktheta)];
    kp(:,i) = Kp_temp;
end

dkval = diff(kp,1,2);
grad = dkval/gamma/dt;
grad = [zeros(3,2),grad];
rf = [zeros(8,1),rf];


localSAR = localSARcom(rf*RFA*1e6,TR,dt);;
end