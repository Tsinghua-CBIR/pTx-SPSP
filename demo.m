clear
close all
clc
addpath('data/');
addpath('src/');


trainID = [1]; %%% serial number of calibration
%%% Tips: put more than one number here for universal pulse training;
%%% For example, trainID = 1:10
%%% Remember to left at least one subject for testing in reality
%%% Here in the demo, as we only have 2 subjects,
%%% Using trainID = [1,5] and testID = 1 just for example

RF_duration = 2.88e-3; %%% duration of the RF pulse, can not be changed the this version
dt = 10e-6; %%% dwell time, s, can not be changed in the demo
TR = 50e-3; %%% repetition time, s, used in SAR-constraint.
RFA = round(ernstAngle(TR)); %%% round: the vendor-provided FA is integer

% design with kT point parameterization 
% [rf,grad,localSAR] = design_pTxSPSP_RF(trainID,RF_duration,dt,TR,RFA,'KT');
% design with SPINS parameterization 
[rf,grad,localSAR] = design_pTxSPSP_RF(trainID,RF_duration,dt,TR,RFA,'SPINS');

showPulse(rf,grad,RF_duration,dt);
offset = 100;%%% To check the robustness to off-resonance
%%% As in the paper, you may change offset to 0 and ±100Hz and even ±200Hz
%%% pTx-SPSP water-excitation pulse is  quite robust to off-resonances 

testID = 1;
showPerform(testID,rf,grad,dt,offset,RFA);

%%% To see the frequency responce if you want
%%% warning: may take some time
% showFrequencyRes(testID,rf,grad,dt,RFA);

%%% Output to ini file
outputfile(rf,grad);


